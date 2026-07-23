#!/usr/bin/env python3
import sys
import os
import subprocess
import tempfile
import argparse
import pybedtools

def main():
    parser = argparse.ArgumentParser(
        description="Pipeline2: HMM-guided genomic search via MMseqs2 + dummer",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Examples:\n"
               "  python3 pipeline2.py profile.hmm msa.msa genome.fa 8\n"
               "  python3 pipeline2.py profile.hmm msa.msa genome.fa 8 --target-db-pad /path/to/targetDB_pad --query-db /path/to/queryDB\n"
    )
    parser.add_argument("hmm_file", help="HMM profile file")
    parser.add_argument("msa_file", help="MSA file (.msa for Stockholm, otherwise used as-is for query DB). Ignored when --query-db is provided.")
    parser.add_argument("fa_file", help="Genome FASTA file")
    parser.add_argument("cpus", help="Number of CPUs")
    parser.add_argument("--target-db-pad", dest="target_db_pad", default=None,
                        help="Path to an existing padded target DB (skips createdb + makepaddedseqdb)")
    parser.add_argument("--query-db", dest="query_db", default=None,
                        help="Path to an existing query profile DB (skips convertmsa + msa2profile)")
    parser.add_argument("--skip-dummer", action="store_true",
                        help="Skip running dummer (only generate debug.fa)")
    parser.add_argument("--output-fa", dest="output_fa", default=None,
                        help="Save debug FASTA to this path (persistent copy)")
    parser.add_argument("--max", action="store_true",
                        help="Max sensitivity: disable heuristic windowing (pad=full contig)")
    parser.add_argument("--no-seeds", action="store_true",
                        help="Omit seed annotations from FASTA headers (dummer runs without seed gating)")
    parser.add_argument("--dummer-bin", dest="dummer_bin", default=None,
                        help="Path to dummer binary (overrides default dummer)")
    parser.add_argument("--prefilter-mode", type=int, default=None, choices=[3],
                        help="MMseqs2 --prefilter-mode 3 (GPU combined ungapped+gapped). "
                             "Seeds use only endpoint positions (point seeds).")

    args = parser.parse_args()

    hmm_file = args.hmm_file
    msa_file = args.msa_file
    fa_file = args.fa_file
    cpus = args.cpus

    script_dir = os.path.dirname(os.path.realpath(__file__))
    dummer_exec = args.dummer_bin or os.path.join(script_dir, "../cmake-build-release/dummer")

    # ---------------------------------------------------------
    # 1. Parse sequences and HMMs
    # ---------------------------------------------------------
    hmm_lens = {}
    with open(hmm_file, 'r') as f:
        curr_acc, curr_name = None, None
        for line in f:
            if line.startswith("NAME"):
                curr_name = line.split()[1].strip()
            elif line.startswith("ACC"):
                curr_acc = line.split()[1].strip()
            elif line.startswith("LENG"):
                length = int(line.split()[1].strip())
                if curr_acc: hmm_lens[curr_acc] = length
                if curr_name: hmm_lens[curr_name] = length
                curr_acc, curr_name = None, None

    dna_lens = {}
    dna_seqs = {}
    with open(fa_file, 'r') as f:
        curr_id = None
        curr_seq = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if curr_id:
                    dna_lens[curr_id] = len("".join(curr_seq))
                    dna_seqs[curr_id] = "".join(curr_seq)
                curr_id = line[1:].split()[0]
                curr_seq = []
            elif curr_id:
                curr_seq.append(line)
        if curr_id:
            dna_lens[curr_id] = len("".join(curr_seq))
            dna_seqs[curr_id] = "".join(curr_seq)

    # total search space: all genome lengths, doubled for both strands
    tot_seq_len = sum(dna_lens.values()) * 2

    # ---------------------------------------------------------
    # --max mode: skip mmseqs2, run dummer directly on all profiles x contigs
    # ---------------------------------------------------------
    if args.max:
        merged_fa_path = os.path.join(tempfile.gettempdir(), "dummer_max.fa")
        trans = str.maketrans("ACGTacgt", "TGCAtgca")

        with open(merged_fa_path, "w") as fout:
            for chrom, seq in dna_seqs.items():
                L = len(seq)
                rc_seq = seq.translate(trans)[::-1]
                for prof in sorted(hmm_lens.keys()):
                    fout.write(f">{chrom}/1-{L} length={L} profile={prof} plus_strand\n")
                    fout.write(f"{seq}\n")
                    fout.write(f">{chrom}/1-{L} length={L} profile={prof} minus_strand_revcomp\n")
                    fout.write(f"{rc_seq}\n")

        print(f"# Max-mode FASTA written to: {merged_fa_path}")

        if not args.skip_dummer:
            try:
                subprocess.run(
                    [dummer_exec, hmm_file, merged_fa_path, '-T 8', '-W 0.1', '-N', str(tot_seq_len)],
                    env=os.environ.copy(), check=True,
                )
            except subprocess.CalledProcessError as e:
                print(f"Error: dummer encountered an issue (Exit status: {e.returncode})")
                sys.exit(1)

        if args.output_fa:
            import shutil
            shutil.copy2(merged_fa_path, args.output_fa)
            print(f"# Debug FASTA saved to: {args.output_fa}")
        else:
            os.remove(merged_fa_path)

        return

    with tempfile.TemporaryDirectory(prefix="mmseqs_tmp_", delete=True) as tmpdir:
        print(f"# Temporary directory is: {tmpdir}")

        # ---------------------------------------------------------
        # 3. MMseqs Protein-Protein Search
        # ---------------------------------------------------------
        db_dir = os.path.join(tmpdir, "db")
        os.makedirs(db_dir)

        ali_file = os.path.join(db_dir, f"result.ali")
        tmp_file = os.path.join(db_dir, f"tmp")

        if args.target_db_pad:
            target_db_pad = args.target_db_pad
            print(f"# Using existing target_db_pad: {target_db_pad}")
        else:
            prot_fa_path = os.path.join(tmpdir, "translated_6frame.pfa")
            subprocess.run([
                "seqkit", "translate", "-f", "6", "-F",
                "--threads", cpus,
                "-o", prot_fa_path, fa_file
            ], check=True)

            target_db_pad = os.path.join(db_dir, "targetDB_pad")
            subprocess.run(["mmseqs", "createdb", prot_fa_path, target_db_pad, "--gpu", "1"], check=True, stdout=subprocess.DEVNULL)

        if args.query_db:
            query_db = args.query_db
            print(f"# Using existing query_db: {query_db}")
        else:
            query_db = os.path.join(db_dir, "queryDB")
            msa_db = os.path.join(db_dir, "msa_db")
            subprocess.run(["mmseqs", "convertmsa", msa_file, msa_db, "--identifier-field", "0"], check=True, stdout=subprocess.DEVNULL)
            subprocess.run(["mmseqs", "msa2profile", msa_db, query_db], check=True, stdout=subprocess.DEVNULL)

        mmseqs_cmd = [
            "mmseqs", "search", query_db, target_db_pad, ali_file, tmpdir,
            "--gpu", "1",
            "--threads", cpus,
            #"-e", "10000",
            "-e", "1000",
            #"--min-ungapped-score", "0",
            #"--num-iterations", "3",
            "--alignment-mode", "2",
        ]
        if args.prefilter_mode is not None:
            mmseqs_cmd.extend(["--prefilter-mode", str(args.prefilter_mode)])
            mmseqs_cmd[mmseqs_cmd.index("--alignment-mode") + 1] = "1"

        # mmseqs_cmd = [
        #     "mmseqs", "search", query_db, target_db_pad, ali_file, tmpdir,
        #     "--threads", cpus,
        #     #"-e", "10000",
        #     "-e", "10000",
        #     "-s", "10.5",
        #     "--alignment-mode", "1",
        # ]

        subprocess.run(mmseqs_cmd, check=True)
        subprocess.run(["mmseqs", "convertalis", query_db, target_db_pad, ali_file, tmp_file], check=True, stdout=subprocess.DEVNULL)

        # ---------------------------------------------------------
        # 4. Map Amino Acid Hits -> Genomic DNA Windows
        # ---------------------------------------------------------
        #
        # Coordinate conventions:
        #
        #   mmseqs convertalis output (BLAST-tab fields):
        #     fields[6]=qstart, fields[7]=qend  — 1-indexed inclusive profile positions
        #     fields[8]=tstart, fields[9]=tend  — 1-indexed inclusive target protein positions
        #
        #   Six-frame translation:
        #     Forward frame N (N=1,2,3): protein[1] ↔ genome nucleotide [N-1] (0-indexed)
        #       protein position p ↔ genome nucleotide (N-1) + 3*(p-1)   (1st nt of codon)
        #       protein position p ↔ genome nucleotide (N-1) + 3*(p-1)+2 (3rd nt of codon)
        #     Reverse frame N (N=1,2,3): translation runs over reverse-complement of genome.
        #       RC position i (0-indexed) ↔ original genome position L-1-i
        #       protein position p in frame N ↔ 1st nt of codon at original position:
        #         L - N - 3*p + 3   (0-indexed, inclusive)
        #
        #   Seed format "seed=prof_s,dna_s,prof_e,dna_e" (semicolons between seeds):
        #     All values are 0-indexed inclusive.
        #     prof_s, prof_e — positions in the HMM profile (query)
        #     dna_s, dna_e   — positions within the EXTRACTED FASTA window
        #       Forward strand: dna position 0 = genome position (start-1), left-to-right
        #       Reverse strand: bedtools getfasta -s produces reverse-complemented FASTA;
        #         dna position 0 = genome position (end-1), right-to-left
        #     dna_s = first  nucleotide of the codon at prof_s
        #     dna_e = first  nucleotide of the codon at prof_e
        #
        #   Genomic overlap check (0-indexed half-open genome intervals):
        #     gen_s  — 0-indexed inclusive start of alignment in genome
        #     gen_e  — 0-indexed exclusive end   of alignment in genome
        #     Window — [start-1, end)  (0-indexed, half-open)
        #     Overlap: gen_e > (start-1) AND gen_s < end
        #
        hits_by_window = {}

        def parse_mmseqs_to_intervals(filepath):
            for line in open(filepath):
                if line.startswith("#"):
                    continue

                fields = line.split()
                if len(fields) < 12:
                    continue

                query_acc = fields[0]
                target_full = fields[1]
                q_start = int(fields[6])
                q_end   = int(fields[7])
                t_start = int(fields[8])
                t_end = int(fields[9])
                hmm_len = hmm_lens.get(query_acc, 0)
                if args.prefilter_mode == 3:
                    p_pos = max(1, t_end - (hmm_len // 2))
                else:
                    p_pos = max(1, (t_start + t_end) // 2)
                e_value = fields[10]
                bitscore = fields[11]
                
                *target_parts, frame_str = target_full.rsplit('_', 1)
                target_base = '_'.join(target_parts)
                frame_val = int(frame_str.split('=')[1])
                strand, frame = ('F', frame_val) if frame_val > 0 else ('R', abs(frame_val))

                if args.prefilter_mode == 3:
                    hits_by_window.setdefault((target_base, query_acc, strand), []).append({
                        'q_start': q_end,
                        'q_end':   q_end,
                        't_start': t_end,
                        't_end':   t_end,
                        'frame':   frame,
                        'bitscore': float(bitscore),
                    })
                else:
                    hits_by_window.setdefault((target_base, query_acc, strand), []).append({
                        'q_start': q_start,
                        'q_end':   q_end,
                        't_start': t_start,
                        't_end':   t_end,
                        'frame':   frame,
                        'bitscore': float(bitscore),
                    })

                L = dna_lens.get(target_base, 0)
                if L == 0: 
                    continue
                pad = 3 * hmm_lens.get(query_acc, 0)

                base_pos = (frame - 1) + 3 * (p_pos - 1)

                if strand == 'F':
                    start = max(0, base_pos - pad)
                    end   = min(L, base_pos + 1 + pad)
                    fake_chrom = f"{target_base}|{query_acc}|+"
                    yield pybedtools.Interval(fake_chrom, start, end, query_acc, bitscore, '+')
                else:
                    start = max(0, L - base_pos - 1 - pad)
                    end   = min(L, L - base_pos + pad)
                    fake_chrom = f"{target_base}|{query_acc}|-"
                    yield pybedtools.Interval(fake_chrom, start, end, query_acc, bitscore, '-')

        def unpack_intervals(feature):
            real_chrom, query_acc, strand = feature.chrom.split('|')
            return pybedtools.Interval(real_chrom, feature.start, feature.end, query_acc, ".", strand)

        merged_bed = pybedtools.BedTool(parse_mmseqs_to_intervals(tmp_file)) \
            .sort() \
            .merge() \
            .each(unpack_intervals)

        # ---------------------------------------------------------
        # 5. Extract Final Genomic FASTA (Bedtools)
        # ---------------------------------------------------------
        merged_bed_path = os.path.join(tmpdir, "merged.bed")
        merged_bed.saveas(merged_bed_path)
        
        raw_fa_path = os.path.join(tmpdir, "raw_ext.fa")
        merged_fa_path = os.path.join(tmpdir, f"debug.fa")

        subprocess.run(["bedtools", "getfasta", "-fi", fa_file, "-bed", merged_bed_path, "-s", "-name+", "-fo", raw_fa_path], check=True)

        with open(raw_fa_path, 'r') as fin, open(merged_fa_path, 'w') as fout:
            for line in fin:
                if line.startswith(">"):
                    header_data = line[1:].strip()
                    query, coord_part = header_data.split("::") if "::" in header_data else ("UNKNOWN", header_data)
                    chrom, rest = coord_part.split(":")
                    pos, strand_part = rest.split("(")
                    start_str, end_str = pos.split("-")
                    strand_sign = strand_part[0]
                    
                    start, end = int(start_str), int(end_str)
                    strand_label = "plus_strand" if strand_sign == '+' else "minus_strand_revcomp"

                    seed_str = ""
                    if not args.no_seeds:
                        strand_key = 'F' if strand_sign == '+' else 'R'
                        hits = hits_by_window.get((chrom, query, strand_key), [])
                        if hits:
                            seeds = []
                            for h in hits:
                                frame   = h['frame']
                                t_start = h['t_start']
                                t_end   = h['t_end']
                                prof_s  = h['q_start'] - 1
                                prof_e  = h['q_end']   - 1
                                L = dna_lens.get(chrom, 0)
                                if L > 0:
                                    win_start = start - 1
                                    win_end   = end
                                    if strand_sign == '+':
                                        gen_s = (frame - 1) + 3 * (t_start - 1)
                                        gen_e = (frame - 1) + 3 * (t_end   - 1) + 3
                                        if gen_e <= win_start or gen_s >= win_end:
                                            continue
                                        dna_s = gen_s - win_start
                                        dna_e = gen_e - 3 - win_start
                                    else:
                                        gen_s = L - frame - 3 * t_end + 1
                                        gen_e = L - frame - 3 * t_start + 4
                                        if gen_e <= win_start or gen_s >= win_end:
                                            continue
                                        dna_s = (end - 1) - (L - frame - 3 * t_start + 3)
                                        dna_e = (end - 1) - (L - frame - 3 * t_end   + 3)
                                    win_size = end - start + 1
                                    dna_s = max(0, min(dna_s, win_size - 1))
                                    dna_e = max(0, min(dna_e, win_size - 1))
                                    seeds.append(f"{prof_s},{dna_s},{prof_e},{dna_e}")
                            seed_str = f" seed={';'.join(seeds)}"

                    fout.write(f">{chrom}/{start+1}-{end} length={dna_lens.get(chrom, 0)} profile={query} {strand_label}{seed_str}\n")
                else:
                    fout.write(line)

        if args.output_fa:
            import shutil
            shutil.copy2(merged_fa_path, args.output_fa)
            print(f"# Debug FASTA saved to: {args.output_fa}")

        # ---------------------------------------------------------
        # 6. Run Dummer
        # ---------------------------------------------------------
        if not args.skip_dummer:
            custom_env = os.environ.copy()
            #custom_env["ASAN_OPTIONS"] = "detect_container_overflow=1:strict_memcmp=1"
            
            try:
                subprocess.run([dummer_exec, hmm_file, merged_fa_path, '-T 8', '-W 0.5' if not args.max else '-W 10', '-N', str(tot_seq_len)], env=custom_env, check=True)
            except subprocess.CalledProcessError as e:
                print(f"Error: dummer encountered an issue (Exit status: {e.returncode})")
                sys.exit(1)

if __name__ == "__main__":
    main()