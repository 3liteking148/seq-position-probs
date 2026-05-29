#!/usr/bin/env python3
import sys
import os
import subprocess
import tempfile
import pybedtools

COMPLEMENT_TABLE = str.maketrans("ATCGatcgNn", "TAGCtagcNn")

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def translate(seq):
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    }
    protein = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3].upper()
        protein.append(table.get(codon, 'X'))
    return "".join(protein)

def six_frame_translate(input_fasta, output_fasta):
    with open(input_fasta, 'r') as f, open(output_fasta, 'w') as out:
        header = None
        seq = []
        
        def process_seq(h, s):
            full_seq = "".join(s)
            rc_seq = reverse_complement(full_seq)
            
            # Forward frames
            for frame in range(3):
                prot = translate(full_seq[frame:])
                out.write(f"{h}_F{frame+1}\n{prot}\n")
            
            # Reverse frames
            for frame in range(3):
                prot = translate(rc_seq[frame:])
                out.write(f"{h}_R{frame+1}\n{prot}\n")

        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if header:
                    process_seq(header, seq)
                header = line
                seq = []
            else:
                seq.append(line)
        if header:
            process_seq(header, seq)

def main():
    usage = (
        "Usage: python3 run_pipeline.py <hmm_file> <msa_file> <fa_file> <cpus>\n"
        "Example: python3 run_pipeline.py profile.hmm msa.msa genome.fa 8"
    )

    if len(sys.argv) != 5:
        print(usage)
        sys.exit(1)

    hmm_file = sys.argv[1]
    msa_file = sys.argv[2]
    fa_file = sys.argv[3]
    cpus = sys.argv[4]

    script_dir = os.path.dirname(os.path.realpath(__file__))
    dummer_exec = os.path.join(script_dir, "dummer")

    # ---------------------------------------------------------
    # 1. Parse sequences and HMMs
    # ---------------------------------------------------------
    print("Parsing sequences and HMMs")
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
    with open(fa_file, 'r') as f:
        curr_id = None
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                curr_id = line[1:].split()[0]
                dna_lens[curr_id] = 0
            elif curr_id:
                dna_lens[curr_id] += len("".join(line.split()))

    with tempfile.TemporaryDirectory(prefix="mmseqs_tmp_", delete=False) as tmpdir:
        print(f"Temporary directory is: {tmpdir}")

        print("Generating 6-frame protein translation...")
        prot_fa_path = os.path.join(tmpdir, "translated_6frame.pfa")
        six_frame_translate(fa_file, prot_fa_path)

        # ---------------------------------------------------------
        # 3. MMseqs Protein-Protein Search
        # ---------------------------------------------------------
        print("Setting up MMseqs protein databases...")
        db_dir = os.path.join(tmpdir, "db")
        os.makedirs(db_dir)

        target_db = os.path.join(db_dir, "targetDB")
        target_db_pad = os.path.join(db_dir, "targetDB_pad")
        query_db = os.path.join(db_dir, "queryDB")
        ali_file = os.path.join(db_dir, f"result.ali")
        tmp_file = os.path.join(db_dir, f"tmp")

        subprocess.run(["mmseqs", "createdb", prot_fa_path, target_db], check=True, stdout=subprocess.DEVNULL)
        subprocess.run(["mmseqs", "makepaddedseqdb", target_db, target_db_pad], check=True, stdout=subprocess.DEVNULL)

        if True: # if original msa file
            msa_db = os.path.join(db_dir, "msa_db")
            subprocess.run(["mmseqs", "convertmsa", msa_file, msa_db, "--identifier-field", "0"], check=True, stdout=subprocess.DEVNULL)
            subprocess.run(["mmseqs", "msa2profile", msa_db, query_db], check=True, stdout=subprocess.DEVNULL)
        else:
            # TODO
            pass

        print(f"Running MMseqs-GPU")
        mmseqs_cmd = [
            "mmseqs", "search", query_db, target_db_pad, ali_file, tmpdir,
            "--max-seqs", "100000", 
            "-e", "1000",
            "--gpu", "1",
            "--spaced-kmer-mode", "1",
            "--prefilter-mode", "1",
            "--alignment-mode", "1", # Pure GPU scoring mode (TODO: check if placebo or real)
            "--k-score", "60" # from nail but lower
        ]
        subprocess.run(mmseqs_cmd, check=True, stdout=subprocess.DEVNULL)
        subprocess.run(["mmseqs", "convertalis", query_db, target_db_pad, ali_file, tmp_file], check=True, stdout=subprocess.DEVNULL)

        # ---------------------------------------------------------
        # 4. Map Amino Acid Hits -> Genomic DNA Windows
        # ---------------------------------------------------------
        print("Merging overlapping alignment windows...")

        # todo: double check ts
        def parse_mmseqs_to_intervals(filepath):
            for line in open(filepath):
                if line.startswith("#"):
                    continue

                fields = line.split()
                if len(fields) < 12:
                    continue

                query_acc = fields[0]
                target_full = fields[1]
                t_start = int(fields[8]) # 0 in align-mode 1
                t_end = int(fields[9])
                hmm_len = hmm_lens.get(query_acc, 0)
                p_pos = max(1, t_end - (hmm_len // 2))

                bitscore = fields[11]
                
                *target_parts, strand_frame = target_full.rsplit('_', 1)
                target_base = '_'.join(target_parts)
                strand, frame = strand_frame[0], int(strand_frame[1])

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
        print("Extracting final padded genomic sequences...")
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
                    fout.write(f">{chrom}/{start+1}-{end} length={dna_lens.get(chrom, 0)} profile={query} {strand_label}\n")
                else:
                    fout.write(line)

        # ---------------------------------------------------------
        # 6. Run Dummer
        # ---------------------------------------------------------
        print("Executing dummer analysis...")
        custom_env = os.environ.copy()
        #custom_env["ASAN_OPTIONS"] = "detect_container_overflow=1:strict_memcmp=1"
        
        try:
            subprocess.run(["time", dummer_exec, hmm_file, merged_fa_path], env=custom_env, check=True)
            print("Pipeline completed successfully.")
        except subprocess.CalledProcessError as e:
            print(f"Error: dummer encountered an issue (Exit status: {e.returncode})")
            sys.exit(1)

if __name__ == "__main__":
    main()