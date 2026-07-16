#time bash bin/pipeline.sh AMRFinder_filtered.hmm originals300.fasta
#time bash bin/pipeline.sh AMRFinder_filtered.hmm originals300.fasta
#\\time cmake-build-release/dummer -v MET-test.hmm met-optimized.fa
source ~/.venv/bin/activate
#time python3 bin/pipeline2.py MET-test.hmm MET.msa MET-target.fa 1
time python3 bin/pipeline2.py MET-test.hmm MET.msa MET-target.fa 1 --no-seeds --prefilter-mode 3
#time bash bin/pipeline.sh MET-test.hmm MET-target.fa
#time bin/dummer -v /home/harvsftw/Documents/bioinfo/BATH-paper/transmark/test-001.AA.hmm /mnt/tmp/test.txt