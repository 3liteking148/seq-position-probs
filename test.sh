#make -j2 && time bash bin/pipeline.sh AMRFinder_filtered.hmm originals300.fasta
make -j2 && time bin/dummer -v MET-test.hmm met-optimized.fa
#make -j2 && time bin/dummer -v MET-test.hmm MET-target.fa
