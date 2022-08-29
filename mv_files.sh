find data/ -name "*fastq.gz" -print0 | xargs -0 -I {} mv {} data/
