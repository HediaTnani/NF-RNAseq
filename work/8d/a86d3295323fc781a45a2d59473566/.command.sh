#!/bin/bash -ue
mkdir index_dir

STAR --runThreadN 2 \
  --runMode genomeGenerate \
  --genomeDir index_dir/ \
  --genomeFastaFiles genome.fa \
  --genomeSAindexNbases 10 \
  --sjdbGTFfile annotation.gtf
