#!/bin/bash

if [ ! -d "../data/fastq" ]; then
    mkdir -p ../data/fastq
fi

while read -r accession; do
    echo "Downloading $accession..."
    prefetch "$accession"
    fasterq-dump "$accession" -O "../data/fastq/"
    echo "$accession downloaded."
done < ../data/accessions.txt
