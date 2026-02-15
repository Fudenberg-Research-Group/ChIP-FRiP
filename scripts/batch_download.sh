#!/bin/bash

while read -r accession; do
    echo "Downloading $accession..."
    prefetch "$accession"
    fasterq-dump "$accession"
    echo "$accession downloaded."
done < accessions.txt
