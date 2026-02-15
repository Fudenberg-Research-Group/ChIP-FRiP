#!/usr/bin/env python
from optparse import OptionParser
import pandas as pd
import pysam
import deeptools.countReadsPerBin as crpb


def main():
    usage = "usage: %prog [options] <bam_file1> <bed_file>"
    parser = OptionParser(usage)

    parser.add_option(
        "--nproc",
        dest="nproc",
        type="int",
        default=1,
        help="How many processes to use for calculation. Default: nproc=1",
    )

    (options, args) = parser.parse_args()

    if len(args) != 1:
        parser.error(
            "Must provide a pathway to metadata table (txt file) that includes four columns [prefix of output file] [output-directory] [BAM file pathway]  and [BED file pathway]."
        )
    else:
        metadata = pd.read_table(args[0])

    frips = []
    for i, row in metadata.iterrows():
        bam = row["BAM"]
        bed = row["BED"]

        reads_counter = crpb.CountReadsPerBin(
            [bam], bedFile=bed, numberOfProcessors=options.nproc
        )
        reads_at_peaks = reads_counter.run()
        total_reads_at_peaks = reads_at_peaks.sum(axis=0)

        alignment = pysam.AlignmentFile(bam)
        frip = float(total_reads_at_peaks[0]) / alignment.mapped
        frips.append(frip)

    metadata["FRiP"] = frips
    metadata.to_csv(args[0], sep="\t", index=False)


if __name__ == "__main__":
    main()
