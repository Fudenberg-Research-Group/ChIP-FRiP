import pysam
import pandas as pd
import numpy as np
import os
import subprocess
import json
import bioframe as bf
from multiprocessing import Pool
from functools import partial

pysam.set_verbosity(0)


####### Helper function #######
def fetch_metadata(accession):
    # Use subprocess to run ffq and capture the output
    result = subprocess.run(["ffq", accession], capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error fetching data for {accession}: {result.stderr}")
        return None
    return json.loads(result.stdout)

def _reads_in_chrom_peak(bam, peaks, chrom):
        alignments = bf.read_alignments(bam, chrom)
        alignments = alignments.rename(columns={"POS": "start"})
        alignments["chrom"] = chrom
        alignments["end"] = alignments["start"] + alignments["SEQ"].apply(len)

        try:
             result = bf.count_overlaps(alignments, peaks).iloc[:,-1].tolist()
        except TypeError:
             result = np.array([])
        
        return result

def count_reads_in_peak(bam, peaks, nproc=2):
    chromosomes = peaks['chrom'].unique()

    if nproc > len(chromosomes):
        nproc = len(chromosomes)
    
    with Pool(processes=nproc) as pool:
        results = pool.map(partial(_reads_in_chrom_peak, bam, peaks), chromosomes)
        results = np.concatenate(results)

    return np.sum(results > 0)

def calculate_frip(bam, peaks_table, blacklist=None, nproc=40):
    alignment = pysam.AlignmentFile(bam)
    read = next(alignment.fetch())
    read_length = read.query_length

    total_reads = alignment.mapped

    if blacklist is not None:
        # Merging ensures we don't double-count reads in overlapping blacklist regions.
        blacklist_merged = bf.merge(blacklist)

        reads_in_blacklist = 0
        for _, row in blacklist_merged.iterrows():
            # pysam.count is efficient for fetching specific regions
            reads_in_blacklist += alignment.count(
                contig=row['chrom'], 
                start=row['start'], 
                stop=row['end']
            )
        total_reads = total_reads - reads_in_blacklist
        if total_reads <= 0:
            raise ValueError("No reads remain after blacklist filtering")

    total_reads_at_peaks = count_reads_in_peak(bam, peaks_table, nproc=nproc)
    frip = float(total_reads_at_peaks) / total_reads # floa can avoid automatical round

    # Calculate total flank region width within which a read can still be considered overlapping with the peak
    distance_between_peaks = bf.closest(peaks_table, None)
    distance_between_peaks[distance_between_peaks['distance'].isna()] = read_length # NA is because there is no peaks within the same scaffold or chromosome, so we set the distance to read_length for the below process
    distance_between_peaks = distance_between_peaks['distance'].to_numpy()
    distance_between_peaks[distance_between_peaks >= read_length] = read_length - 1
    flank_regions_for_reads_overlap_peaks = distance_between_peaks.sum()

    return frip, total_reads_at_peaks, total_reads, flank_regions_for_reads_overlap_peaks

def create_frip_table_from_bed(
    samples_metadata,
    path_to_bed,
    path_to_data,
    genome_size,
    species,
    nproc,
    peak_protein_srun="",
    customized_metadata=False,
    blacklist=None
):

    if customized_metadata:
        samples_metadata = samples_metadata[samples_metadata["BAM"] != ""]
        bams = samples_metadata["BAM"].to_list()
    else:
        sruns = samples_metadata["SRUN"].to_list()
        if os.path.exists(f"{path_to_data}/{sruns[0]}/{sruns[0]}.q30.{species}.sort.bam"):
            suffix = f"{species}.sort"
        else:
            suffix = "dedup"
        bams = [f"{path_to_data}/{sample}/{sample}.q30.{suffix}.bam" for sample in sruns]

    peaks_table = bf.read_table(path_to_bed, schema="bed").iloc[:, :3]
    if blacklist is not None:
        peaks_table = bf.subtract(peaks_table, blacklist)
        # Ensure no invalid intervals remain (when start >= end) after subtraction
        peaks_table = peaks_table[peaks_table['end'] > peaks_table['start']].copy()
        if peaks_table.empty:
            raise ValueError("No peaks remain after blacklist filtering")
    num_peaks = [peaks_table.shape[0]] * len(bams)
    # count basepairs within peaks
    total_bp_in_peaks = (peaks_table["end"] - peaks_table["start"]).sum() + len(peaks_table) # include the basepair at start

    total_reads = []
    frip_enrich = []
    samples_frips = []
    for idx, bam in enumerate(bams):
        print(f"### bams: {idx + 1}/{len(bams)}")
        result = calculate_frip(bam, peaks_table, blacklist=blacklist, nproc=nproc)
        samples_frips.append(result[0])
        total_reads.append(result[2])
        flank_regions_for_reads_overlap_peaks = result[3]
        expected_prob_read_in_peak = (total_bp_in_peaks + flank_regions_for_reads_overlap_peaks) / genome_size
        frip_enrich.append(result[0] / expected_prob_read_in_peak) # This is the ratio of observed reads in peaks / expected reads in peaks

    frip_df = pd.DataFrame({"FRiP": samples_frips})
    extra_df = pd.DataFrame(
        {
            "FRiP enrichment": frip_enrich,
            "#Peaks": num_peaks,
            "Total #basepairs in peaks": total_bp_in_peaks,
            "Total #reads": total_reads,
        }
    )
    frip_df = pd.concat([frip_df, samples_metadata, extra_df], axis=1)
    
    # replace GSM_accession column with peaks-SRA
    frip_df["GSM_accession"] = [peak_protein_srun] * len(frip_df)
    frip_df = frip_df.rename(columns={"GSM_accession": "Peak-SRUN"})

    return frip_df
