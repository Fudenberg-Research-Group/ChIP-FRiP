import pysam
import pandas as pd
import os
import subprocess
import json
import bioframe as bf
import deeptools.countReadsPerBin as crpb

pysam.set_verbosity(0)


# Helper function
def fetch_metadata(accession):
    # Use subprocess to run ffq and capture the output
    result = subprocess.run(["ffq", accession], capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error fetching data for {accession}: {result.stderr}")
        return None
    return json.loads(result.stdout)


def calculate_frip(bam, bed, nproc=40):
    alignment = pysam.AlignmentFile(bam)
    read = next(alignment.fetch())
    read_length = read.query_length
    peaks_table = bf.read_table(bed, schema='bed').iloc[:,:3]

    distance_between_peaks = bf.closest(peaks_table, None)
    distance_between_peaks[distance_between_peaks['distance'].isna()] = read_length # NA is because there is no peaks within the same scaffold or chromosome, so we set the distance to read_length for the below process
    distance_between_peaks = distance_between_peaks['distance'].to_numpy()
    distance_between_peaks[distance_between_peaks >= read_length] = read_length - 1
    outside_regions_for_reads_overlap_peaks = distance_between_peaks.sum()

    reads_counter = crpb.CountReadsPerBin([bam], bedFile=bed, numberOfProcessors=nproc)
    reads_at_peaks = reads_counter.run()
    total_reads = alignment.mapped
    total_reads_at_peaks = reads_at_peaks.sum(axis=0)
    frip = float(total_reads_at_peaks[0]) / total_reads

    return frip, reads_at_peaks, total_reads, outside_regions_for_reads_overlap_peaks


def create_frip_table_from_bed(
    samples_metadata,
    path_to_bed,
    path_to_data,
    genome_size,
    species,
    nproc,
    peak_protein_srun="",
    customized_metadata=False
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

    peaks = bf.read_table(path_to_bed, schema="bed").iloc[:, :3]
    num_peaks = [peaks.shape[0]] * len(bams)
    peaks_width = peaks["end"] - peaks["start"]
    total_bp_in_peaks = [peaks_width.sum()] * len(bams)


    total_reads = []
    frip_enrich = []
    samples_frips = []
    for bam in bams:
        result = calculate_frip(bam, path_to_bed, nproc=nproc)
        samples_frips.append(result[0])
        total_reads.append(result[2])
        outside_regions_for_reads_overlap_peaks = result[3]
        expected_prob_read_in_peak = (total_bp_in_peaks[0] + outside_regions_for_reads_overlap_peaks) / genome_size
        frip_enrich.append(result[0] / expected_prob_read_in_peak) # This is the ratio of observed reads in peaks / expected reads in peaks
        print(bam, f"frip calculation done:")
        print(frip_enrich[-1])

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
    frip_df = frip_df.rename(columns={"GSM_accession": "peaks-SRA"})

    return frip_df
