import bioframe

def get_region_snippet(ctcf_peaks, ctcf_motifs, region):
    """
    Extracts CTCF peaks and motif overlaps within a specified genomic region.

    Parameters:
    - ctcf_peaks (str or file-like): Path to the BED file containing CTCF peak coordinates.
    - ctcf_motifs (str or file-like): Path to the file containing CTCF motif data.
    - region (str): Genomic region string in the format 'chrom:start-end' (e.g., 'chr1:10000-20000').

    Returns:
    - pandas.DataFrame: DataFrame containing the overlapping peaks and motifs with the following columns:
        - 'chrom': Chromosome name.
        - 'start': Start coordinate of the region.
        - 'end': End coordinate of the region.
        - 'strand': Strand information ('+' or '-').
        - 'mid': Midpoint of the region (calculated as (start + end) / 2).
    """
    # Read CTCF peaks and motifs into DataFrames
    ctcf = bioframe.read_table(ctcf_peaks, schema='bed')
    motif = bioframe.read_table(ctcf_motifs)

    # Rename columns in the motif DataFrame for consistency
    motif = motif.rename(columns={0: 'chrom', 1: 'start', 2: 'end', 
                                   3: 'name', 4: 'score', 5: 'pval', 6: 'strand'})

    # Find overlaps between CTCF peaks and motifs
    peaks_motifs = bioframe.overlap(ctcf, motif, how='inner')

    # Rename overlapping columns to avoid conflicts
    peaks_motifs = peaks_motifs.rename(columns={'strand': 'strand_', 'strand_': 'strand'})

    # Extract relevant columns
    dataframe = peaks_motifs[['chrom', 'start', 'end', 'strand']]

    # Parse the specified region and filter the DataFrame for that region
    region_start = bioframe.parse_region_string(region)[1]  # Extract start coordinate of the region
    region_dataframe = bioframe.select(dataframe, region, cols=['chrom', 'start', 'end'])

    # Add a 'mid' column with the midpoint of the region
    region_dataframe['mid'] = (region_dataframe.end + region_dataframe.start) / 2

    # Reset index for a clean DataFrame
    region_dataframe = region_dataframe.reset_index(drop=True)

    return region_dataframe
