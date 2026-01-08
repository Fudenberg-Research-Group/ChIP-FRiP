# Import necessary libraries

import bioframe
import numpy as np

# Define utility functions

# Function: Convert CTCF occupancy based on scores
def convert_ctcf_occupancy(ctcf_bed_df, score_bed_path, max_occup=0.9):
    """
    Convert CTCF occupancy from bed file scores.

    Parameters:
    - ctcf_bed_df (DataFrame): DataFrame with CTCF binding regions. Example:
        chrom  | start  | end    | mid     | strand
        -------|--------|--------|---------|-------
        chr1   | 340000 | 341000 | 340500  | +
    - score_bed_path (str): Path to BED file containing score values.
        Example: '/path/to/CTCF_score.bed'
    - max_occup (float): Maximum occupancy to scale scores to. Default is 0.9.

    Returns:
    - DataFrame: Updated DataFrame with occupancy values added.
    """
    bed_file = bioframe.read_table(score_bed_path, schema='bed')
    dataframe_ctcf = bioframe.overlap(ctcf_bed_df, bed_file, how='inner')[
        ['chrom', 'start', 'end', 'mid', 'strand', 'score_']
    ]
    max_score = dataframe_ctcf['score_'].max()
    dataframe_ctcf['occupancy'] = (dataframe_ctcf['score_'] / max_score) * max_occup
    return dataframe_ctcf


# Function: Make CTCF occupancy (predicted or converted)
def make_ctcf_occupancy(data_frame, bed_df, predicted=True):
    """
    Generate occupancy data for CTCF regions.

    Parameters:
    - data_frame (DataFrame): CTCF binding data. See `convert_ctcf_occupancy` for structure.
    - bed_df (DataFrame): BED file data. Example:
        chrom  | start  | end    | score_
        -------|--------|--------|--------
        chr1   | 340000 | 341000 | 0.85
    - predicted (bool): Whether to use prediction or convert from scores. Default is True.

    Returns:
    - DataFrame: DataFrame with occupancy data added.
    """
    if predicted:
        # Placeholder for ML prediction function
        bed_with_occupancy = predict_ctcf_occupancy(data_frame, ctcf_frequency)
    else:
        bed_with_occupancy = convert_ctcf_occupancy(data_frame)
    
    return bed_with_occupancy


# Function: Refine occupancy for overlapping sites
def get_refined_occupancy(region_dataframe, region, lattice_size=250):
    """
    Refine occupancy data for overlapping lattice locations.

    Parameters:
    - region_dataframe (DataFrame): DataFrame with occupancy data. Example:
        chrom  | start  | end    | mid     | strand | occupancy
        -------|--------|--------|---------|--------|-----------
        chr1   | 340000 | 341000 | 340500  | +      | 0.72
    - region (str): Genomic region string (e.g., 'chr1:34000000-34500000').
    - lattice_size (int): Size of lattice in base pairs. Default is 250.

    Returns:
    - DataFrame: Refined occupancy data. Example:
        chrom  | start  | end    | mid     | strand | lattice_loc | occupancy
        -------|--------|--------|---------|--------|-------------|-----------
        chr1   | 340000 | 341000 | 340500  | +      | 136000      | 0.85
    """
    region_start, region_end = bioframe.parse_region_string(region)[1:3]
    number_lattice = int((region_end - region_start) / lattice_size)
    region_dataframe['lattice_loc'] = ((region_dataframe['mid'] - region_start) // lattice_size).astype('int')
    region_dataframe['lattice_loc'] = np.clip(region_dataframe['lattice_loc'], a_min=None, a_max=number_lattice-1)

    result_c = region_dataframe.groupby(['lattice_loc', 'strand'])['predicted_occupancy'].apply(
        lambda x: 1 - ((1 - x).prod())
    ).reset_index()
    result = result_c.merge(region_dataframe.drop_duplicates(['lattice_loc', 'strand']),
                            on=['lattice_loc', 'strand'], how='left')
    result = result.rename(columns={'predicted_occupancy_x': 'predicted_occupancy'})
    result = result[['chrom', 'start', 'end', 'mid', 'strand', 'lattice_loc', 'predicted_occupancy']]
    return result


# Function: Calculate lattice size
def get_lattice_size(region, lattice_site=250):
    """
    Calculate lattice size for a given genomic region.

    Parameters:
    - region (str): Genomic region string (e.g., 'chr1:34000000-34500000').
    - lattice_site (int): Size of each lattice site in base pairs. Default is 250.

    Returns:
    - int: Number of lattice sites within the region.
    """
    region_start = bioframe.parse_region_string(region)[1]
    region_end = bioframe.parse_region_string(region)[2]
    lattice_size = (region_end - region_start) // lattice_site
    return lattice_size


# Function: Generate CTCF lists for simulation
def get_ctcf_list(target_dataframe, parameters, insert_on='bound_time'):
    """
    Extract CTCF positions and related parameters for simulations.

    Parameters:
    - target_dataframe (DataFrame): DataFrame with occupancy and lattice locations.
        chrom  | start  | end    | mid     | strand | lattice_loc | predicted_occupancy
        -------|--------|--------|---------|--------|-------------|-----------
        chr1   | 340000 | 341000 | 340500  | +      | 136000      | 0.85
    - parameters (dict): Dictionary of simulation parameters. Example:
        {
            'CTCF_offtime': [180],
            ...
        }

    Returns:
    - Tuple: (CTCF_left_positions, CTCF_right_positions, ctcf_loc_list, ctcf_lifetime_list, ctcf_offtime_list)
        - CTCF_left_positions (np.array): Positions of '+' strand CTCFs.
        - CTCF_right_positions (np.array): Positions of '-' strand CTCFs.
        - ctcf_loc_list (list): Sorted positions of all CTCFs.
        - ctcf_lifetime_list (np.array): Lifetimes of CTCFs.
        - ctcf_offtime_list (np.array): Off-times of CTCFs.
    """
    ctcfleftlist = list(target_dataframe[target_dataframe['strand'] == '+']['lattice_loc'])
    ctcfrightlist = list(target_dataframe[target_dataframe['strand'] == '-']['lattice_loc'])
    ctcf_loc_list = list(np.sort(ctcfrightlist + ctcfleftlist))
    CTCF_right_positions = np.array(ctcfrightlist)
    CTCF_left_positions = np.array(ctcfleftlist)
    ctcf_occup_list = list(target_dataframe['predicted_occupancy'])
    ctcf_occup_array = np.array(ctcf_occup_list)
    if insert_on=='bound_time':
        CTCF_offtime_val = parameters['CTCF_offtime'][0]
        ctcf_lifetime_list = CTCF_offtime_val * ctcf_occup_array / (1 - ctcf_occup_array)
        ctcf_offtime_list = np.ones(len(ctcf_loc_list)) * CTCF_offtime_val
    else:
        CTCF_lifetime_val = parameters['CTCF_lifetime'][0]
        ctcf_lifetime_list = np.ones(len(ctcf_loc_list)) * CTCF_lifetime_val
        ctcf_offtime_list = CTCF_lifetime_val * (1-ctcf_occup_array)/ctcf_occup_array
        
    return CTCF_left_positions, CTCF_right_positions, ctcf_loc_list, ctcf_lifetime_list, ctcf_offtime_list
