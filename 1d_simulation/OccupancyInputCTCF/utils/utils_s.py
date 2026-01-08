import bioframe

import numpy as np

def region_data_frame(dataframe, region, lattice_size=250):
    """
    Extracts and processes a specified genomic region from a dataframe.

    Parameters:
    ----------
    dataframe : DataFrame
        The input dataframe containing genomic data with columns such as 'chrom', 'start', 'end', and 'mid'.
    region : str
        A string specifying the genomic region in the format 'chrom:start-end'.
    lattice_size : int, optional
        The size of each lattice or bin to segment the region, default is 250.

    Returns:
    -------
    DataFrame
        A dataframe filtered to the specified region with an added column 'lattice_loc' indicating
        the lattice location of each row based on 'mid' position and lattice size.
    """
    region_start = bioframe.parse_region_string(region)[1]
    region_dataframe = bioframe.select(dataframe, region, cols=['chrom', 'start', 'end'])
    region_dataframe['mid']=(region_dataframe.end+region_dataframe.start)/2
    region_dataframe['lattice_loc'] = ((region_dataframe['mid'] - region_start) // lattice_size).astype('int')
    region_dataframe = region_dataframe.reset_index(drop=True)
    return region_dataframe

def make_region_occupancy(file):
    df = pandas.read_csv(file)
    result_c = df.groupby(['lattice_loc', 'strand'])['predicted_occupancy'].apply(lambda x: 1-((1 - x).prod())).reset_index()
    result = result_c.merge(df.drop_duplicates(['lattice_loc', 'strand']), on=['lattice_loc', 'strand'], how='left')
    result = result.rename(columns={'predicted_occupancy_x':'predicted_occupancy'})
    result = result[['chrom','start','end','mid','strand','lattice_loc','predicted_occupancy']]
    return result 



def contact_map_from_lefs(dset, sites_per_replica):
    
    lef_array = np.mod(dset.reshape((-1, 2)), sites_per_replica)
    lef_array = lef_array[lef_array[:,1] > lef_array[:,0]]
    
    lef_map = np.histogram2d(lef_array[:,0], lef_array[:,1], np.arange(sites_per_replica+1))[0]
    
    return (lef_map + lef_map.T)

def chip_seq_from_lef(lef_positions, site_number_per_replica, min_time=0):
    lef_lefts = lef_positions[min_time:,:,0].flatten()
    lef_rights = lef_positions[min_time:,:,1].flatten() 
    lef_positions_aray = np.hstack((lef_lefts,lef_rights))
    hist,hist_ary = np.histogram(  np.mod( lef_positions_aray , site_number_per_replica ), np.arange(0, site_number_per_replica+1, 1))
    return hist

def peak_positions(boundary_list, window_sizes=[1]):
    """
    Calculate peak positions based on a boundary_list within window_sizes.

    Args:
        boundary_list (list): List of boundary values.
        window_sizes (list, optional): List of window sizes. Defaults to [1].

    Returns:
        np.ndarray: Array containing peak positions.
    """
    peak_monomers = np.array([])

    for i in window_sizes:
        inds_to_add = [boundary + i for boundary in boundary_list]
        peak_monomers = np.hstack((peak_monomers, inds_to_add))

    return peak_monomers.astype(int)

def FRiP(num_sites_t, lef_positions, peak_positions ):
    
    hist,edges = np.histogram(  lef_positions  , np.arange(num_sites_t+1) )
    return np.sum(hist[peak_positions] )/len(lef_positions)

def chip_seq_from_ctcf(ctcf_array_right, ctcf_array_right_sites, ctcf_array_left, ctcf_array_left_sites, site_number_per_replica):
    ctcfrightary = np.concatenate([arr.flatten()*ctcf_array_right_sites for arr in ctcf_array_right if arr.size > 0])
    ctcfleftary = np.concatenate([arr.flatten()*ctcf_array_left_sites for arr in ctcf_array_left if arr.size >0])
    ctcfs = np.concatenate([ctcfrightary[ctcfrightary>0], ctcfleftary[ctcfleftary>0]])
    ctcfhist, hist_array = np.histogram(ctcfs, np.arange(0, site_number_per_replica+1, 1))
    common_list = np.intersect1d(ctcf_array_right_sites, ctcf_array_left_sites)
    for elements in common_list:
        ctcfhist[elements] = ctcfhist[elements]/2
    return ctcfhist
