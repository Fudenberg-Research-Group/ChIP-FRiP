from .lattice_translocators import LEFTranslocator, LEFTranslocatorDynamicBoundary
import numpy as np

def make_site_array(site_types, 
                    values, 
                    at_ids=None, 
                    number_of_replica=1, 
                    **kwargs):
    """
    Constructs an array based on site types and assigns values to specific positions.

    Parameters:
    - site_types (np.ndarray): Array defining site types across the lattice.
    - values (list): Values corresponding to each site type.
    - at_ids (np.ndarray, optional): Indices specifying where values should be applied.
    - number_of_replica (int): Number of replicas to tile the array.

    Returns:
    - np.ndarray: Array with values assigned based on site types and positions.
    """
    assert site_types.max() < len(values), (
        'Number of values (%d) incompatible with number of site types (%d)' 
        % (len(values), site_types.max())
    )
    prop_array = np.zeros(len(site_types), dtype=np.double)
    
    for i, value in enumerate(values):
        prop_array[site_types == i] = value
        
    if isinstance(at_ids, np.ndarray):
        mask = np.zeros(len(site_types), dtype=bool)
        mask[at_ids] = True
        prop_array[~mask] = 0
        
    return np.tile(prop_array, number_of_replica)

def make_CTCF_arrays(site_types,
                     CTCF_left_positions,
                     CTCF_right_positions,
                     CTCF_facestall,
                     CTCF_backstall,
                     **kwargs):
    """
    Creates arrays for CTCF stalling probabilities.

    Parameters:
    - site_types (np.ndarray): Lattice site type array.
    - CTCF_left_positions (np.ndarray): Indices of left-facing CTCF sites.
    - CTCF_right_positions (np.ndarray): Indices of right-facing CTCF sites.
    - CTCF_facestall (list): Stalling probabilities for face-on collisions.
    - CTCF_backstall (list): Stalling probabilities for back collisions.

    Returns:
    - list: Two arrays [stall_left_array, stall_right_array] representing stalling probabilities.
    """
    if len(CTCF_left_positions) == 0:
        CTCF_left_positions = None
    if len(CTCF_right_positions) == 0:
        CTCF_right_positions = None
        
    stall_left_array = make_site_array(site_types, CTCF_facestall, at_ids=CTCF_left_positions, **kwargs)
    stall_right_array = make_site_array(site_types, CTCF_facestall, at_ids=CTCF_right_positions, **kwargs)
    
    stall_left_array += make_site_array(site_types, CTCF_backstall, at_ids=CTCF_right_positions, **kwargs)
    stall_right_array += make_site_array(site_types, CTCF_backstall, at_ids=CTCF_left_positions, **kwargs)
    
    return [stall_left_array, stall_right_array]

def make_CTCF_dynamic_arrays_site_specific(site_types,
                                           CTCF_loc_list,
                                           CTCF_lifetime_list,
                                           CTCF_offtime_list,
                                           sites_per_monomer,
                                           velocity_multiplier,
                                           number_of_replica,
                                           **kwargs):
    """
    Generates dynamic arrays for CTCF birth and death rates at specific locations.

    Parameters:
    - site_types (np.ndarray): Lattice site type array.
    - CTCF_loc_list (list): Indices of CTCF locations.
    - CTCF_lifetime_list (list): Lifetimes for each CTCF site.
    - CTCF_offtime_list (list): Off-times for each CTCF site.
    - sites_per_monomer (int): Number of sites per monomer.
    - velocity_multiplier (float): Multiplier to adjust simulation speed.
    - number_of_replica (int): Number of replicas to tile the array.

    Returns:
    - list: Arrays [CTCF_death_array, CTCF_birth_array] for death and birth rates.
    """
    CTCF_lifetime_array = np.zeros(len(site_types))
    CTCF_offtime_array = np.zeros(len(site_types))

    for i, value in enumerate(CTCF_lifetime_list):
        CTCF_lifetime_array[CTCF_loc_list[i]] = CTCF_lifetime_list[i]
        CTCF_offtime_array[CTCF_loc_list[i]] = CTCF_offtime_list[i]
    
    CTCF_death_array = 1. / CTCF_lifetime_array / sites_per_monomer
    CTCF_birth_array = 1. / CTCF_offtime_array / sites_per_monomer

    return [np.tile(CTCF_death_array, number_of_replica), 
            np.tile(CTCF_birth_array, number_of_replica)]

def make_LEF_arrays(site_types,
                    LEF_lifetime,
                    LEF_stalled_lifetime,
                    LEF_birth,
                    LEF_pause,
                    sites_per_monomer,
                    velocity_multiplier,
                    **kwargs):
    """
    Creates arrays for LEF properties such as lifetime, birth, and pausing.

    Parameters:
    - site_types (np.ndarray): Lattice site type array.
    - LEF_lifetime (list): Lifetimes of LEFs.
    - LEF_stalled_lifetime (list): Stalled lifetimes of LEFs.
    - LEF_birth (list): Birth probabilities for LEFs.
    - LEF_pause (list): Pausing probabilities for LEFs.
    - sites_per_monomer (int): Number of sites per monomer.
    - velocity_multiplier (float): Multiplier to adjust simulation speed.

    Returns:
    - list: Arrays [death_array, stalled_death_array, birth_array, pause_array].
    """
    lifetime_array = make_site_array(site_types, LEF_lifetime, **kwargs)
    stalled_lifetime_array = make_site_array(site_types, LEF_stalled_lifetime, **kwargs)
    birth_array = make_site_array(site_types, LEF_birth, **kwargs)
    pause_array = make_site_array(site_types, LEF_pause, **kwargs)
    
    death_array = 1. / lifetime_array / sites_per_monomer
    stalled_death_array = 1. / stalled_lifetime_array / sites_per_monomer

    return [death_array, stalled_death_array, birth_array, pause_array]

def make_translocator(extrusion_engine, 
                      site_types,
                      CTCF_left_positions,
                      CTCF_right_positions,
                      CTCF_loc_list,
                      CTCF_lifetime_list,
                      CTCF_offtime_list,
                      **kwargs):
    """
    Initializes a LEFTranslocator or LEFTranslocatorDynamicBoundary with appropriate arrays.

    Parameters:
    - extrusion_engine (callable): Engine class for creating the translocator.
    - site_types (np.ndarray): Lattice site type array.
    - CTCF_left_positions (np.ndarray): Indices of left-facing CTCF sites.
    - CTCF_right_positions (np.ndarray): Indices of right-facing CTCF sites.
    - CTCF_loc_list (list): Indices of dynamic CTCF locations.
    - CTCF_lifetime_list (list): Lifetimes for dynamic CTCF sites.
    - CTCF_offtime_list (list): Off-times for dynamic CTCF sites.
    - kwargs (dict): Additional parameters (e.g., LEF_separation, velocity_multiplier, etc.).

    Returns:
    - LEFTranslocator or LEFTranslocatorDynamicBoundary: Configured translocator object.
    """
    LEF_separation = kwargs['LEF_separation']    
    velocity_multiplier = kwargs['velocity_multiplier'] 
    sites_per_monomer = kwargs['sites_per_monomer'] 
    number_of_replica = kwargs['number_of_replica'] 
    monomers_per_replica = kwargs['monomers_per_replica'] 

    number_of_monomers = number_of_replica * monomers_per_replica
    number_of_LEFs = number_of_monomers // LEF_separation
    sites_per_replica = monomers_per_replica * sites_per_monomer

    assert len(site_types) == sites_per_replica, (
        "Site type array (%d) doesn't match replica lattice size (%d)"
        % (len(site_types), sites_per_replica)
    )

    # Create arrays
    LEF_arrays = make_LEF_arrays(site_types, **kwargs)
    CTCF_arrays = make_CTCF_arrays(site_types, CTCF_left_positions, CTCF_right_positions, **kwargs)
    CTCF_dynamic_arrays = make_CTCF_dynamic_arrays_site_specific(site_types, CTCF_loc_list, 
                                                                  CTCF_lifetime_list, 
                                                                  CTCF_offtime_list, **kwargs)

    LEFTran = extrusion_engine(number_of_LEFs, *LEF_arrays, *CTCF_arrays, *CTCF_dynamic_arrays)

    if not isinstance(LEFTran, LEFTranslocatorDynamicBoundary):
        LEFTran.stallProbLeft = 1 - (1 - LEFTran.stallProbLeft) ** (1. / velocity_multiplier)
        LEFTran.stallProbRight = 1 - (1 - LEFTran.stallProbRight) ** (1. / velocity_multiplier)

    return LEFTran
