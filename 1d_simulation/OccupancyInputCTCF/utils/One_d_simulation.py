import numpy as np
import matplotlib.pylab as plt
import os
import h5py 
import time
import sys

from .lattice_translocators import LEFTranslocator, LEFTranslocatorDynamicBoundary
sys.path.append('/home1/rahmanin/start/polychrom/projects/Site_wise_occupancy/OccupancyInputCTCF/')
from . import simfuncs as rfuncs

def Perform_1d_simulation(paramdict, ctcf_params, trajectory_length, output_directory):
    """
    Perform a 1D lattice simulation of LEF translocation dynamics.

    Parameters:
    - paramdict (dict): Dictionary containing simulation parameters. Expected keys:
        - 'monomers_per_replica': Number of monomers per replica (int).
          Example: 500.
        - 'sites_per_monomer': Number of sites per monomer (int).
          Example: 2.
        - 'number_of_replica': Total number of replicas (int).
          Example: 10.
        - 'LEF_separation': Separation between LEFs in lattice units (int).
          Example: 10.
    - ctcf_params (tuple): Parameters for CTCF positions and dynamics. Contains:
        - CTCF_left_positions (np.array): Array of left CTCF positions (lattice indices).
        - CTCF_right_positions (np.array): Array of right CTCF positions (lattice indices).
        - ctcf_loc_list (list): List of all CTCF positions sorted.
        - ctcf_lifetime_list (np.array): Lifetimes of CTCF binding at positions.
        - ctcf_offtime_list (np.array): Off-times for CTCF unbinding.
    - trajectory_length (int): Number of simulation steps to perform.
      Example: 100.
    - output_directory (str): Path to save simulation output files.
      Example: "/output/simulation_results".

    Returns:
    - None: Writes simulation results to an HDF5 file in `output_directory`.

    Output Files:
    - `LEFPositions.h5`: HDF5 file containing LEF positions, CTCF positions, and related data.
    """
    # Extract parameters for lattice and monomer configurations
    monomers_per_replica = paramdict['monomers_per_replica']
    sites_per_monomer = paramdict['sites_per_monomer']
    sites_per_replica = monomers_per_replica * sites_per_monomer
    typedict = {'A': 0}  # Simple monomer type mapping (expandable)
    monomer_types = np.zeros(monomers_per_replica, dtype=int)  # Initialize monomer types
    site_types = np.repeat(monomer_types, sites_per_monomer)  # Repeat for lattice sites
    
    # Unpack CTCF parameters
    CTCF_left_positions, CTCF_right_positions, ctcf_loc_list, ctcf_lifetime_list, ctcf_offtime_list = ctcf_params

    ########### Lattice simulation setup ###########
    num_dummy_steps = int(0.25 * trajectory_length)  # Equilibration steps for LEFs
    blocksteps = 5  # Number of steps per block for data writing
    bins = np.linspace(0, trajectory_length, blocksteps, dtype=int)  # Define simulation bins
    N = paramdict['monomers_per_replica'] * paramdict['number_of_replica']
    LEFNum = N // paramdict['LEF_separation']  # Compute LEF number based on separation

    # Initialize LEF translocator with dynamic boundaries
    translocator = rfuncs.make_translocator(
        LEFTranslocatorDynamicBoundary, 
        site_types,
        CTCF_left_positions,
        CTCF_right_positions, 
        ctcf_loc_list,
        ctcf_lifetime_list,
        ctcf_offtime_list,
        **paramdict
    )

    # Create HDF5 file to store simulation results
    with h5py.File(f"{output_directory}/LEFPositions.h5", mode='w') as myfile:
        # Create datasets for LEF positions and CTCF site information
        dset = myfile.create_dataset(
            "positions", 
            shape=(trajectory_length, LEFNum, 2), 
            dtype=np.int32, 
            compression="gzip"
        )
        dset_ctcf_sites_right = myfile.create_dataset(
            "CTCF_sites_right",
            shape=(len(CTCF_right_positions)), 
            compression="gzip", 
            data=CTCF_right_positions.copy()
        )
        dset_ctcf_sites_left = myfile.create_dataset(
            "CTCF_sites_left",
            shape=(len(CTCF_left_positions)), 
            compression="gzip",
            data=CTCF_left_positions.copy()
        )
        dset_ctcf_positions_right = myfile.create_dataset(
            "CTCF_positions_right",
            shape=(trajectory_length, len(CTCF_right_positions), 1), 
            compression="gzip"
        )
        dset_ctcf_positions_left = myfile.create_dataset(
            "CTCF_positions_left",
            shape=(trajectory_length, len(CTCF_left_positions), 1), 
            compression="gzip"
        )
        
        # Run initial dummy steps to equilibrate LEFs
        translocator.steps(num_dummy_steps)

        # Perform simulation in blocks and write results
        for st, end in zip(bins[:-1], bins[1:]):
            cur = []
            ctcf_right_cur = []
            ctcf_left_cur = []
            for i in range(st, end):
                translocator.step()  # Simulate one step
                cur.append(translocator.LEFs.copy())

                if len(CTCF_right_positions) > 0:
                    ctcf_positions_right = (translocator.stallProbRight)[CTCF_right_positions] * 1
                    ctcf_right_cur.append(ctcf_positions_right.reshape(len(ctcf_positions_right), 1))
                else:
                    ctcf_right_cur.append(np.empty((0, 1)))  # Append an empty (0,1) array for consistency

                if len(CTCF_left_positions) > 0:
                    ctcf_positions_left = (translocator.stallProbLeft)[CTCF_left_positions] * 1
                    ctcf_left_cur.append(ctcf_positions_left.reshape(len(ctcf_positions_left), 1))
                else:
                    ctcf_left_cur.append(np.empty((0, 1)))
                    
            # Convert collected data to arrays and write to HDF5
            cur = np.array(cur)
            ctcf_right_cur = np.array(ctcf_right_cur)
            ctcf_left_cur = np.array(ctcf_left_cur)
            dset[st:end] = cur
            dset_ctcf_positions_right[st:end] = ctcf_right_cur
            dset_ctcf_positions_left[st:end] = ctcf_left_cur
        
        # Add attributes to HDF5 file for metadata
        myfile.attrs["N"] = N * paramdict['sites_per_monomer']
        myfile.attrs["LEFNum"] = LEFNum

