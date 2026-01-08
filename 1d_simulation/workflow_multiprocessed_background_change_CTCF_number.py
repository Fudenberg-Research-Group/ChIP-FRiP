import os
import json
import numpy as np
np.random.seed(42)
import ctypes
import argparse

# Import utility modules
import OccupancyInputCTCF.utils.convert as convert
import OccupancyInputCTCF.utils.makeparams as params
import OccupancyInputCTCF.utils.One_d_simulation as simulation
# ================ INPUT PARAMETERS =====================

parser = argparse.ArgumentParser(description="LEF simulation")
parser.add_argument("config_path", type=str, help="Path to the configuration file.")
parser.add_argument("ctcf_num", type=str, help="Relative number of CTCF")


args = parser.parse_args()

# Simulation parameters
with open(args.config_path, 'r') as json_file:
    paramdict = json.load(json_file)

# Define genomic region and data files
region = 'chr1:25000000-26000000' # use 1Mb

# Flags for simulation
run_md_sim = True  # Flag for MD simulation
use_predicted_occupancy = True  # Use predicted occupancy

output_prefix = './frip_simulation_change_CTCF/'

# =================== WORKFLOW ===========================
print("Starting workflow...")

CTCF_left_positions = np.load('./data/different_CTCF_number_simulation_configs/ctcf_params/CTCF_left_positions.npy')
CTCF_right_positions = np.load('./data/different_CTCF_number_simulation_configs/ctcf_params/CTCF_right_positions.npy')
ctcf_loc_list = np.load('./data/different_CTCF_number_simulation_configs/ctcf_params/ctcf_loc_list.npy')
ctcf_lifetime_list = np.load(f'./data/different_CTCF_number_simulation_configs/ctcf_params/ctcf_lifetime_list({args.ctcf_num}).npy')
ctcf_offtime_list = np.load(f'./data/different_CTCF_number_simulation_configs/ctcf_params/ctcf_offtime_list({args.ctcf_num}).npy')

# Running 1D Simulation
import multiprocessing as mp
from functools import partial
print("Running 1D Simulation...")
output_path = f'{output_prefix}/cohesin_loc/'

############## Parameters #########################
trajectory_length = 4000
n = 100 # number of simulations in multiprocessing
###################################################

file_name = params.paramdict_to_filename(paramdict)
output_directory = output_path + 'folder_' + f'ctcfnum{args.ctcf_num}' + file_name.split('file_')[1]
os.makedirs(output_directory, exist_ok=True)

paramdict['monomers_per_replica'] = convert.get_lattice_size(region, lattice_site=250)//10 # 1Mb
ctcf_params = CTCF_left_positions, CTCF_right_positions, ctcf_loc_list, ctcf_lifetime_list, ctcf_offtime_list

outputt_dirs = []
for sim_id in range(1, n+1):
    file_name = f"simulation_{sim_id}"
    output_directory_partial = os.path.join(output_directory, f"{file_name}")
    os.makedirs(output_directory_partial, exist_ok=True)
    outputt_dirs.append(output_directory_partial)
    
def Perform_1d_simulation(output_directory, paramdict, ctcf_params, trajectory_length):
    return simulation.Perform_1d_simulation(paramdict, ctcf_params, trajectory_length, output_directory)
# use partial to set up some parameters of the function and leave output dir later
partially_set_Perform_1d_simulation = partial(Perform_1d_simulation, paramdict=paramdict, ctcf_params=ctcf_params, trajectory_length=trajectory_length)


master_seed = 42
num_processes = mp.cpu_count()
worker_id_counter = mp.Value(ctypes.c_int, 0)
initializer_args = (worker_id_counter, master_seed)
def worker_init_seeded(shared_counter, master_seed_value):
    # Atomically get a unique ID for this worker
    with shared_counter.get_lock():
        worker_id = shared_counter.value
        shared_counter.value += 1  # Increment for the next worker

    # Derive a unique but deterministic seed for this worker
    worker_seed = master_seed_value + worker_id
    # Seed this worker's NumPy random number generator
    np.random.seed(worker_seed)

with mp.Pool(processes=mp.cpu_count(), initializer=worker_init_seeded, initargs=initializer_args) as pool:
    pool.map(
        partially_set_Perform_1d_simulation, outputt_dirs
    )
    
print("1D simulation complete. Output at", output_path)
 