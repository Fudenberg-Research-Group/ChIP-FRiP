import numpy as np
import pandas as pd
import os
import yaml
import argparse
from utils import fetch_metadata

parser = argparse.ArgumentParser(description='Fetch a metadata table')
parser.add_argument('config_path', type=str, help='Path to the configuration file.')
args = parser.parse_args()

with open(args.config_path, "r") as f:
    config = yaml.safe_load(f)

########################################################################
############################ Parameters ################################
########################################################################
    
dataset = config["parameters"]["dataset"]
geo_accession = config["parameters"]["geo_accession"]

path_to_accessions = config["input"]["path_to_accessions"]

output_dir = config["output"]["output_directory"]

########################################################################

accessions = np.loadtxt(path_to_accessions, dtype="str")

### fetch title from SRUN accession
experi_infos = []
gsm = []
for i, a in enumerate(accessions):
    while True:
        try:
            data = fetch_metadata(a)
            title = data[a]["title"].split(";")
        except TypeError:
            continue
        break
    experi_info = title[1].split(":")
    gsm.append(experi_info[0])
    experi_infos.append(experi_info[1])
    print(i, a, "title has been fetched")

df = pd.DataFrame(
    {"SRUN": accessions, "Experiment": experi_infos, "GSM_accession": gsm}
)

### fetch attributes from GSM accession
condition = []
antibody = []
celltype = []
organism = []

KEY_WORDS = ["gen", "cell", "organism", "anti"]
attribute_keys = {}
for i, row in df.iterrows():
    a = row["GSM_accession"].replace(" ", "")

    # Step1: fetch the dictionary of attributes
    while True:
        try:
            data = fetch_metadata(a)
            attributes = data[a]["samples"][list(data[a]["samples"].keys())[0]]["attributes"]
        except (TypeError, KeyError):
            continue
        break

    # Step2: fetch attributes' keys in dictionary
    if len(attribute_keys) == 0:
        key_df = pd.DataFrame({"Keys": attributes.keys()})
        for kw in KEY_WORDS:
            try:
                attribute_keys[kw] = key_df[key_df["Keys"].str.contains(kw, case=False)].iloc[0, 0]
            except IndexError:
                print(f"There is no attributes match with the substring {kw}, so we add its title here, and you can check manually")
                attribute_keys[kw] = 'title'
    
    # Step3: add entries to each attribute column
    for j, col in enumerate([condition, celltype, organism, antibody]):
        try:
            col.append(attributes[attribute_keys[KEY_WORDS[j]]])
        except KeyError:
            col.append(row["Experiment"])

    print(i, accessions[i], "attributes have been fetched")

# create a dataframe
df["Condition"] = condition
df["Antibody"] = antibody
df["Celltype"] = celltype
df["Organism"] = organism

df["Peak_ChIP"] = ["N/A"] * len(df)
df["author_year"] = [dataset] * len(df)
df["GEO"] = [geo_accession] * len(df)

# reform dataframe
cols = [
    "Organism",
    "Celltype",
    "Condition",
    "Antibody",
    "Peak_ChIP",
    "author_year",
    "SRUN",
    "GSM_accession",
    "GEO",
    "Experiment",
]

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

df[cols].to_csv(output_dir + '/metadata.txt', sep="\t", index=False)
