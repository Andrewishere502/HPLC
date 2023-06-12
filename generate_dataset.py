"""Generate a final dataset from the accumulated Report01.csv files.
Run after accumulate_report01s.py.
"""
import os
import config

import pandas as pd
import numpy as np

from logger import Logger


def get_chem_id(ret_time, chemmeta_df):
    """Return the chem_id associated with the retention range ret_time
    falls in.
    """
    for line in chemmeta_df.values:
        begin_ret_time, end_ret_time, chem_id = list(line)
        if ret_time >= float(begin_ret_time) and ret_time <= float(end_ret_time):
            return chem_id
    return  # no id found, return None


def get_meta_data(parent_directory, filename, meta_df):
    """Return the meta data for a file name found in a specific
    directory.
    """
    # Narrow down the meta_df by the parent folder name
    dir_df = meta_df[meta_df["ParentFolderName"] == parent_directory]
    # Now filter out all rows that do not match the file name
    file_df = dir_df[dir_df["FileName"] == filename]

    # Convert those values into a list. If there were no matching rows
    # found, catch the IndexError and return an list full of None.
    try:
        meta_data = list(file_df.values[0])
    except IndexError:
        meta_data = ["NA"] * len(meta_df.columns)
    return dict(zip(meta_df.columns, meta_data))

# Create an object to log errors that occur. Don't erase the contents
# when run to preserve error logs from accumulate_report01s.py
logger = Logger("", "log.txt", erase_on_init=False)

# Define a folder for where all intermediate input data should be.
# PIF stands for processed input folder.
PIF = config.PROCESSED_INPUT_FOLDER
# Raise an error if the ProcessedInput folder doesn't exist because
# that means the user hasn't run accumulate_report01s.py yet.
if not os.path.exists(PIF):
    raise OSError("The directory ProcessedInput does not exist. Please run accumulate_report01s.py to create and populate it.")

# Define a folder where the output will go. OF stands for output folder
OF = config.OUTPUT_FOLDER
# Make the output folder if it doesn't exist
if not os.path.exists(OF):
    print("Output directory did not exist: Created Output.")
    os.mkdir("Output")

# Get the chemical meta data for reference later
chemmeta_df = pd.read_csv("chemmeta.csv")
# Get a list of all unique chemical ids
all_chem_ids = list(set(chemmeta_df["ChemicalID"]))

# Get the sample meta data for reference later
meta_df = pd.read_csv("supermeta.csv")
# Sort the meta data for easier reference
meta_df = meta_df.sort_values(by=["FileName"])
# Get a list of the columns in meta_df
meta_columns = list(meta_df.columns)

# DataFrame where all the output will be located
output_df = pd.DataFrame(columns=(meta_columns + sorted(all_chem_ids)))

# Get all csv files by name into a list.
processed_filenames = [fn for fn in os.listdir(PIF) if fn[-4:] == ".csv"]
for processed_filename in processed_filenames:
    # Read the processed csv file into a DataFrame
    df = pd.read_csv(f"{PIF}/{processed_filename}")
    # Remove the acc_ and .csv parts from the processed_filename
    # to get the parent directory name of the files referenced within
    # the csv.
    parent_directory = processed_filename[4:].replace(".csv", "")
    print(len(set([*df["FileName"].values])))
    for filename in set([*df["FileName"].values]):
        # Get the meta data for this file
        meta_data = get_meta_data(parent_directory, filename, meta_df)

        # Assemble the output row for output_df
        output_row = meta_data
        output_row.update({chem_id: None for chem_id in sorted(all_chem_ids)})

        # Get a DataFrame for all the retention times and areas for
        # this file name
        chem_data = df[df["FileName"] == filename][["RetTime", "Area"]].values
        for ret_time, area in chem_data:
            chem_id = get_chem_id(ret_time, chemmeta_df)
            if chem_id != None:
                output_row.update({chem_id: area})
        
        # Update the output DataFrame
        output_df = output_df.append(output_row, ignore_index=True)

# Replace all na chem ids with 1
for chem_id in all_chem_ids:
    output_df[chem_id].replace(np.NaN, 1, inplace=True)

# Replace any nan in ProportionSyriaca with "na"
output_df["ProportionSyriaca"].replace(np.NaN, "na", inplace=True)

# Save the final dataset if it is not empty
if len(output_df.values) > 0:
    output_df.to_excel(f"{OF}/finalDataset.xlsx")
