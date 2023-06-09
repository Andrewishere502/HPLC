"""Generate a final dataset from the accumulated Report01.csv files.
Run after accumulate_report01s.py.
"""
import os

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


def get_meta_data(directory, filename, meta_df):
    """Return the meta_data associated with the given filename.
    Return 1 if this whole directory is missing meta data.
    Return 2 if this specific file is missing meta data.
    """
    # Get all of the meta entries associated the given directory
    dir_df = meta_df[meta_df["ParentFolderName"] == directory]
    if len(dir_df.values) == 0:
        logger.log(f"Directory {directory} does not have any meta data.")
        return 1
    
    # Remove the .csv suffix from filename if it is there
    filename = filename.replace(".csv", "")

    # Grab the meta data in this directory that matches the given filename
    results = dir_df[dir_df["FileName"] == filename].values
    if len(results) == 1:
        # Get the metadata and label each value with it's column
        meta_data = dict(zip(dir_df.columns, results[0]))
        return meta_data

    # Got no matches
    elif len(results) < 1:
        # This doesn't need to be a fatal error because we can just
        # skip the file.
        logger.log(f"There is no meta data for file named {filename} in {directory}.")
        # If we get here, then there was no good result from querrying the
        return 2

     # Got more than 1 match
    else:
        # As of yet, there is no way to tell which file is the one we
        # want, so this should be a fatal error.
        raise ValueError(f"There are {len(results)} files in {directory} with name {filename}. Duplicates are not allowed.")
    

# Create an object to log errors that occur. Don't erase the contents
# when run to preserve error logs from accumulate_report01s.py
logger = Logger("", "log.txt", erase_on_init=False)

# Define a folder for where all intermediate input data should be.
# PIF stands for processed input folder.
PIF = "ProcessedInput"
if not os.path.exists(PIF):
    raise OSError("The directory ProcessedInput does not exist. Please run accumulate_report01s.py to create and populate it.")

# Define a folder where the output will go. OF stands for output folder
OF = "Output"
if not os.path.exists(OF):
    print("Output directory did not exist: Created Output.")
    os.mkdir("Output")

# Get the chemical meta data for reference later
chemmeta_df = pd.read_csv(f"chemmeta.csv")
# Get a list of all unique chemical ids
all_chem_ids = list(set(chemmeta_df["ChemicalID"]))

# Get the sample meta data for reference later
meta_df = pd.read_csv(f"supermeta.csv")
# Sort the meta data for easier reference
meta_df = meta_df.sort_values(by=["FileName"])
# Get a list of the columns in meta_df
meta_columns = list(meta_df.columns)


# Get a list of all the processed data folders, i.e. they have been
# output by accumulate_report01s.py
processed_data_folders = os.listdir(PIF)
# Exclude the meta data files and hidden files like .DS_Store
processed_data_folders = [folder_name for folder_name in processed_data_folders
                          if folder_name[-4:] != ".csv" and folder_name[0] != "."]
for folder_name in processed_data_folders:
    # This is the folder where all of the input data should be
    accumulated_folder_path = f"{PIF}/{folder_name}"

    # Get a list of all the file names containing raw data
    filenames = [filename for filename in os.listdir(accumulated_folder_path)
                 if filename[-6:-4] == ".D"]

    # Create the output dataframe
    columns = meta_columns + sorted(all_chem_ids)
    output_df = pd.DataFrame(columns=columns)

    # Loop through each filename containing peak, retention time and area
    # data.
    for filename in filenames:
        # Get the data from the file
        df = pd.read_csv(f"{accumulated_folder_path}/{filename}")

        # Get the meta data associated with this file via filename,
        # excluding the 'acc_' portion.
        meta_data = get_meta_data(folder_name[4:], filename, meta_df)
        
        # Catch non-fatal errors from get_meta_data call
        if meta_data == 1:
            # Skip this whole folder and move on. Check log.txt for error log.
            break
        elif meta_data == 2:
            # Skip this file and move on. Check log.txt for error log.
            continue

        # Put together a row summarizing this whole file to be output in
        # output_df.
        output_row = meta_data
        output_row.update({chem_id: None for chem_id in sorted(all_chem_ids)})
        for row in df.values:
            peak_num, ret_time, area = list(row)
            # Get the chem meta data associated with the retention time
            chem_id = get_chem_id(ret_time, chemmeta_df)
            # If a chem_id was found, record the area for it in output_row
            if chem_id != None:
                output_row.update({chem_id: area})

        # Update the output_df with the new row
        output_df = output_df.append(output_row, ignore_index=True)


    # Replace missing data for all chem_id columns
    for chem_id in all_chem_ids:
        output_df[chem_id].replace(np.NaN, 1, inplace=True)

    # Save the final dataset if it actually has data in it
    if len(output_df.values) > 0:
        output_df.to_excel(f"{OF}/{folder_name}_finalDataset.xlsx")
