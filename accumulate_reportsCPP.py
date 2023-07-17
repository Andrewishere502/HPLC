"""Take the peak number, retention time, and area from each of the
Report01.csv files scattered throughout an HPLC Chromatogram directory,
and put that data in a file named after the directory the data was
found in.
"""
# Standard library imports
import os

# Third party imports
import pandas as pd

# Import from our app
from logger import Logger
import config
import wcidNEW


# Create an object to log errors that occur
logger = Logger("", "log.txt", erase_on_init=True)

# RawInput is the folder to put all of the raw data into. RIF stands
# for raw input folder
# RIF = "RawInput"
RIF = config.RAW_INPUT_FOLDER
if not os.path.exists(RIF):
    raise OSError("RawInput folder not found. Please ensure it exists in your path and is named correctly.")

# Define a folder for where all processed data should go. PIF stands
# for processed input folder.
PIF = config.PROCESSED_INPUT_FOLDER
if not os.path.exists(PIF):
    print(f"{PIF} directory did not exist: Created {PIF}.")
    os.mkdir(PIF)

# Get a list of all the raw data folders. Parse non-folders by excluding
# anything with a file extension.
raw_data_folders = [folder for folder in os.listdir(RIF)
                    if len(folder.split(".")) == 1]
folder_counter = 1
reports_of_interest = [1,2,3,4,5] # maybe go into report.txt so we don't have to hard code this 
for folder_name in raw_data_folders:
    
    print(f"Starting folder {folder_counter}/{len(raw_data_folders)} {folder_name}")
    # Get all of the folders which have relevant Report01.csv
    #  files,
    # sorted alphanumerically.
    filenames = sorted([fn for fn in os.listdir(f"{RIF}/{folder_name}") if fn[-2:] == ".D"])

    acc_reports = []
    for report_num in reports_of_interest:
        # Make a DataFrame to store all of the Report01 files found within
        # a directory's sub-directories.
        folder_df = pd.DataFrame(columns=["FileName", "PeakNum", "RetTime", "Area", "AreaPerc"])

        print(f"\t~Processing REPORT0{report_num} raw data files ({len(filenames)}):")
        file_counter = 1
        for filename in filenames:
            print(f"\t\t-> {file_counter}/{len(filenames)} {filename}")
            # Catch error thrown if no report file is found
            try:
                # Get the contents of the old file
                with open(f"{RIF}/{folder_name}/{filename}/REPORT0{report_num}.csv", "r", encoding="utf-16") as file:
                    # Read lines then close to save memory
                    lines = file.readlines()

                # Store the file name, peak number, retention time, and
                # area in each line of the .csv file as a row in the
                # dataframe.
                for line in lines:
                    # Get relevant data from the line, exclude
                    # irrelevant data with _
                    peak_num, ret_time, _, _, area, _, area_perc = line.split(",")
                    # Add as the last row in folder_df
                    folder_df.loc[len(folder_df.index)] = (filename, peak_num, ret_time, area, area_perc)
                # Increment file counter
                file_counter += 1
                    
            except FileNotFoundError:
                # This is a non-fatal error, log it and skip this file
                logger.log(f"There was no report0{report_num} file for {RIF}/{folder_name}/{filename}.")
                # Increment file counter
                file_counter += 1
                continue  # skip to next in filenames

        # Save the folder's DataFrame if it is not empty.
        if len(folder_df.values) > 0:
            acc_reports.append(f"{PIF}/acc_rep{report_num}_{folder_name}.csv")
            folder_df.to_csv(f"{PIF}/acc_rep{report_num}_{folder_name}.csv", index=False, index_label=False)
            print("\t~Raw data processed")
        else:
            print("\t~Folder was empty")
        
    #combine the two accumulated reports and use their differences in area to determine ChemID
    print("\t~Combining REPORTs and identifying chemicals")
    wcid = wcidNEW.WCID(acc_reports, config.LAM_LIST)
    wcid.the_works().to_csv(f"{PIF}/acc_{folder_name}.csv", index=False, index_label=False)
    for i in range(len(acc_reports)):
        os.remove(acc_reports[i])
    print("\t~Folder completed and saved")
    print(f"saved as: {PIF}/acc_{folder_name}.csv")
    # Update the folder counter
    folder_counter += 1
