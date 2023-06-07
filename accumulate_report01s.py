"""Take the peak number, retention time, and area from each of the
Report01.csv files scattered throughout an HPLC Chromatogram directory,
and put that data in a file named after the directory the data was
found in.
"""
import os

from logger import Logger


# Create an object to log errors that occur
logger = Logger("", "log.txt", erase_on_init=True)

raw_input_folder = "RawInput"
processed_input_folder = "ProcessedInput"

# Get a list of all the raw data folders. Parse non-folders by excluding
# anything with a file extension.
raw_data_folders = [folder for folder in os.listdir(raw_input_folder)
                    if len(folder.split(".")) == 0]
for folder_name in raw_data_folders:

    # Get all of the folders which have relevant Report01.csv files,
    # sorted alphanumerically
    filenames = sorted([fn for fn in os.listdir(f"{raw_input_folder}/{folder_name}") if fn[-2:] == ".D"])

    # Make the folder for accumulated data
    accumulated_folder_name = "acc_" + folder_name
    try:
        os.mkdir(f"{processed_input_folder}/{accumulated_folder_name}")
    except FileExistsError:
        pass

    # Convert the Report01.csv files to only include peak, ret time
    # and area.
    for filename in filenames:
        # Catch the error thrown if no Report01.csv file is found
        try:
            # Get the contents of the old file
            with open(f"{raw_input_folder}/{folder_name}/{filename}/Report01.csv", "r", encoding="utf-16") as file:
                # NOTE: the Report01.csv files have no headers, so don't skip
                # first line.
                lines = file.readlines()
        except FileNotFoundError:
            # This is a non-fatal error, log it and skip this file
            logger.log(f"There was no Report01.csv file for {raw_input_folder}/{folder_name}/{filename}.")
            continue  # skip to next in filenames

        # Process the data and write it to a file named after the
        # lowest level parent folder the data was found in
        with open(f"{processed_input_folder}/{accumulated_folder_name}/{filename}.csv", "w") as file:
            header = "Peak,RetTime,Area\n"
            file.write(header)

            # Write only the peak number, retention time, and the area
            # under the peak to the new .csv file
            for line in lines:
                peak_num, ret_time, _, _, area, *_ = line.split(",")
                new_line = f"{peak_num},{ret_time},{area}\n"
                file.write(new_line)
