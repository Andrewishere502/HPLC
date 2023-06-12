# Absolute or relative path to the txt file containing cardenolide peak
# retention times.
CARD_TXT = "/Volumes/OlsonLab/card.txt"

# Absolute or relative path to the txt file containing phenylpropanoid
# peak retention times.
PP_TXT = "/Volumes/OlsonLab/pp.txt"

# Margin for Cardenolide buckets
CARD_MARGIN = 0.20

# Margin for Phenylpropanoid buckets
PP_MARGIN = 0.25

# The absolute or relative path to the .csv file that stores the
# chemical meta data. This must exist before running
# generate_dataset.py. Must be a csv file.
CHEM_META_FILE = "chemmeta.csv"

# The absolute or relative path to the .csv file that stores the
# sample (super) meta data. This must exist before running
# generate_dataset.py. Must be a csv file.
SUPER_META_FILE = "/Volumes/OlsonLab/supermeta.csv"

# The absolute or relative path to the directory where raw input is
# stored. This must exist before running accumulate_report01s.py.
# Do not include a '/' at the end.
RAW_INPUT_FOLDER = "/Volumes/OlsonLab/RawInput"

# The absolute or relative path to the directory where processed input
# will be stored. If it doesn't exist, it will be created for you.
# Do not include a '/' at the end.
PROCESSED_INPUT_FOLDER = "/Volumes/OlsonLab/ProcessedInput"

# The absolute or relative path to the directory where the final output
# will be stored. If it doesn't exist, it will be created for you.
# Do not include a '/' at the end.
OUTPUT_FOLDER = "/Volumes/OlsonLab/Output"
