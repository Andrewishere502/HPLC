# Absolute or relative path to the txt file containing cardenolide peak
# retention times.
CARD_TXT = "card.txt"

# Absolute or relative path to the txt file containing phenylpropanoid
# peak retention times.
PP_TXT = "pp.txt"

#The wavelength used by the HPLC to emphasize cardenolides
LAM_CARD = 350

#The wavelength used by the HPLC to emphasize phenylpropanoids
LAM_PP = 330

# Margin for buckets
BUCKET_MARGIN = 0.1

# Whether or not chemmeta should delete ambiguous labels
DELETE_AMBIG = True

# Max and minimum range for a compound to be recognized
MAX_RANGE = 1.3
MIN_RANGE = 0.0

# The absolute or relative path to the .csv file that stores the
# chemical meta data. This must exist before running
# generate_dataset.py. Must be a csv file.
CHEM_META_FILE = f"chemmeta_{BUCKET_MARGIN}_Ambig{not DELETE_AMBIG}.csv"

# The absolute or relative path to the .csv file that stores the
# sample (super) meta data. This must exist before running
# generate_dataset.py. Must be a csv file.
SUPER_META_FILE = "supermeta.csv"

# The absolute or relative path to the directory where raw input is
# stored. This must exist before running accumulate_report01s.py.
# Do not include a '/' at the end.
RAW_INPUT_FOLDER = "RawInput"

# The absolute or relative path to the directory where processed input
# will be stored. If it doesn't exist, it will be created for you.
# Do not include a '/' at the end.
PROCESSED_INPUT_FOLDER = "ProcessedInput"

# The absolute or relative path to the directory where the final output
# will be stored. If it doesn't exist, it will be created for you.
# Do not include a '/' at the end.
OUTPUT_FOLDER = "Output"
