# Absolute or relative path to the txt file containing cardenolide peak
# retention times.
CARD_TXT = "card.txt"

# Absolute or relative path to the txt file containing phenylpropanoid
# peak retention times.
PP_TXT = "pp.txt"

# Margin for Cardenolide buckets
CARD_MARGIN = 0.2

# Margin for Phenylpropanoid buckets
PP_MARGIN = 0.25

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
