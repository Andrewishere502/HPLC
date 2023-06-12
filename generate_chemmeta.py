from math import sqrt

import pandas as pd

import config

class Bucket:
    def __init__(self, peaks=[]):
        self.peaks = peaks
        return

    @property
    def mean(self):
        """Return the mean for self.peaks."""
        return sum(self.peaks) / len(self.peaks)

    @property
    def std(self):
        """Calculate the population standard deviation for self.peaks."""
        u = self.mean
        N = len(self.peaks)
        std = sqrt(sum([(x - u) ** 2 for x in self.peaks]) / N)
        return std

    @property
    def min(self):
        return min(self.peaks)

    @property
    def max(self):
        return max(self.peaks)

    def add_peak(self, peak):
        """Add a peak to the population."""
        self.peaks.append(peak)
        return

    def info(self):
        """Return a string with info on size of self.peaks, and the
        mean and standard deviation, and min and max
        """
        return f"size: {len(self.peaks)}\nmean: {self.mean}\nStd: {self.std}\nMinimum: {self.min}\nMaximum: {self.max}"

    def is_outlier(self, peak):
        """Return True if peak is more than 3 standard deviations from
        the mean.
        """
        out_of_lower_bound = peak <= self.mean - 3 * self.std
        out_of_upper_bound = peak >= self.mean + 3 * self.std
        return out_of_lower_bound or out_of_upper_bound

    def outliers(self):
        """Return the outliers for this bucket. If there is only one
        peak in this bucket, return an empty list.
        """
        if len(self.peaks) > 1:
            return [p for p in self.peaks if self.is_outlier(p)]
        else:
            return []  # self.std is 0 if there is one peak

    def __getitem__(self, i):
        return self.peaks[i]
    

def load_peaks_txt(filename):
    """Return a sorted float list of the peaks in a given text file."""
    with open(f"{filename}", "r",) as file:
        #named jacobs journal in honor of the person who ran the HPLC for us :)
        jacobs_journal = [float(n) for n in file.read().split(",")]  
    return jacobs_journal

def group_peaks(peaks, margin):
    """Iteratively group peaks that are within margin of each other until
    no more peaks are within margin of each other.
    """
    # Get all unique peaks and sorted them
    unique_peaks = sorted(set(peaks))
    
    # Group the peaks into buckets
    buckets = [Bucket(peaks=[unique_peaks[0]])]  # Put first peak in bucket to start
    for peak in unique_peaks[1:]:  # skip the first peak
        if buckets[-1][-1] + margin >= peak:
            buckets[-1].add_peak(peak)
        else:
            buckets.append(Bucket([peak]))
    return buckets

def name_chem_range(id_prefix, mean):
    """Return the name for a chemical range's mean."""
    return f"{id_prefix}{round(mean, 1)}_Area"

def write_buckets(buckets, chemmeta_file, id_prefix=""):
    for bucket in buckets:
        line = f"{bucket.min},{bucket.max},{name_chem_range(id_prefix, bucket.mean)}\n"
        chemmeta_file.write(line)
    return

def add_chem(df, name, start, end):
    """Add a single new chem to the df."""
    df.loc[name] = round(start, 1), round(end, 1)
    return

def change_chem_times(df, name1, name2, start1, end1, start2, end2):
    """Replace old start and end retention times with new values."""
    add_chem(df, name1, start1, end1)
    add_chem(df, name2, start2, end2)
    return

def combine_names(name1, name2):
    """Return the combined name for two chem ids."""
    name = name1.split("_")[0] + name2.split("_")[0] + "_Area"
    return name

def merge_chem_times(df, cpp1, cpp2):
    """Add a new row to the df with the combined name of the
    cardenolide and phenylpropanoid.
    """
    add_chem(df, name_chem_range("CPP", (cpp1 + cpp2) / 2), cpp1, cpp2)
    return


# Set up the chemical meta data file, erasing its contents if it
# already exists.
with open("chemmeta.csv", "w") as chemmeta_file:
    chemmeta_file.write("BeginRetTime,EndRetTime,ChemicalID\n")

    # Get the cardenolide specific buckets
    card_buckets = group_peaks(load_peaks_txt(config.CARD_TXT), config.CARD_MARGIN)
    # Write the card buckets to chemmeta
    write_buckets(card_buckets, chemmeta_file, id_prefix="C")
    print(f"Number of cardenolide buckets: {len(card_buckets)}")

    # Get the phenylpropanoid specific buckets
    pp_buckets =  group_peaks(load_peaks_txt(config.PP_TXT), config.PP_MARGIN)
    # Write teh card buckets to chemmeta
    write_buckets(pp_buckets, chemmeta_file, id_prefix="PP")
    print(f"Number of phenylpropanoid buckets: {len(pp_buckets)}")


# Get all unique overlapping card and pp pairs
with open("chemmeta.csv", "r") as chemmeta_file:
    lines = chemmeta_file.readlines()
hits = set()
header = lines.pop(0)
for line in lines:
    line = line.split(",")
    begin, end, chem_id = line
    #if chem_id[0] == 'P':
        #break
    begin = float(begin)
    end = float(end)
    for line in lines:
        line = line.split(",")
        begin2, end2, chem_id2 = line
        if chem_id == chem_id2:
            continue
        begin2 = float(begin2)
        end2 = float(end2)
        if ((begin >= float(begin2) and begin <= float(end2))
            or (end >= float(begin2) and end <= float(end2))):
            hits.add(tuple(sorted((chem_id[:-1], chem_id2[:-1]))))
            #hits.add((chem_id, chem_id2))


# Fix all overlaps
delete_me = set()
df = pd.read_csv("chemmeta.csv", index_col=2)
for c_name, pp_name in hits:
    c1, c2 = df.loc[c_name]
    p1, p2 = df.loc[pp_name]

    # Start-End overlap
    if c1 < p1 and c2 < p2:
        # Modify C and PP
        change_chem_times(df, c_name, pp_name, c1, p1, c2, p2)
        # Add a CPP
        merge_chem_times(df, p1 + 0.1, c2 - 0.1)
    elif p1 < c1 and p2 < c2:
        # Modify C and PP
        change_chem_times(df, c_name, pp_name, p2, c2, p1, c1)
        # Add a CPP
        merge_chem_times(df, c1 + 0.1, p2 - 0.1)

    # Contains
    elif c1 < p1 and c2 > p2:
        # Add two Cs
        c_name1 = name_chem_range("C", (c1 + p1) / 2)
        c_name2 = name_chem_range("C", (p2 + c2) / 2)
        change_chem_times(df, c_name1, c_name2, c1, p1, p2, c2)
        # Remove old C
        delete_me.add(c_name)
        # Add a CPP
        merge_chem_times(df, p1 + 0.1, p2 - 0.1)
    elif c1 > p1 and c2 < p2:
        # Add two PPs
        p_name1 = name_chem_range("PP", (p1 + c1) / 2)
        p_name2 = name_chem_range("PP", (c2 + p2) / 2)
        change_chem_times(df, p_name1, p_name2, p1, c1, c2, p2)
        # Remove old PP
        delete_me.add(pp_name)
        # Add a CPP
        merge_chem_times(df, c1 + 0.1, c2 - 0.1)

    # Start-End same
    elif c2 == p1 and c1 < p2:
        # Just shift the p range over by 0.1
        change_chem_times(df, c_name, pp_name, c1, c2, p1 + 0.1, p2)
    elif p2 == c1 and p1 < c2:
        # Just shift the c range over by 0.1
        change_chem_times(df, c_name, pp_name, c1 + 0.1, c2, p1, p2)
    
    # Start-Start same
    elif c1 == p1 and c2 < p2:
        # Add a new PP range
        add_chem(df, name_chem_range("PP", (c2 + p2) / 2), c2, p2)
        # Delete old C range
        delete_me.add(c_name)
        # Add a CPP
        merge_chem_times(df, c1, c2 - 0.1)
    elif c1 == p1 and c2 > p2:
        # Add a new C range
        add_chem(df, name_chem_range("C", (p2 + c2) / 2), p2, c2)
        # Delete old PP range
        delete_me.add(pp_name)
        # Add a CPP
        merge_chem_times(df, p1, p2 - 0.1)

    # End-End same
    elif c2 == p2 and c1 < p1:
        # Add a new C range
        add_chem(df, name_chem_range("C", (c1 + p1) / 2), c1, p1)
        # Delete old PP range
        delete_me.add(pp_name)
        # Add cpp
        merge_chem_times(df, p1 + 0.1, p2)

    elif c2 == p2 and c1 > p1:
        # Add a new C range
        add_chem(df, name_chem_range("PP", (p1 + c1) / 2), p1, c1)
        # Delete old PP range
        delete_me.add(c_name)
        # Add cpp
        merge_chem_times(df, c1 + 0.1, c2)

    # Start-Start and End-End same
    elif c1 == p1 and c2 == p2:
        # Remove the old C range
        delete_me.add(c_name)
        # Remove the old PP range
        delete_me.add(pp_name)
        # Merge
        merge_chem_times(df, c1, p1)

    # Safety catch-all
    else:
        raise ValueError(f"Oh crap... {c1} {c2} and {p1} {p2} should overlap but didn't?")


# Drop all necessary rows
df.drop(labels=list(delete_me), axis=0, inplace=True)

# Save the csv
df.to_csv("chemmeta.csv")
