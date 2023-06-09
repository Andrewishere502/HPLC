from math import sqrt


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
    with open(f"MakeChemMeta/{filename}", "r",) as file:
        #named jacobs journal in honor of the person who ran the HPLC for us :)
        jacobs_journal = [float(n) for n in file.read().split(",")]  
    return jacobs_journal

def group_peaks(peaks):
    """"""
    unique_peaks = sorted(set(peaks))

    margin = 0.15
    buckets = [Bucket(peaks=[unique_peaks[0]])]  # Put first peak in bucket to start
    for peak in unique_peaks[1:]:  # skip the first peak
        if buckets[-1][-1] + margin >= peak:
            buckets[-1].add_peak(peak)
        else:
            buckets.append(Bucket([peak]))
    print(f"Number of buckets: {len(buckets)}")
    return buckets


chemmeta_file = open("MakeChemMeta/chemmeta.csv", "w")
chemmeta_file.write("BeginRetTime,EndRetTime,ChemicalID\n")

cards = group_peaks(load_peaks_txt("card.txt"))
pps =  group_peaks(load_peaks_txt("pp.txt"))

for bucket in cards:
    line = f"{bucket.min},{bucket.max},C{round(bucket.mean,1)}_Area\n"
    chemmeta_file.write(line)

for bucket in pps:
    line = f"{bucket.min},{bucket.max},PP{round(bucket.mean,1)}_Area\n"
    chemmeta_file.write(line)

chemmeta_file.close()
