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



jacobs_journal = [12.28,16.80,17.10,18.00,19.62,22.70,14.20,16.80,17.09,
                  18.00,19.40,19.60,22.70,17.08,22.70,17.10,19.40,22.80,
                  14.25,17.12,22.76,14.24,17.11,12.28,14.24,14.49,17.08,
                  17.99,19.13,19.40,19.60,22.71,14.12,16.97,14.13,16.98,
                  14.12,16.97,12.26,14.22,17.07,22.69,26.25,12.26,17.06,
                  17.97,19.59,20.45,12.25,16.79,17.06,17.96,19.59,22.68,
                  14.18,17.05,14.16,17.02,14.13,14.39,17.00,14.14,17.00,
                  14.10,16.96,14.08,16.92,16.01,22.50,14.06,16.93,18.96,
                  19.73,22.53,14.10,16.96,14.11,16.97,18.99,22.57,12.20,
                  16.97,17.87,19.49,22.57,26.11,16.95,22.56,16.93,19.44,
                  22.52,12.15,14.08,16.93,12.15,16.92,17.82,19.44,22.51,
                  26.05,16.92,17.82,19.45,22.53,12.15,16.92,17.82,19.44,
                  22.52,14.07,16.93,16.70,16.96,17.86,19.48]

# Sort unique_peaks just for output-readability
unique_peaks = sorted(set(jacobs_journal))

margin = 0.15
buckets = [Bucket(peaks=[unique_peaks[0]])]  # Put first peak in bucket to start
for peak in unique_peaks[1:]:  # skip the first peak
    if buckets[-1][-1] + margin >= peak:
        buckets[-1].add_peak(peak)
    else:
        buckets.append(Bucket([peak]))
print(f"Number of buckets: {len(buckets)}")


chemmeta_file = open("MakeChemMeta/chemmeta.csv", "w")
chemmeta_file.write("BeginRetTime,EndRetTime,ChemicalID")
for bucket in buckets:
    line = f"{bucket.min},{bucket.max},PP{round(bucket.mean,1)}_Area\n"
    chemmeta_file.write(line)

chemmeta_file.close()
