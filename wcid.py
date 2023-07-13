import pandas as pd
import config

class WCID:
    def __init__(self, acc_reports, lam_list):
        self.lam_list = lam_list
        self.lam_dic = {}
        for i in range(len(acc_reports)):
            self.lam_dic[lam_list[i]] = pd.read_csv(acc_reports[i])

    def split_by_sample(self):
        self.samples_dic = {}
        for lam in self.lam_list:
            ###Each sample has two dataframes (one to emphasize pp, and the other to emphasize card).
            ###Each dataframe will have a slightly different rettime, and a significantly different area for each peak
            ###If the pp-emphasis area is higher, that peak is a pp, same is true for a card
            ###but to compare the areas in each df, we have to figure out what rettimes are talking about the same peak
            #I chose to round them to one decimal; thus if they were essentially equal, we have now made them equal
            self.lam_dic[lam]["RoundRet"] = [round(peak, 1) for peak in self.lam_dic[lam]["RetTime"]]

            ###split that df into a df for each sample
            #this is makes a grouped df that clumps all of the rows with the same file name 
            df_grouped = self.lam_dic[lam].groupby(self.lam_dic[lam].FileName)
            #this makes a list (array?) of each unique file name. This tells us how many samples are in this folder.
            sample_names = df_grouped.FileName.unique()
            #Finally, we make a list of dfs. Each df reprents one HPLC graph. It's just a bunch of peaks for one sample.
            sample_list = []
            for sample_name in sample_names:
                sample_list.append(df_grouped.get_group(str(sample_name)[2:-2]))
            self.samples_dic[lam] = sample_list

        #Check to see if the pp and card dfs have been split into the same number of samples
        for samples in self.samples_dic.values():
            if len(list(self.samples_dic.values())[0]) != len(samples):
                print(f"Error: Nonstandard number of samples across wavelength. Literally don't know how that's possible.")
        pseudo = self.samples_dic
        #return pseudo

        ###Now that we have the pp and card dfs split into individual samples, 
        ###we can combine the pp and card dfs for each sample,
        ###but first we have to clean up the dfs
    
    def clean_dfs(self):
        for samples in list(self.samples_dic.values()):
            #for every sample in this experiment
            for i in range(len(samples)):

                ###The cleaning and organization starts here. We need to make these two dfs similar enough to compare
                ###it is important that we do this after we split up are orignial two dfs into two groups of sample dfs
                ###we wouldn't want stuff from from other samples to affect how this sample is cleaned
                
                #Get rid of any duplicates with one df
                #Duplicates are peaks that were rounded to the same value. 
                #They are close enough, that they should be treated as one, so we drop the duplicate
                samples[i] = samples[i].drop_duplicates(subset="RoundRet")

                #At this point, every peak has a unique rounded ret time, so we breifly use this as our index 
                samples[i] = samples[i].set_index("RoundRet")
        pseudo = self.samples_dic
        #return pseudo

    def combine_lams(self):
        #this is the list where we will store our finished, combined dfs
        self.combo_dfs = []
        #for every sample in this experiment
        for i in range(len(self.samples_dic[330])):
            #an empty data frame for our combined data to land
            combo_df = pd.DataFrame({"FileName": [], "RoundRet": []})
            for lam in self.lam_list:
                combo_df[f"{lam}Area"] = []
                combo_df[f"{lam}AreaPerc"] = []
            

            #find out what df is longest. They should be about the same size, 
            #but sometimes there are more peaks in one because they are integrated differently at different wavelengths
            longlam = 0
            maxlen = 0
            for lam in self.lam_list:
                if len(self.samples_dic[lam][i]) > maxlen:
                    maxlen = len(self.samples_dic[lam][i])
                    longlam = lam


            #Now that we are organized, we can combine our dfs


            #longer_df.index is a list (array?) of RoundRet. 
            #RoundRet was set to the index so it would be easier to comare these dfs
            for peak in self.samples_dic[longlam][i].index:
                #This means a peak has to be in ALL dfs for it to be considered.
                peakInAll = True
                for lam in self.lam_list:
                    if peak not in self.samples_dic[lam][i].index:
                        peakInAll = False

                if (peakInAll):
                    #if the peak is in all, add its data to the combo_df
                    new_row = {
                        "FileName": self.samples_dic[longlam][i]["FileName"].loc[peak],
                        "RoundRet": float(peak)
                    }
                    for lam in self.lam_list: 
                        new_row[f"{lam}Area"] = self.samples_dic[lam][i]["Area"].loc[peak]
                        new_row[f"{lam}AreaPerc"] = self.samples_dic[lam][i]["AreaPerc"].loc[peak]
            
                    combo_df.loc[len(combo_df)] = new_row

            #once we have complete a combo_df, we add it to this list
            #combo_df = combo_df.loc[(combo_df["330Area"] > 800) & (combo_df["350Area"] > 800)]
            self.combo_dfs.append(combo_df)
        pseudo = self.combo_dfs
        #return pseudo

    def assign_ids(self):
        #every sample has a combo_df, so this essentially says "for every sample"
        for df in self.combo_dfs:
            #for every peak in that sample
            for i in range(len(df)):
                
            
                #is the pp-emphasized peak bigger than the c-emphasized peak? If so, it's a pp

                if ((df.loc[i, "330Area"] > df.loc[i, "350Area"])
                    & (df.loc[i, "240Area"] > df.loc[i, "250Area"])
                    ):
                    #this creates and populates the "ID" column
                    df.loc[i, "ID"] = "PP"
                #is the c-emphasized peak bigger than the pp-emphasized peak? If so, it's a c        
                if ((df.loc[i, "330Area"] < df.loc[i, "350Area"])
                    & (df.loc[i, "240Area"] < df.loc[i, "250Area"])
                    ):
                    df.loc[i, "ID"] = "C"
                #there is no way a peak would have the same area if it was truly a pp or a c.
                #if it has the same area, it's neither, so let's drop it.
                ###This code could be improved by dropping peaks that have essentially to same area. 
                # if df.loc[i, "330Area"] == df.loc[i, "350Area"]:
                #     df.drop(i)
                #     print(f"Error: Peak at {df.loc[i, 'RoundRet']} has same area at both wavelengths. Peak has been dropped")
        pseudo = self.combo_dfs
        #return pseudo
            
    def concat_samples(self):
        #each of our combo_dfs will end up in here so they they can be exported to a csv
        #printing_df = pd.DataFrame(columns=["FileName", "RoundRet", "330Area", "350Area"])
        self.printing_df = pd.DataFrame(columns=self.combo_dfs[0].columns)

        #every sample has a combo_df, so this essentially says "for every sample"
        for df in self.combo_dfs:
            #once every peak have been named, add that sample to the printing_df
            self.printing_df = pd.concat([self.printing_df, df])
        self.printing_df = self.printing_df.loc[(self.printing_df["330Area"] > 800) & (self.printing_df["350Area"] > 800)]
        self.printing_df = self.printing_df.loc[self.printing_df["RoundRet"] > 4]

        #printing_df.to_csv("cpp.csv", index=False, index_label=False)
        #print("cpp.csv was created")

    def the_works(self):
        self.split_by_sample()
        self.clean_dfs()
        self.combine_lams()
        self.assign_ids()
        self.concat_samples()
        return self.printing_df
    
    def group_ids(self):
        pp = list(self.printing_df[self.printing_df["ID"] == "PP"].RoundRet.values)
        card = list(self.printing_df[self.printing_df["ID"] == "C"].RoundRet.values)
        diff = max([len(pp), len(card)]) - min([len(pp), len(card)])
        if len(pp) < len(card):
            for i in range(diff):
                pp.append("")
        if len(card) < len(pp):
            for i in range(diff):
                card.append("")
        df = pd.DataFrame()
        df["Phenylpropanoids"] = pp
        df["Cardenolides"] = card
        return df
