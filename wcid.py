import pandas as pd
import config

class WCID:
    def __init__(self, card_filename, pp_filename):
        self.card_filename = card_filename
        self.pp_filename = pp_filename
        self.card_df = pd.read_csv(f"{self.card_filename}")
        self.pp_df = pd.read_csv(f"{self.pp_filename}")

    def split_clean_dfs(self):
        ###Each sample has two dataframes (one to emphasize pp, and the other to emphasize card).
        ###Each dataframe will have a slightly different rettime, and a significantly different area for each peak
        ###If the pp-emphasis area is higher, that peak is a pp, same is true for a card
        ###but to compare the areas in each df, we have to figure out what rettimes are talking about the same peak
        #I chose to round them to one decimal; thus if they were essentially equal, we have now made them equal
        self.card_df["RoundRet"] = [round(peak, 1) for peak in self.card_df["RetTime"]]
        self.pp_df["RoundRet"] = [round(peak, 1) for peak in self.pp_df["RetTime"]]

        ###split that df into a df for each sample
        #this is makes a grouped df that clumps all of the rows with the same file name 
        card_df_grouped = self.card_df.groupby(self.card_df.FileName)
        #this makes a list (array?) of each unique file name. This tells us how many samples are in this folder.
        sample_names = card_df_grouped.FileName.unique()
        #Finally, we make a list of dfs. Each df reprents one HPLC graph. It's just a bunch of peaks for one sample.
        card_sample_dfs = []
        for sample_name in sample_names:
            card_sample_dfs.append(card_df_grouped.get_group(str(sample_name)[2:-2]))

        #see above, we just do the same thing of pps. I should probs make this a method :/
        pp_df_grouped = self.pp_df.groupby(self.pp_df.FileName)
        sample_names = pp_df_grouped.FileName.unique()
        pp_sample_dfs = []
        for sample_name in sample_names:
            pp_sample_dfs.append(pp_df_grouped.get_group(str(sample_name)[2:-2]))


        #Check to see if the pp and card dfs have been split into the same number of samples
        if len(pp_sample_dfs) != len(card_sample_dfs):
            print(f"Error: there are {len(card_sample_dfs)} samples in the card df and {len(pp_sample_dfs)} in the pp df.")
    
        ###Now that we have the pp and card dfs split into individual samples, 
        ###we can combine the pp and card dfs for each sample,
        ###but first we have to clean up the dfs

        #for every sample in this experiment
        for i in range(len(pp_sample_dfs)):

            ###The cleaning and organization starts here. We need to make these two dfs similar enough to compare
            ###it is important that we do this after we split up are orignial two dfs into two groups of sample dfs
            ###we wouldn't want stuff from from other samples to affect how this sample is cleaned

            #Make sure we are dealing with the pp and card dfs that correspond to the same sample
            if pp_sample_dfs[i]["FileName"].values[0] != card_sample_dfs[i]["FileName"].values[0]:
                print("Error: File names do not match.")
            
            #Get rid of any duplicates with one df
            #Duplicates are peaks that were rounded to the same value. 
            #They are close enough, that they should be treated as one, so we drop the duplicate
            pp_sample_dfs[i] = pp_sample_dfs[i].drop_duplicates(subset="RoundRet")
            card_sample_dfs[i] = card_sample_dfs[i].drop_duplicates(subset="RoundRet")

            #At this point, every peak has a unique rounded ret time, so we breifly use this as our index 
            pp_sample_dfs[i] = pp_sample_dfs[i].set_index("RoundRet")
            card_sample_dfs[i] = card_sample_dfs[i].set_index("RoundRet")

        return card_sample_dfs, pp_sample_dfs

    def combine_dfs(self, card_sample_dfs, pp_sample_dfs):
        #this is the list where we will store our finished, combined dfs
        combo_dfs = []
        #for every sample in this experiment
        for i in range(len(pp_sample_dfs)):
            #an empty data frame for our combined data to land
            combo_df = pd.DataFrame(
                {
                "FileName": [],
                "RoundRet": [],
                "330Area": [],
                "350Area": []
                }
            ) 

            #find out what df is longer. They should be about the same size, 
            #but sometimes there are more peaks in one because they are integrated differently at different wavelengths
            if len(pp_sample_dfs[i]) > len(card_sample_dfs[i]):
                longer_df = pp_sample_dfs[i]
                shorter_df = card_sample_dfs[i]


                #Now that we are organized, we can combine our dfs


                #longer_df.index is a list (array?) of RoundRet. 
                #RoundRet was set to the index so it would be easier to comare these dfs
                for peak in longer_df.index:
                    #This means a peak has to be in BOTH dfs for it to be considered.
                    if peak in shorter_df.index:
                        #if the peak is in both, add its data to the combo_df
                        new_row = {
                            "FileName": longer_df["FileName"].loc[peak],
                            "RoundRet": peak,
                            "330Area": longer_df["Area"].loc[peak],
                            "350Area": shorter_df["Area"].loc[peak]
                            }
                        combo_df.loc[len(combo_df)] = new_row

            #see above, this probably should have been a method :/            
            else:
                longer_df = card_sample_dfs[i]
                shorter_df = pp_sample_dfs[i]

                for peak in longer_df.index:
                    if peak in shorter_df.index:
                        new_row = {
                            "FileName": longer_df["FileName"].loc[peak],
                            "RoundRet": peak,
                            "330Area": shorter_df["Area"].loc[peak],
                            "350Area": longer_df["Area"].loc[peak]
                            }
                        combo_df.loc[len(combo_df)] = new_row

            #once we have complete a combo_df, we add it to this list
            combo_dfs.append(combo_df)
        return combo_dfs

    def assign_ids(self, combo_dfs):
        #every sample has a combo_df, so this essentially says "for every sample"
        for df in combo_dfs:
            #for every peak in that sample
            for i in range(len(df)):

                #is the pp-emphasized peak bigger than the c-emphasized peak? If so, it's a pp
                if df.loc[i, "330Area"] > df.loc[i, "350Area"]:
                    #this creates and populates the "ID" column
                    df.loc[i, "ID"] = "PP"
                #is the c-emphasized peak bigger than the pp-emphasized peak? If so, it's a c        
                if df.loc[i, "330Area"] < df.loc[i, "350Area"]:
                    df.loc[i, "ID"] = "C"
                #there is no way a peak would have the same area if it was truly a pp or a c.
                #if it has the same area, it's neither, so let's drop it.
                ###This code could be improved by dropping peaks that have essentially to same area. 
                if df.loc[i, "330Area"] == df.loc[i, "350Area"]:
                    df.drop(i)
                    print(f"Error: Peak at {df.loc[i, 'RoundRet']} has same area at both wavelengths. Peak has been dropped")
        return combo_dfs
            
    def concat_samples(self, combo_dfs):
        #each of our combo_dfs will end up in here so they they can be exported to a csv
        printing_df = pd.DataFrame(columns=["FileName", "RoundRet", "330Area", "350Area"])

        #every sample has a combo_df, so this essentially says "for every sample"
        for df in combo_dfs:
            #once every peak have been named, add that sample to the printing_df
            printing_df = pd.concat([printing_df, df])
        return printing_df

        #printing_df.to_csv("cpp.csv", index=False, index_label=False)
        #print("cpp.csv was created")

    def the_works(self):
        card_sample_dfs, pp_sample_dfs = self.split_clean_dfs()
        combo_dfs = self.combine_dfs(card_sample_dfs, pp_sample_dfs)
        combo_dfs = self.assign_ids(combo_dfs)
        return self.concat_samples(combo_dfs)

#WCID(pd.read_csv(f"{config.PROCESSED_INPUT_FOLDER}/acc_rep3_New2022-06-2910-10-34.csv"), pd.read_csv(f"{config.PROCESSED_INPUT_FOLDER}/acc_rep4_New2022-06-2910-10-34.csv")).the_works().to_csv("zz.csv", index=False, index_label=False)



            

