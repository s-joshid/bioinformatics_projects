# -*- coding: utf-8 -*-
"""
Author: Dhruvi Joshi 

This script will download Total RNA seq data 
from human donor samples from the ENCODE database

Prior to running, please ensure you have 
https://www.encodeproject.org/report.tsv?type=Experiment&control_type!=*&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&biosample_ontology.classification=tissue&status=released
downloaded and saved as exp_report.tsv and 
https://www.encodeproject.org/metadata/?control_type%21=%2A&status=released&perturbed=f[…]iles.analyses.status=released&files.preferred_default=true
downloaded, this is already properly named

Additionally, please ensure you have a data_directory subdirectory prior

"""

#import statements
import pandas as pd 
import os
import subprocess
import time

#read in metadata 
metadata = pd.read_table("metadata.tsv")
print(metadata.head)

#read in experimental report
exp_report = pd.read_table("experiment_report.tsv",
                           header = 1,
                           usecols=['Accession', 'Life stage', 'Biosample age'])



#filter metadata 
final_metadata = metadata.loc[(metadata["Assay"] == "total RNA-seq") &
                    (metadata["Biosample organism"] == "Homo sapiens") &
                        (metadata["File format"] == "tsv") &
                        (metadata["Biosample type"] == "tissue") &
                        (metadata["Audit WARNING"].isna()) &
                        (metadata["Audit ERROR"].isna()) &
                        (metadata["Audit NOT_COMPLIANT"].isna())]
#renaming and casting to str
final_metadata = final_metadata.rename(columns = {"Experiment accession":"Accession"})
final_metadata = final_metadata.rename(columns = {"Biosample term name":"tissue_source"})
final_metadata["Accession"] = final_metadata["Accession"].astype(str)
exp_report["Accession"] = exp_report["Accession"].astype(str)

#joining dfs
final_metadata = final_metadata.merge(exp_report, on = "Accession", how = "left")

#formatting rows
final_metadata["Biosample age"] = final_metadata["Biosample age"].str.replace(" ", "_")
final_metadata["Biosample age"] = final_metadata["Biosample age"].str.replace(",", "_")
final_metadata["tissue_source"] = final_metadata["tissue_source"].str.replace(" ", "_")
final_metadata["tissue_source"] = final_metadata["tissue_source"].str.replace("'", "")
final_metadata["Life stage"] = final_metadata["Life stage"].str.replace(",", "_")

#export this df into a csv for refence 
final_metadata.to_csv("filtered_metadata.csv")

#data frame for URLs and download names 
final_metadata['file_name'] = final_metadata['Accession'] + '_' + final_metadata["tissue_source"] + "_" + final_metadata["Biosample age"]+ "_" + final_metadata["Life stage"] + ".tsv"

#create df to download URLS
downloads = final_metadata[["file_name", "File download URL"]]
downloads["command"] = "curl -L " + downloads["File download URL"] + " > ./data_directory/" + downloads["file_name"]
downloads = downloads.dropna(subset="command")


#setting up variables and directories 
dataframes = []
input_dir = "./data_directory/"
output_file = "./data_directory/combined_total_RNA_seq.tsv"

#downloading all files
#code will retry each file up to 3 times before skipping 
def execute_commands(df, retries = 3, delay = 5 ):
    for index, row in df.iterrows():
        command = row['command']
        success = False 
        for attempt in range(retries):
            
            try:
                #dowloading each file
                print(f"Executing command {index + 1}/{len(df)}")
                subprocess.run(command, shell=True, check=True)
                print(f"Command {index+1} completed successfully.")
                
                #processing each file and adding it to master dataframes list 
                file_path = os.path.join(input_dir, row['file_name'])
                curr_df = pd.read_csv(file_path, sep="\t", header=None)
                curr_df["file_name"] = row["file_name"]
                dataframes.append(curr_df)
                
                #deleting current file
                os.remove(file_path)
                success = True 
                break #exiting retry loop after successful download 
            #handeling misdownloaded files
            except subprocess.CalledProcessError as e: 
                print(f"Error executing command {index+1}") 
                print(f"Attempt {attempt + 1} failed: {e}")
                if attempt < retries - 1:
                    print(f"Retrying in {delay} seconds")
                    time.sleep(delay)
                else:
                    print(f"All {retries} attempts failed for command {index+1}.")
                
            if not success:
                print(f"skipping command {index +1} due to repeated failures.")
                continue
#calling method
execute_commands(downloads)


#concating dfs
combined_df = pd.concat(dataframes, ignore_index=True)

# Add column labels 
columns = ["gene_id", "transcript_id(s)", "length",	"effective_length " , "expected_count", "TPM", "FPKM",	"posterior_mean_count",	"posterior_standard_deviation_of_count", "pme_TPM", "pme_FPKM", "TPM_ci_lower_bound", "TPM_ci_upper_bound", "TPM_coefficient_of_quartile_variation",	"FPKM_ci_lower_bound", "FPKM_ci_upper_bound",	"FPKM_coefficient_of_quartile_variation", "file_name"][:combined_df.shape[1]]
combined_df.columns = columns

#add in age, life stage, and tissue data 
filtered_meta_sub = final_metadata[["Accession", "Biosample age", "tissue_source", "Life stage", "file_name"]]
combined_df = combined_df.merge(filtered_meta_sub, on = "file_name", how = "left")
# Save the combined DataFrame to a TSV file
combined_df.to_csv(output_file, sep="\t", index=False)
print(f"Combined file saved to: {output_file}")












