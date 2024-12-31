# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 22:40:51 2024

@author: Dhruvi Joshi

Script to perform pliminary data analysis for
total RNA seq data from human donors in the ENCODE database

please run get_total_rna_Seq.py prior to this script and 
ensure this script is located in the data_directory 
along with the combined_total_RNA_seq.tsv output from the previosu script
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# Read in total RNA seq tsv 
rna_Seq_df = pd.read_table("combined_total_RNA_seq.tsv")

print("running")
print(rna_Seq_df.head())


# Convert TPM column to numeric, coerce errors into NaN
rna_Seq_df["TPM"] = pd.to_numeric(rna_Seq_df["TPM"], errors="coerce")
print("TPM converted to numeric")
# Drop rows where TPM is nan
rna_Seq_df = rna_Seq_df.dropna(subset=["TPM"])
print("Nan values dropped")

#exploring the longest gene 
rna_Seq_df["length"] = pd.to_numeric(rna_Seq_df["length"], errors="coerce")
rna_Seq_df = rna_Seq_df.dropna(subset=["length"])
largest_value_index = rna_Seq_df['length'].idxmax()
largest_row = rna_Seq_df.loc[largest_value_index]

#break df into the 3 stages

emb_df  = rna_Seq_df[rna_Seq_df["Life stage"] == "embryonic"]
print("emb_df made")

child_df = rna_Seq_df[rna_Seq_df["Life stage"] == "child"]
print("child_df made")

adult_df = rna_Seq_df[rna_Seq_df["Life stage"] == "adult"]
print("adult_df made")

# Define a function to get most commonly expressed genes
def get_most_common_genes(df):
    non_zero_counts = (
        df[df["TPM"] > 0]
        .groupby("gene_id")["TPM"]
        .count()
        .reset_index(name="non_zero_count")
    )
    max_count = non_zero_counts["non_zero_count"].max()
    most_common_genes = non_zero_counts[non_zero_counts["non_zero_count"] == max_count]
    return set(most_common_genes["gene_id"])

# Get most common genes for each life stage
most_common_emb_genes = get_most_common_genes(emb_df)
most_common_child_genes = get_most_common_genes(child_df)
most_common_adult_genes = get_most_common_genes(adult_df)

# Find unique most common genes
unique_emb_common_genes = most_common_emb_genes - (most_common_child_genes | most_common_adult_genes)
unique_child_common_genes = most_common_child_genes - (most_common_emb_genes | most_common_adult_genes)
unique_adult_common_genes = most_common_adult_genes - (most_common_emb_genes | most_common_child_genes)

# Print results
print(f"Number of unique most commonly expressed genes in embryonic stage: {len(unique_emb_common_genes)}")
print(f"Unique most commonly expressed genes in embryonic stage: {unique_emb_common_genes}")

print(f"Number of unique most commonly expressed genes in child stage: {len(unique_child_common_genes)}")
print(f"Unique most commonly expressed genes in child stage: {unique_child_common_genes}")

print(f"Number of unique most commonly expressed genes in adult stage: {len(unique_adult_common_genes)}")
print(f"Unique most commonly expressed genes in adult stage: {unique_adult_common_genes}")

#exploring p53
df_p53_ad = adult_df[adult_df["gene_id"].str.contains("ENSG00000141510", na = False)]
df_p53_em = emb_df[emb_df["gene_id"].str.contains("ENSG00000141510", na = False)]
df_p53_ch = child_df[child_df["gene_id"].str.contains("ENSG00000141510", na = False)]

#reconfiguring adult into youner older

df_p53_ad["Biosample age"] = pd.to_numeric(df_p53_ad["Biosample age"].str.split("_").str[0])



# Calculate average TPM
emb_avg_tpm = df_p53_em["TPM"].sum() / len(df_p53_em)
child_avg_tpm = df_p53_ch["TPM"].sum() / len(df_p53_ch)

# Split adults into two groups
adult_under_50 = df_p53_ad[df_p53_ad["Biosample age"].astype(float) <= 50]
adult_over_50 = df_p53_ad[df_p53_ad["Biosample age"].astype(float) > 50]

adult_under_50_avg_tpm = adult_under_50["TPM"].sum() / len(adult_under_50)
adult_over_50_avg_tpm = adult_over_50["TPM"].sum() / len(adult_over_50)

# Prepare data for plotting
groups = ["Embryonic", "Child", "Adult ≤50", "Adult >50"]
avg_tpm_values = [emb_avg_tpm, child_avg_tpm, adult_under_50_avg_tpm, adult_over_50_avg_tpm]

# Plot the bar graph
plt.figure(figsize=(8, 6))
plt.bar(groups, avg_tpm_values, color=["blue", "green", "orange", "red"])
plt.title("Average TPM of P53 by Age Group")
plt.ylabel("Average TPM")
plt.xlabel("Life Stage / Age Group")
plt.xticks(rotation=15)
plt.tight_layout()
plt.show()

#plotting p53 over diff tissue sources

# Calculate average TPM by tissue source for each group
avg_tpm_em = df_p53_em.groupby("tissue_source")["TPM"].mean()
avg_tpm_ch = df_p53_ch.groupby("tissue_source")["TPM"].mean()
avg_tpm_ad_under_50 = adult_under_50.groupby("tissue_source")["TPM"].mean()
avg_tpm_ad_over_50 = adult_over_50.groupby("tissue_source")["TPM"].mean()

# Combine results into a single DataFrame
avg_tpm_df = pd.DataFrame({
    "Embryonic": avg_tpm_em,
    "Child": avg_tpm_ch,
    "Adult ≤50": avg_tpm_ad_under_50,
    "Adult >50": avg_tpm_ad_over_50
}).fillna(0)  # Fill missing values with 0

# Plot grouped bar chart
tissue_sources = avg_tpm_df.index
age_groups = avg_tpm_df.columns
bar_width = 0.2
x = np.arange(len(tissue_sources))

plt.figure(figsize=(12, 6))

# Create bars for each age group
for i, age_group in enumerate(age_groups):
    plt.bar(
        x + i * bar_width,
        avg_tpm_df[age_group],
        width=bar_width,
        label=age_group
    )

# Formatting
plt.xticks(x + bar_width * 1.5, tissue_sources, rotation=90, ha="right")
plt.title("Average TPM for p53 by Tissue Source and Age Group")
plt.xlabel("Tissue Source")
plt.ylabel("Average TPM")
plt.legend(title="Age Group")
plt.tight_layout()

plt.show()

