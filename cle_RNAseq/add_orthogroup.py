# Name: add orthogroup
# Author: EY
# Date: May 15 2024
# Version: Python 3.10
# Description: add Arabidopsis orthogroup to DESeq output

import os
import pandas as pd

def process_files(directory, output_filename):
    """function will read in all .csv files in a directory, add arabipopsis orthogroups, and output as a single excel
    file with different sheets"""
    dataframes = {}
    orthogroups = pd.read_csv("Aug30_Orthogroups_longest_iso.csv", sep=",")
    column_to_search = 'Helianthus'
    values_lists = []

    for filename in os.listdir(directory):
        if filename.endswith('.csv'):
            filepath = os.path.join(directory, filename)
            df = pd.read_csv(filepath)
            df_subset = df[df['padj'] < 0.05]
            dataframes[filename] = df_subset

            values_list = []
            for index, row in df_subset.iterrows():
                value_to_check = row['Unnamed: 0']
                value_found = False
                for separate_index, separate_row in orthogroups.iterrows():
                    try:
                        if pd.notna(separate_row[column_to_search]) and value_to_check in separate_row[column_to_search]:
                            separate_value = separate_row['Arabidopsis']
                            values_list.append(separate_value)
                            value_found = True
                            break
                    except TypeError:
                        values_list.append('nan')
                if not value_found:
                    values_list.append('nan')
            values_lists.append(values_list)

    new_column_name = 'Arabidopsis_orthologs'

    for df, values_list in zip(dataframes.values(), values_lists):
        df[new_column_name] = values_list

    excel_writer = pd.ExcelWriter(output_filename, engine='xlsxwriter')

    for df_name, df in dataframes.items():
        df.to_excel(excel_writer, sheet_name=df_name, index=False)

    excel_writer.save()
    excel_writer.close()

# Process files for 5h
process_files("cle_results/x5h/", "cle_results/CLE_5h.xlsx")

# Process files for 24h
process_files("cle_results/x24h/", "cle_results/CLE_24h.xlsx")



