
import pandas as pd
import numpy as np 
import csv

working_directory = r'C:/Users/15156/Desktop/motif code smaller validation'
c_motif_path = r"C:\Users\15156\Desktop\Motif score\c_term_motif.csv"
n_motif_path = r"C:\Users\15156\Desktop\Motif score\n_term_motif.csv"
seq_list_path = r"C:\Users\15156\Desktop\Motif score\validation_peptide_list - Copy.csv"
c_motif_list = pd.read_csv(c_motif_path)
n_motif_list = pd.read_csv(n_motif_path)
seq_list = pd.read_csv(seq_list_path)



theo_seq_storage = []
n_motif_storage = []
c_motif_storage = []
ion_matches = []


#create df to store matches
matches = pd.DataFrame()
for ind in range(0,len(seq_list)):
        sequence = seq_list['Sequence '].iloc[[ind]]
        scan = seq_list['Peptide/Scan'].iloc[[ind]]
matches_df = pd.DataFrame(columns=['Theoretical Sequence', 'Motif', 'Scan'])

matches_df_theo_seq = []
matches_df_motif = []
matches_df_scan = []

for k in range(0,len(seq_list['Sequence '])):
    for _, row in n_motif_list.iterrows():
        n_motif = row['Sequence']
        theoretical_seq = seq_list['Sequence '][k]
        s_motif = seq_list['Peptide/Scan'][k]

        if theoretical_seq.startswith(n_motif):
            
            # Create a dictionary with the matching results
            matches_df_theo_seq.append(theoretical_seq)
            matches_df_motif.append(n_motif)
            matches_df_scan.append(s_motif)
        else:
            pass




for k in range(0,len(seq_list['Sequence '])):
    for _, row in c_motif_list.iterrows():
        c_motif = row['Sequence']
        theoretical_seq = seq_list['Sequence '][k]
        s_motif = seq_list['Peptide/Scan'][k]
        print('c motif',c_motif,theoretical_seq,s_motif)
        if theoretical_seq.endswith(c_motif):
            # Create a dictionary with the matching results
            matches_df_theo_seq.append(theoretical_seq)
            matches_df_motif.append(c_motif)
            matches_df_scan.append(s_motif)

matches_df = pd.DataFrame()
matches_df['Theoretical Sequences'] = matches_df_theo_seq
matches_df['Motif'] = matches_df_motif
matches_df['Scan'] = matches_df_scan
            
matches_df.to_csv('matches_0522.csv', index=False)

    