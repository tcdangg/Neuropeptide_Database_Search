# -*- coding: utf-8 -*-
"""
Created on Mon May 22 18:30:55 2023

@author: 15156
"""

import pandas as pd 
import numpy as np 
import csv 
from matching_motif import matches_df

#reading csv 
working_directory = r'C:/Users/15156/Desktop/Motif score'
seq_file = matches_df
base_file_path = r'C:/Users/15156/Desktop/Motif score/for_Tina/perfect_spectra/fragment_matches'

seq_w_PTM_path = r"C:\Users\15156\Desktop\Motif score\sequence_list_w_PTM.csv"
seq_w_PTM = pd.read_csv(seq_w_PTM_path)
motif_seq = seq_file['Motif'].values.tolist()


seq_file = pd.merge(seq_file,seq_w_PTM, left_on=['Theoretical Sequences','Scan'],right_on=['Name (no PTM)','Scan'])

seq_input_list = seq_file['Theoretical Sequences'].values.tolist()
scan_input_list = seq_file['Scan'].values.tolist()
modded_input_list = seq_file['Name'].values.tolist()


final_report_sequence = []
final_report_motif = []
final_report_final_score = []
final_report_scan = [] 

for a in range(0,len(seq_input_list)):
    motif_input = motif_seq[a]
    sequence_input = seq_input_list[a]
    scan = scan_input_list[a]
    
    mod_name = modded_input_list[a]
    exp_file_path = base_file_path + '\\' +mod_name + '_' + str(scan) + '_fragment_report.csv'
    exp_file = pd.read_csv(exp_file_path)
    
    #store b ions 
    b_ion = [] 
    #matchesdf
    matches = pd.DataFrame()
    #extracting column out 
    motif_list = seq_file['Motif'].values.tolist()
    # y-ions: based on motif length
    seq_file['Y_Ion'] = [list(range(1, len(motif) + 1)) for motif in motif_list]
    
    # Modify 'Y_Ion' column
    for i in range(len(seq_file['Y_Ion'])):
        y_ions = seq_file['Y_Ion'][i]
        modified_y_ions = [f'y{ion}' for ion in y_ions]  #Add prefix 'y' 
        seq_file.at[i, 'Y_Ion'] = modified_y_ions
        
    #b-ions 
    for b in seq_file['Theoretical Sequences']:
        sequence_length = len(b)
        motif_length = len(motif_input)
        end_position = []
        for j in range(sequence_length - motif_length + 1):
            end = j + motif_length
            #add b in front
            end_position.append(f'b{end}')  
    
        b_ion.append(end_position)
    seq_file['B_Ions'] = b_ion
    
    #going into exp specra file 
    #extracting out ion columns 
    exp_list = exp_file['ion'].values.tolist()
    print(exp_list)
    y_ion = seq_file['Y_Ion'].values.tolist()
    
    #making into df to use isin funct 
    y_ion_df = pd.DataFrame({'Y_Ions': y_ion})
    b_ion_df = pd.DataFrame({'B_Ions': b_ion})
    #Split values in "ions" column into sep rows
    b_ion_df = b_ion_df.explode('B_Ions').reset_index(drop=True)
    y_ion_df = y_ion_df.explode('Y_Ions').reset_index(drop=True)
    
    # Extract first three characters from list2 elementsex
    exp_list_sliced = []
    
    for aa in exp_list:
        aa = aa.replace('-H2O','')
        aa = aa.replace('-NH3','')
        exp_list_sliced.append(aa)

    for u in y_ion_df['Y_Ions']:
        for o in exp_list_sliced:
            if u == o:
                matches = matches.append({'Y_Ions': u, 'exp_y_ions': u}, ignore_index=True)
    for t in b_ion_df['B_Ions']:
        for p in exp_list_sliced:
            if t == p:
                matches = matches.append({'B_Ions': t, 'exp_b_ions': p}, ignore_index=True)  
                
                
            
    matches = matches.drop_duplicates()
    matches['Motif'] = motif_input
    matches['Sequence'] = sequence_input
    matches['Scan'] = scan
    #new column to store the reference motif score db 
    matches['Motif Score'] = matches.apply(lambda row: len(row['Motif']) / len(matches['Sequence']), axis=1)
    matches['Sequence_Length'] = matches['Sequence'].str.len()
    matches['not yet'] = matches['Motif Score'] * np.sqrt(matches['Sequence_Length'])

    y_theo = matches['Y_Ions'].count()
    y_exp = matches['exp_y_ions'].count()
    b_theo = matches['B_Ions'].count() if 'B_Ions' in matches.columns else 0
    b_exp = matches['exp_b_ions'].count() if 'exp_b_ions' in matches.columns else 0
    
    if b_theo == 0 or b_exp == 0:
        matches['Sequence Coverage'] = y_exp / y_theo
    else:
        matches['Sequence Coverage'] = (y_exp + b_exp) / (y_theo + b_theo)  

    matches['Final Score'] = matches['not yet']*(matches['Sequence Coverage'])
    final_score = max(matches['Final Score'].values.tolist())
    
    final_report_sequence.append(sequence_input)
    final_report_motif.append(motif_input)
    final_report_final_score.append(final_score)
    final_report_scan.append(scan)

final_report = pd.DataFrame()
final_report['Sequence'] = final_report_sequence
final_report['Motif'] = final_report_motif
final_report['Final Score'] = final_report_final_score
final_report['Scan'] = final_report_scan

final_report = final_report.sort_values(by='Final Score',ascending=False)
final_final_report = final_report.drop_duplicates(subset=['Sequence','Scan'])

final_report.to_csv('final_report.csv', index=False)   
final_final_report.to_csv('final_final_report.csv', index=False)   
 
