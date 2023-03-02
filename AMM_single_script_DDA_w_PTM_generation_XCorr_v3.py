# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 13:53:10 2022

@author: lawashburn
"""



import csv
import pandas as pd
import re
import os
from itertools import permutations
from Bio.SeqIO.FastaIO import SimpleFastaParser
import random
import collections
import time
import numpy as np
from scipy import signal
start = time.time()

output_folder = r"C:\Users\lawashburn\Documents\DBpep_v2\XCorr_Opt\v24_large_scale_test" #folder in which all output directories will be generated
raw_converter_path =  r"C:\Users\lawashburn\Documents\DBpep_v2\results_log\formatted_MS2\CoG2_formatted.txt" #path to formatted RawConverter output
db_path = r"C:\Users\lawashburn\Desktop\ALC50_Mass_Search_Files\duplicate_removed_crustacean_database_validated_formatted20220725.fasta" #database fasta path
sample_name = 'SG2'
precursor_error_cutoff = 20 #ppm
fragment_error_cutoff = 0.02
precursor_charges = [1,2,3,4,5,6,7,8]
fragment_charges = [1,2,3,4]
h_mass = 1.00784
bin_size = 1
number_of_steps = 7
min_seq_coverage = 25
normalization=50


amidation = True
oxidation_M_status = False
pyroglu_E_status = False
pyroglu_Q_status = False
sulfo_Y_status = False
max_modifications = 1

### generating database of mods selected ###
mods = []
mod_dict = {}
if oxidation_M_status == True:
    mods.append('M(Oxidation)')
    mod_dict['M'] = 'M(Oxidation)'
if pyroglu_E_status == True:
    mods.append('E(Pyro-glu)')
    mod_dict['E'] = 'E(Pyro-glu)'
if pyroglu_Q_status == True:
    mods.append('Q(Pyro-glu)')
    mod_dict['Q'] = 'Q(Pyro-glu)'
if sulfo_Y_status == True:
    mods.append('Y(Sulfo)')
    mod_dict['Y'] = 'Y(Sulfo)'
else:
    pass

modded_aas = []
for a in mods:
    modded_aas.append(a[0])
    
### generate output folder ###
peptide_report_output = output_folder+'\\peptide_reports'
if not os.path.exists(peptide_report_output):
    os.makedirs(peptide_report_output)

### Theoretical fragment calculator ###
proton_mass = 1.00727646688
charge = 1 #fragment charge
H = 1.0078250352
O = 15.99491463
C = 12.0000000
N = 14.003074
P = 30.973762
S = 31.9720707

aa_masses = {
    'G' : C*2  + H*3  + N   + O,
    'A' : C*3  + H*5  + N   + O,
    'S' : C*3  + H*5  + N   + O*2,
    'P' : C*5  + H*7  + N   + O,
    'V' : C*5  + H*9  + N   + O,
    'T' : C*4  + H*7  + N   + O*2,
    'C' : C*3  + H*5  + N   + O   + S,
    'L' : C*6  + H*11 + N   + O,
    'I' : C*6  + H*11 + N   + O,
    'N' : C*4  + H*6  + N*2 + O*2,
    'D' : C*4  + H*5  + N   + O*3,
    'Q' : C*5  + H*8  + N*2 + O*2,
    'K' : C*6  + H*12 + N*2 + O ,
    'E' : C*5  + H*7  + N   + O*3 ,
    'M' : C*5  + H*9  + N   + O   + S ,
    'H' : C*6  + H*7  + N*3 + O ,
    'F' : C*9  + H*9  + N   + O ,
    'R' : C*6  + H*12 + N*4 + O ,
    'Y' : C*9  + H*9  + N   + O*2 ,
    'W' : C*11 + H*10 + N*2 + O ,
    'O' : C*5  + H*12 + N*2 + O*2,
    'C(Pyro-glu)' : C*3  + H * 2 + O + S,
    'Q(Pyro-glu)' : C*5  + H*5  + N + O*2,
    'E(Pyro-glu)' : C*5  + H*4 + O*3,
    'M(Oxidation)' : C*5  + H*9  + N   + O*2   + S,
    'Y(Sulfo)' :  C*9  + H*9  + N   + O*5 + S
    }

termini = {'Standard' : H * 2 + O,
'(Amidated)' : N + H * 3}
PTMs = {'C(Pyro-glu)' : H * -3 - N,
        'Q(Pyro-glu)' : H * -3 - N,
        'E(Pyro-glu)' : H * -3 - N,
        'M(Oxidation)' : O,
        'Y(Sulfo)' : S + O * 3
        }
adducts = {
    'H2O' : H * 2 + O,
    'NH3' : N + H * 3}

def check_termini_Key(dict, key):
    if key in dict.keys():
        return dict[key]
    else:
        return termini['Standard']

def check_PTM_Key(dict, key):
    if key in dict.keys():
        return dict[key]
    else:
        return 0

def check_termini(pot_mod_peptide):
    if '(' in pot_mod_peptide:
        term_start = (pot_mod_peptide.rindex('('))
        termini_ID = pot_mod_peptide[term_start:]
        termini_mass_change = check_termini_Key(termini, termini_ID)
        return termini_mass_change
    else:
        return termini['Standard']

def check_PTM(pot_mod_peptide):
    number_of_mods = pot_mod_peptide.count('(')
    if number_of_mods > 0:
        current_peptide = []
        mass_change_collection = []
        current_peptide.append(pot_mod_peptide)
        for a in range(0,number_of_mods):
            peptide_less_mod = current_peptide[-1]
            ptm_start = (peptide_less_mod.index('('))-1
            ptm_end = (peptide_less_mod.index(')'))+1
            ptm_ID = peptide_less_mod[ptm_start:ptm_end]
            ptm_mass_change = check_PTM_Key(PTMs, ptm_ID)
            mass_change_collection.append(ptm_mass_change)
            peptide_less_mod2 = peptide_less_mod[:ptm_start] + peptide_less_mod[ptm_end:]
            current_peptide.append(peptide_less_mod2)
            
        ptm_mass_change = sum(mass_change_collection)
        return ptm_mass_change
    else:
        ptm_mass_change = 0
        return ptm_mass_change

def list_of_residues(pot_mod_peptide):
    list_of_res = []
    pep_update = []
    pep_update.append(pot_mod_peptide)
    no_mods = pot_mod_peptide.count('(')
    if no_mods > 1:
        for c in range(0,no_mods+1):
            pep_of_interest = pep_update[-1]
            if '(' in pep_of_interest:
                first_ptm_start = pep_of_interest.index('(')
                first_ptm_end = pep_of_interest.index(')')

                first_residues = pep_of_interest[:(first_ptm_start-1)]
                for a in first_residues:
                    list_of_res.append(a)
                ptm_residue = pep_of_interest[(first_ptm_start-1):(first_ptm_end+1)]
                list_of_res.append(ptm_residue)
                remaining_pep = pep_of_interest[(first_ptm_end+1):]
                pep_update.append(remaining_pep)  
            else:
                for d in pep_of_interest:
                    list_of_res.append(d)
    elif no_mods == 1:
        for c in range(0,1):
            pep_of_interest = pep_update[-1]
            if '(' in pep_of_interest:
                first_ptm_start = pep_of_interest.index('(')
                first_ptm_end = pep_of_interest.index(')')
                if first_ptm_start == 1:
                    ptm_residue =  pep_of_interest[0] + (pep_of_interest[(first_ptm_start):(first_ptm_end+1)])
                    list_of_res.append(ptm_residue)
                    remaining_pep = pep_of_interest[(first_ptm_end+1):]
                    for d in remaining_pep:
                        list_of_res.append(d)
                if first_ptm_start != 1:
                    first_residues = pep_of_interest[:(first_ptm_start-1)]
                    for a in first_residues:
                        list_of_res.append(a)
                    ptm_residue = pep_of_interest[(first_ptm_start-1):(first_ptm_end+1)]
                    list_of_res.append(ptm_residue)
                    remaining_pep = pep_of_interest[(first_ptm_end+1):]
                    for d in remaining_pep:
                        list_of_res.append(d)              
            else:
                for d in pep_of_interest:
                    list_of_res.append(d) 
    elif no_mods == 0:
        for c in pot_mod_peptide:
            list_of_res.append(c)
    return list_of_res
### end of theoretical fragment calculator ###
### start of monoisotopic mass calculator ###
def monoisotopic_mass_calculator(peptide_from_fasta):
        plain_peptide = re.sub("[\(\[].*?[\)\]]", "",peptide_from_fasta)
        res_list_for_fragment = list_of_residues(peptide_from_fasta)
        mass_of_residues = []
        for residue in plain_peptide:
            residue_mass = aa_masses[residue]
            mass_of_residues.append(residue_mass)
        peptide_mass = (sum(mass_of_residues)) + check_termini(peptide_from_fasta) + check_PTM(peptide_from_fasta)
        mass_to_charge = (peptide_mass + (proton_mass * charge))/charge
        return mass_to_charge
### end of monoisotopic mass calculator ###
### start of convert fasta to dataframe with monoisotopic masses ###
fasta_to_df = []
with open(db_path) as fasta_file:  # Will close handle cleanly
    for title, sequence in SimpleFastaParser(fasta_file):
        fasta_to_df.append(sequence)
um_database_list=fasta_to_df
###start of decoy db generation function###
#um_database_list = db['Sequence'].to_list()
modified_db_list = []

for db in um_database_list:
    while '(' in db:
        db = db.replace(db[db.index('('):db.index(')')+1],'')
    modified_db_list.append(db)
   
database = modified_db_list 
###


class Identity_Threshold:
    def ID_threshold_check(decoy_seqeunce):
        """
        ### Local alignment with the decoy sequence and all possible target sequences in the database.
        
        Args:
            decoy_sequence (str): decoy sequence generated from the target-decoy method
        Return:
            True: if the decoy sequence has less than 50% similarity with all target sequences in the database
            False: if the decoy sequence has more than 50% similarity with any of target sequence in the database
        """
        filtered_target = []
        if len(decoy_seqeunce) == 2:    ### bypassing sequence length with 2
            return True
        for t in database:  ### loops through all database sequence that has the same length as decoy sequence
            if len(decoy_seqeunce) == len(t):
                filtered_target.append(t)
        for target_sequence in filtered_target:
            alignment_score = 0
            for amino_acid in list(zip(target_sequence,decoy_seqeunce)):
                if amino_acid[0] == amino_acid[1]:
                    alignment_score+=1
            if alignment_score/len(decoy_seqeunce) < 0.5:
                continue
            else:
                return False
        return True
    
    def max_ID_threshold(decoy_sequence):
        """
        ### Find the maximum identity threshold percentage from the decoy sequence taht matches to target sequences
        
        Args:
            decoy_sequence (str): decoy sequence generated from the target-decoy method
        Return:
            list:   max_ID_T (float): maximum identity threshold percentage
                    target_index (int): index of corresponding target sequence that yields the maximum identity threshold percentage
        """
        max_ID_T = 0
        target_list_all = database
        target_list = []
        for t in target_list_all:
            if len(t) == len(decoy_sequence):
                target_list.append(t)
        for target_index in range(len(target_list)):
            target_sequence = target_list[target_index]
            alignment_score = 0
            for amino_acid in list(zip(target_sequence,decoy_sequence)):
                if amino_acid[0] == amino_acid[1]:
                    alignment_score+=1
            alignment_ID_threshold = alignment_score/len(decoy_sequence)
            if alignment_ID_threshold > max_ID_T:
                max_ID_T = alignment_ID_threshold
        return [max_ID_T, target_index]
    
    def alignment_sequence_index(decoy_sequence):
        """
        ### Gives the matched residues index from target and decoy sequences
        
        Args:
            decoy_sequence (str): decoy sequence generated from the target-decoy method
        Return:
            match_index_list (list): list contains the index of matched residues from both target and decoy sequences
        """
        max_num, index = Identity_Threshold.max_ID_threshold(decoy_sequence)
        target_sequence = database[index]
        x = 0
        match_index_list = []
        for amino_acid in list(zip(target_sequence,decoy_sequence)):
            if amino_acid[0] == amino_acid[1]:
                match_index_list.append(x)
            x+=1
        return match_index_list 
    
    
def shuffled_decoy_database():
    """
    ### Residues in each target sequences are randomly shuffled and generated as the decoy sequences
    
    Return:
        shuffled_decoy_list (list): Shuffle decoy database with Identify Threshold applied
    """
    shuffled_decoy_list = []
    
    def shuffled_algorithm(database_sequence):
        """
        ### Shuffles each amino acid residue in each target sequences
        
        Args:
            database_seqeuence (str): target sequence from the database
        
        Return:
            shuffled_database_sequence (str): decoy seqeunce from shuffle decoy method
        """
        shuffled_database_sequence = ''
        for i in random.sample(range(0,len(database_sequence)),len(database_sequence)):
            shuffled_database_sequence += database_sequence[i]
        return shuffled_database_sequence

    shuffled_table = []
    # original shuffled decoy list
    for database_sequence in database:
        shuffled_decoy_sequence = shuffled_algorithm(database_sequence)
        if len(database_sequence) == 2:         ### to make sure sequence with 2 residues have different decoy sequence with original database sequence
            while database_sequence==shuffled_decoy_sequence:
                shuffled_decoy_sequence =shuffled_algorithm(database_sequence) 
                
        shuffled_decoy_list.append(shuffled_decoy_sequence)
        shuffled_table.append([shuffled_decoy_sequence, round(Identity_Threshold.max_ID_threshold(shuffled_decoy_sequence)[0],3)])
        
    iterations = 1
    over_threshold_list = []
    for decoy_index in range(len(shuffled_decoy_list)):
        decoy_sequence = shuffled_decoy_list[decoy_index]
        if len(decoy_sequence) == 2:    #bypassing sequence with length of 2
            continue
        if not Identity_Threshold.ID_threshold_check(decoy_sequence):
            over_threshold_list.append([decoy_index, decoy_sequence])
    
    while iterations < 31:
        if len(over_threshold_list) == 0:
            break
        reshuffled_list = []
        for index,sequence in over_threshold_list:   
            reshuffled = shuffled_algorithm(sequence)
            reshuffled_list.append([index, sequence, reshuffled])
        for index,sequence,reshuffled in reshuffled_list:
            if Identity_Threshold.ID_threshold_check(reshuffled):
                shuffled_decoy_list[index] = reshuffled
                over_threshold_list.remove([index,sequence])
        if len(over_threshold_list) == 0:
            iterations = 30
        iterations+=1
        
    if len(over_threshold_list) > 0:
        all_amino_acid_left = []
        for i, seq in over_threshold_list:
            for s in seq:
                all_amino_acid_left.append(s)
    
        if len(over_threshold_list) > 0:
            for index, sequence in over_threshold_list:
                while not Identity_Threshold.ID_threshold_check(sequence):
                    index_list = Identity_Threshold.alignment_sequence_index(sequence)
                    mutation = sequence.replace(sequence[index_list[0]],random.choice(all_amino_acid_left))        ### double check this
                    if Identity_Threshold.ID_threshold_check(mutation):
                        sequence = mutation
                        shuffled_decoy_list[index] = mutation
                    else:
                        continue
    return shuffled_decoy_list

def target_decoy_dataframe():
    target_decoy_df = pd.DataFrame()
    target_decoy_df['Target DB Sequence'] = database
    target_decoy_df['Decoy DB Sequence'] = shuffled_decoy_database()
    
    decoy_entries = target_decoy_df['Decoy DB Sequence'].values.tolist()
    
    
    target_status = []
    for xx in fasta_to_df:
        target_status.append('Target')
    decoy_entries = target_decoy_df['Decoy DB Sequence'].values.tolist()
    decoy_status = []
    for jj in decoy_entries:
        decoy_status.append('Decoy')
    fasta_to_df.extend(decoy_entries)
    target_status.extend(decoy_status)
    
    database_export_report = pd.DataFrame()
    database_export_report['Sequence'] = fasta_to_df
    database_export_report['Status'] = target_status
   
    file_out_path = output_folder + '\\target_decoy_DB.csv'
    with open(file_out_path,'w',newline='') as filec:
            writerc = csv.writer(filec)
            database_export_report.to_csv(filec,index=False)
    
    return target_decoy_df
def param_log_export():

    out_path_report = 'Ouput directory path: ' + output_folder
    rawconverter_path_report = 'RawConverter output file path: ' + raw_converter_path    
    db_path_report = 'Database path: ' + db_path
    sample_ID_report = 'Sample name: ' + sample_name
    p_err_report = 'Precursor error threshold (ppm): ' + str(precursor_error_cutoff)
    f_err_report = 'Fragment error threshold (Da): ' + str(fragment_error_cutoff)
    max_p_charge_report = 'Maximum precursor charge: +' + str(precursor_charges[-1])
    max_f_charge_report = 'Maximum fragment charge: +' + str(fragment_charges[-1])
    bin_size_report = 'Bin size: ' + str(bin_size)
    no_steps_report = 'Number of steps: ' + str(number_of_steps)
    normalization_report = 'Normalized intensity : ' + str(normalization)    
    max_modifications_report = 'Maximum number of modifications : ' + str(max_modifications)    
    amidation_modifications_report = 'Amidation? : ' + str(amidation) 
    oxidation_M_status_modifications_report = 'Oxidation on M? : ' + str(oxidation_M_status) 
    pyroglu_E_status_modifications_report = 'Pyro-glu on E? : ' + str(pyroglu_E_status) 
    pyroglu_Q_status_modifications_report = 'Pyro-glu on Q? : ' + str(pyroglu_Q_status) 
    sulfo_Y_status_modifications_report = 'Sulfo on Y? : ' + str(sulfo_Y_status) 
    elapsed_time_report = 'Time elapsed : ' + str(elapsed_time) 
    
    param_file_entries = [out_path_report,rawconverter_path_report,db_path_report,sample_ID_report,p_err_report,f_err_report,max_p_charge_report,max_f_charge_report,
                          bin_size_report,no_steps_report,normalization_report,max_modifications_report,amidation_modifications_report,oxidation_M_status_modifications_report,
                          pyroglu_E_status_modifications_report,pyroglu_Q_status_modifications_report,sulfo_Y_status_modifications_report,elapsed_time_report]

    param_file_path = output_folder + '\\' + 'parameter_file.txt'
    with open(param_file_path,'a') as f:
        f.writelines('\n'.join(param_file_entries))
###end of decoy db generation function###
target_decoy_dataframe()

final_seq_list = []
for b in fasta_to_df:
    aas_present = []
    aas_present_index = []
    aas_absent = []
    aas_absent_index = []
    
    seq_len = len(b)
    for c in range(0,seq_len):
        aa_ID = b[c]
        if aa_ID in modded_aas:
            aas_present.append(aa_ID)
            aas_present_index.append(str(c))
            aas_absent.append(aa_ID)
            aas_absent_index.append(str(c))
        else:
            aas_absent.append(aa_ID)
            aas_absent_index.append(str(c))

    number_mods = len(aas_present)
    
    seq_log = pd.DataFrame()
    seq_log['Original Residue'] = aas_present
    seq_log['Index'] = aas_present_index
    seq_log['Res+Index'] = seq_log[['Original Residue', 'Index']].apply(lambda x: ''.join(x), axis=1)
    seq_w_ind = seq_log['Res+Index'].values.tolist()
    
    seq_absent_log = pd.DataFrame()
    seq_absent_log['Original Residue'] = aas_absent
    seq_absent_log['Index'] = aas_absent_index

    for d in range(0,number_mods):
        seq_w_ind.append('X'+(str(d)))
    
    aa_combination_no_dups = []
    
    if number_mods > max_modifications:
        aa_combinations = permutations(seq_w_ind,max_modifications)
        aa_combo_list = list(aa_combinations)
        for a in aa_combo_list:
                if a not in aa_combination_no_dups:
                    aa_combination_no_dups.append(a)
    else:
        aa_combinations = permutations(seq_w_ind,number_mods)
        aa_combo_list = list(aa_combinations)
        for a in aa_combo_list:
                if a not in aa_combination_no_dups:
                    aa_combination_no_dups.append(a)
            
    restyle_aa_list = []
    restyle_aa_ind_list = []
    for z in aa_combination_no_dups:
        restyle_aa = []
        restyle_aa_ind = []
        for y in z:
            restyle_aa.append(y[0])
            restyle_aa_ind.append(y[1:])
        restyle_aa_list.append(restyle_aa)
        restyle_aa_ind_list.append(restyle_aa_ind)
    
    unmodded_seq = pd.DataFrame()
    unmodded_seq['Original Residue'] = aas_absent
    unmodded_seq['Index'] = aas_absent_index
    
    all_modded_seqs = []
    for k in range(0,len(restyle_aa_list)):
        new_seq_log = pd.DataFrame()
        new_seq_log['New Residue'] = restyle_aa_list[k]
        new_seq_log['Index'] = restyle_aa_ind_list[k]
        new_seq_log_filtered = new_seq_log[new_seq_log['New Residue'] != 'X']
        ptm_applied_seq_log = new_seq_log_filtered.replace({'New Residue':mod_dict})
        
        merge_seq_log = unmodded_seq.merge(ptm_applied_seq_log, on='Index', how='left')
        ptm_applied_seq_log = merge_seq_log.sort_values(by='Index')
        ptm_applied_seq_log['New Residue'] = ptm_applied_seq_log['New Residue'].replace('', pd.NA).fillna(ptm_applied_seq_log['Original Residue'])
    
        modified_seq = ptm_applied_seq_log['New Residue'].values.tolist()
        modified_seq_format = "".join([str(item) for item in modified_seq])
        all_modded_seqs.append(modified_seq_format)

    all_modded_seqs_nodups = []
    for l in all_modded_seqs:
        if l not in all_modded_seqs_nodups:
            all_modded_seqs_nodups.append(l)
    
    if amidation == True:
        for m in all_modded_seqs_nodups:
            final_seq_list.append(m)
            amidated_pep = m + '(Amidated)'
            final_seq_list.append(amidated_pep)
    else:
        for m in all_modded_seqs_nodups:
            final_seq_list.append(m)
            

###end of PTM generation###

###end of PTM generation###
###start of mass calculations ###
fasta_monoiso_mass = []

for sequence in final_seq_list:
    mono_mass = monoisotopic_mass_calculator(sequence)
    fasta_monoiso_mass.append(mono_mass)

db = pd.DataFrame()
db['Sequence'] = final_seq_list
db['Monoisotopic Mass'] = fasta_monoiso_mass

### end of convert fasta to dataframe with monoisotopic masses ###
### start of importing spectra output from RawConverter ###
raw_converter = pd.read_csv(raw_converter_path, sep=",",skiprows=[0], names= ["m/z","resolution","charge","intensity", "MS2", "scan_number","precursor_charge","null"])
raw_converter = raw_converter[raw_converter['charge'] != 0]

precursor_mz = []
precursor_z = []
precursor_scan = []

rawconv_mz = raw_converter['MS2'].values.tolist()
rawconv_z = raw_converter['precursor_charge'].values.tolist()
rawconv_scan = raw_converter['scan_number'].values.tolist()

for b in rawconv_mz:
    precursor_mz.append(b)

for d in rawconv_z:
    precursor_z.append(d)

for e in rawconv_scan:
    precursor_scan.append(e)

exp_precursor = pd.DataFrame()
exp_precursor['Precursor actual m/z'] = precursor_mz
exp_precursor['Precursor actual charge'] = precursor_z
exp_precursor['Precursor scan'] = precursor_scan


exp_precursor = exp_precursor.drop_duplicates() 
exp_precursor['Monoisotopic Mass'] =  ((exp_precursor['Precursor actual m/z']) * (exp_precursor['Precursor actual charge']))-(h_mass*(exp_precursor['Precursor actual charge']))

### end of importing spectra output from RawConverter ###

if len(db)<1: #throws an error if database file is empty
    raise ValueError('Database file is empty')

precursor_temp_cutoff = precursor_error_cutoff/100 #rough estimate of ppm to Da to minimize extra search space

### start of precursor AMM ###

precursor_amm_actual_mz = []
precursor_amm_actual_z = []
precursor_amm_actual_scan = []
precursor_amm_actual_monoiso = []
precursor_amm_theoretical_sequence = []
precursor_amm_theoretical_monoiso = []
precursor_amm_err = []

for a in precursor_charges:
    db_sorted = db.sort_values(by = 'Monoisotopic Mass')
    db_sorted = db_sorted.rename(columns={'Monoisotopic Mass':'Theoretical Monoisotopic Mass'})

    exp_precursor_sorted = exp_precursor.sort_values(by = 'Monoisotopic Mass')
    exp_precursor_sorted = exp_precursor_sorted.rename(columns={'Monoisotopic Mass':'Actual Monoisotopic Mass'})

    exp_precursor_z_filter = exp_precursor_sorted[exp_precursor_sorted['Precursor actual charge'] == a]

    merge_match = pd.merge_asof(exp_precursor_z_filter,db_sorted, left_on='Actual Monoisotopic Mass', right_on='Theoretical Monoisotopic Mass',
                                tolerance = precursor_temp_cutoff, allow_exact_matches=True)

    merge_match_filtered = merge_match.dropna(subset=['Sequence','Theoretical Monoisotopic Mass'])
    merge_match_filtered['Precursor error (ppm)'] = ((abs((merge_match_filtered['Theoretical Monoisotopic Mass'])-(merge_match_filtered['Actual Monoisotopic Mass'])))/
                                                      (merge_match_filtered['Theoretical Monoisotopic Mass'])) * 1E6
    merge_match_filtered2 = merge_match_filtered[merge_match_filtered['Precursor error (ppm)'] <= precursor_error_cutoff]

    actual_mz = merge_match_filtered2['Precursor actual m/z'].values.tolist()
    actual_z = merge_match_filtered2['Precursor actual charge'].values.tolist()
    actual_scan = merge_match_filtered2['Precursor scan'].values.tolist()
    actual_monoiso = merge_match_filtered2['Actual Monoisotopic Mass'].values.tolist()
    theoretical_sequence = merge_match_filtered2['Sequence'].values.tolist()
    theoretical_monoiso = merge_match_filtered2['Theoretical Monoisotopic Mass'].values.tolist()
    err = merge_match_filtered2['Precursor error (ppm)'].values.tolist()
    
    for f in actual_mz:
        precursor_amm_actual_mz.append(f)
    for g in actual_z:
        precursor_amm_actual_z.append(g)
    for h in actual_scan:
        precursor_amm_actual_scan.append(h)
    for i in actual_monoiso:
        precursor_amm_actual_monoiso.append(i)
    for j in theoretical_sequence:
        precursor_amm_theoretical_sequence.append(j)
    for k in theoretical_monoiso:
        precursor_amm_theoretical_monoiso.append(k)
    for l in err:
        precursor_amm_err.append(l)

precursor_amm_results = pd.DataFrame()
precursor_amm_results['Precursor Actual m/z'] = precursor_amm_actual_mz
precursor_amm_results['Precursor Actual z'] = precursor_amm_actual_z
precursor_amm_results['Scan'] = precursor_amm_actual_scan
precursor_amm_results['Precursor Actual Monoisotopic'] = precursor_amm_actual_monoiso
precursor_amm_results['Sequence'] = precursor_amm_theoretical_sequence
precursor_amm_results['Precursor Theoretical Monoisotopic'] = precursor_amm_theoretical_monoiso
precursor_amm_results['Precursor error (ppm)'] = precursor_amm_err

secondary_amm_prep = pd.merge(raw_converter,precursor_amm_results,left_on=['MS2','scan_number','precursor_charge'], right_on=['Precursor Actual m/z','Scan','Precursor Actual z'])
secondary_amm_prep_clean = secondary_amm_prep.drop(['MS2','scan_number','precursor_charge','null','resolution'],axis=1)
secondary_amm_prep_clean = secondary_amm_prep_clean.rename(columns={'charge':'Fragment actual charge','m/z':'Fragment actual m/z'})

candidate_sequences_raw = secondary_amm_prep_clean['Sequence'].values.tolist()

candidate_sequences = []

final_report_storage_scan = []
final_report_storage_seq_coverage = []
final_report_storage_peptide = []

for m in candidate_sequences_raw:
    if m not in candidate_sequences:
        candidate_sequences.append(m)
### end of precursor AMM ###

### start of fragment AMM ###
corr_correlation = []
corr_pep = []
corr_scan = []
correlation_max_store = []
correlation_mean_store = []
correlation_array_store = []
log10_corr_store = []
for peptide in candidate_sequences:
        ###pull theoretical masses###
        plain_peptide = re.sub("[\(\[].*?[\)\]]", "",peptide) #removes any modification for mass calculations
        res_list_for_fragment = list_of_residues(peptide)
        mass_of_residues = []
        for residue in plain_peptide:
            residue_mass = aa_masses[residue]
            mass_of_residues.append(residue_mass)
        peptide_mass = (sum(mass_of_residues)) + check_termini(peptide) + check_PTM(peptide) #calculated MH mass
        mass_to_charge = (peptide_mass + (proton_mass * charge))/charge #translates the MH mass to m/z for each charge of interest

        num_ions = len(plain_peptide)-1 #number of expected fragment ions is one less than the number of AAs in the sequence

        b_ions = []
        b_ion_name = []
        
        y_ions = []
        y_ion_name = []
        
        for a in range(0,num_ions):
            
            residue_identity = res_list_for_fragment[a]
            if len(b_ions) == 0:
                ion_mass = aa_masses[residue_identity]
                ion_mz = ion_mass + proton_mass
                b_ions.append(ion_mz)
                ion_name = 'b' + str(a+1)
                b_ion_name.append(ion_name)
            
            elif len(b_ions) > 0:
                ion_mass = (aa_masses[residue_identity]) + b_ions[-1]
                b_ions.append(ion_mass)
                ion_name = 'b' + str(a+1)
                b_ion_name.append(ion_name)
        
        for b in (range(0,num_ions)):
            residue_identity = res_list_for_fragment[b]
            if len(y_ions) == 0:
                ion_mass = mass_to_charge - aa_masses[residue_identity]
                y_ions.append(ion_mass)
                ion_name = 'y' + str((num_ions-b))
                y_ion_name.append(ion_name)
            elif len(y_ions) > 0:
                ion_mass = y_ions[-1] - aa_masses[residue_identity]
                y_ions.append(ion_mass)
                ion_name = 'y' + str((num_ions-b))
                y_ion_name.append(ion_name)

        b_ions_report = pd.DataFrame()
        b_ions_report['ion'] = b_ion_name
        b_ions_report['mass'] = b_ions
        
        b_ions_water_adduct = pd.DataFrame()
        b_ions_water_adduct['ion'] = b_ions_report['ion'] + '-H2O'
        b_ions_water_adduct['mass'] = b_ions_report['mass'] - adducts['H2O']
        
        b_ions_ammonia_adduct = pd.DataFrame()
        b_ions_ammonia_adduct['ion'] = b_ions_report['ion'] + '-NH3'
        b_ions_ammonia_adduct['mass'] = b_ions_report['mass'] - adducts['NH3']
        
        y_ions_report = pd.DataFrame()
        y_ions_report['ion'] = y_ion_name
        y_ions_report['mass'] = y_ions
        
        y_ions_ammonia_adduct = pd.DataFrame()
        y_ions_ammonia_adduct['ion'] = y_ions_report['ion'] + '-NH3'
        y_ions_ammonia_adduct['mass'] = y_ions_report['mass'] - adducts['NH3']

        y_ions_water_adduct = pd.DataFrame()
        y_ions_water_adduct['ion'] = y_ions_report['ion'] + '-H2O'
        y_ions_water_adduct['mass'] = y_ions_report['mass'] - adducts['H2O']
        
        ion_report = pd.DataFrame()
        ion_report = ion_report.append(b_ions_report)
        ion_report = ion_report.append(y_ions_report)
        ion_report = ion_report.append(b_ions_water_adduct)
        ion_report = ion_report.append(b_ions_ammonia_adduct)
        ion_report = ion_report.append(y_ions_ammonia_adduct)
        ion_report = ion_report.append(y_ions_water_adduct)
        ion_report = ion_report.drop_duplicates()
        ion_report = ion_report.rename(columns={'mass':'Fragment theoretical Monoisotopic Mass'})
        ### end of pulling theoretical fragments ###

        filtered_secondary_amm = secondary_amm_prep_clean[secondary_amm_prep_clean['Sequence'] == peptide] #filter amm from before for the sequence we are looking at

        filtered_secondary_amm['Fragment Actual Monoisotopic Mass'] = (filtered_secondary_amm['Fragment actual m/z'] * filtered_secondary_amm['Fragment actual charge']) - (h_mass*filtered_secondary_amm['Fragment actual charge'])     

        scans_present_raw = filtered_secondary_amm['Scan'].values.tolist()
        
        scans_present = []
        for z in scans_present_raw:
            if z not in scans_present:
                scans_present.append(z)
        
        for y in scans_present:
            scan_to_report = y
            final_report_storage_scan.append(y)
            scans_filtered_secondary_amm = filtered_secondary_amm[filtered_secondary_amm['Scan'] == y] #filter to look at just one scan
            scans_filtered_secondary_amm = scans_filtered_secondary_amm.sort_values(by='Fragment Actual Monoisotopic Mass')

            ion_report = ion_report.sort_values(by='Fragment theoretical Monoisotopic Mass')
            prelim_fragment_matches = pd.merge_asof(scans_filtered_secondary_amm,ion_report,left_on='Fragment Actual Monoisotopic Mass',
                                                    right_on='Fragment theoretical Monoisotopic Mass', tolerance=fragment_error_cutoff,allow_exact_matches=True)
            merge_fragment_match_filtered = prelim_fragment_matches.dropna(subset=['ion','Fragment theoretical Monoisotopic Mass'])
            #merge_fragment_match_filtered = prelim_fragment_matches
            
            merge_fragment_match_filtered['Fragment error (Da)'] = merge_fragment_match_filtered['Fragment Actual Monoisotopic Mass'] - merge_fragment_match_filtered['Fragment theoretical Monoisotopic Mass']
            merge_fragment_match_filtered = merge_fragment_match_filtered[merge_fragment_match_filtered['Fragment error (Da)'] <= fragment_error_cutoff]
            
            if len(merge_fragment_match_filtered)>0:
                file_path = peptide_report_output+'\\'+peptide+'_'+str(y)+'_prelim_fragment_matches.csv'
                with open(file_path,'w',newline='') as filec:
                        writerc = csv.writer(filec)
                        merge_fragment_match_filtered.to_csv(filec,index=False) 

            merge_fragment_match_filtered['ion_format'] = merge_fragment_match_filtered['ion'].str.extract('(\d+)', expand=False)
            ion_report['ion_format'] = ion_report['ion'].str.extract('(\d+)', expand=False)
            
            merge_fragment_match_filtered_ion = merge_fragment_match_filtered.drop_duplicates(subset=['ion_format'],keep='first').reset_index(drop=True)
            ion_report_filtered_ion = ion_report.drop_duplicates(subset=['ion_format'],keep='first').reset_index(drop=True)
            
            seq_coverage = (len(merge_fragment_match_filtered_ion)/len(ion_report_filtered_ion))*100
            final_report_storage_seq_coverage.append(seq_coverage)
            final_report_storage_peptide.append(peptide)
            
            scan_filtered_score = merge_fragment_match_filtered[merge_fragment_match_filtered['Scan'] == y]
            scan_filtered_pep = scan_filtered_score['Sequence'].values.tolist()
            actual_peak_intensities = raw_converter[raw_converter['scan_number'] == y]
            
            for l in scan_filtered_pep:
                ###pull theoretical masses again###
                ###start of XCorr###
                ion_report = ion_report.rename(columns={'mass':'Fragment theoretical Monoisotopic Mass'}) #generate theoretical spectra with intensity of 50
                ion_report['Intensity'] = 50

                full_theoretical_spectra = pd.DataFrame()
                
                for e in fragment_charges:
                    ion_report['Fragment z'] = e
                    ion_report['Fragment theoretical m/z'] = (ion_report['Fragment theoretical Monoisotopic Mass'] + (h_mass*e))/e
                    full_theoretical_spectra = full_theoretical_spectra.append(ion_report)


                max_min = []
                max_min.append(full_theoretical_spectra['Fragment theoretical m/z'].max())
                max_min.append(actual_peak_intensities['m/z'].max())
                max_min.append(full_theoretical_spectra['Fragment theoretical m/z'].min())
                max_min.append(actual_peak_intensities['m/z'].min())
                
                true_min = min(max_min)
                true_max = max(max_min)
                
                round_min = (int(true_min))-1
                round_max = (int(true_max))+1
                mz_range = true_max - true_min
                
                step_sizes = []

                for x in range(0,number_of_steps):
                    if len(step_sizes) == 0:
                        first_step = mz_range/number_of_steps
                        step_sizes.append(float(first_step))
                    else:
                        step_sizes.append(float(step_sizes[-1]/number_of_steps))

                confirmed_step_sizes = []
                for y in step_sizes:
                    if y > bin_size:
                        confirmed_step_sizes.append(float(y))
                confirmed_step_sizes.append(bin_size)
                confirmed_step_sizes = sorted(confirmed_step_sizes,reverse=True, key=float)

                initial_index = []
                stepped_index = pd.DataFrame()

                for f in range(round_min,round_max):
                    initial_index.append(f)

                for w in confirmed_step_sizes:
                    w_index = confirmed_step_sizes.index(w)
                    if w_index == 0:
                        numerator = mz_range
                    else:
                        numerator = confirmed_step_sizes[w_index-1]
                    
                    bin_index = []
                    number_bins = numerator/w
                    for o in initial_index:
                        for n in range(0,int(number_bins)+2):
                            index = o + (n * w)
                            index = round(index,2)
                            bin_index.append(index)

                    
                    index_assess = pd.DataFrame()
                    index_assess['Bin #'] = bin_index
                    
                    full_theoretical_spectra = full_theoretical_spectra.sort_values(by='Fragment theoretical m/z')
                    actual_peak_intensities = actual_peak_intensities.sort_values(by='m/z')
                    index_assess = index_assess.sort_values(by='Bin #')
                    
                    theoretical_index = pd.merge_asof(index_assess,full_theoretical_spectra, left_on='Bin #', right_on='Fragment theoretical m/z',
                                                tolerance = w, allow_exact_matches=True)
                    theoretical_index = pd.merge_asof(theoretical_index,actual_peak_intensities, left_on='Bin #', right_on='m/z',
                                                tolerance = w, allow_exact_matches=True)
                    theoretical_index = theoretical_index.drop(columns=['resolution', 'MS2','precursor_charge','null'])
                    
                    theoretical_index = theoretical_index.rename(columns={"ion": "Fragment ion", "Intensity": "Theoretical intensity","Fragment z":"Theoretical fragment z",
                                                                          "m/z":'Actual m/z',"charge":"Actual fragment charge","intensity":"Fragment actual intensity",
                                                                          "scan_number":"Scan #"})
                    
                    theoretical_index = theoretical_index.dropna(subset=['Fragment theoretical m/z', 'Actual m/z'], how='all')
                    theoretical_index['Step Size'] = w
                    theoretical_index['Bin #'] = theoretical_index['Bin #'].round(2)
                    theoretical_index = theoretical_index.drop_duplicates()
                    stepped_index = theoretical_index
                    post_step_index = theoretical_index['Bin #'].values.tolist()
                    initial_index.clear()
                    initial_index.extend(post_step_index)
                
                stepped_index['Theoretical intensity'] = stepped_index['Theoretical intensity'].replace(np.nan, 0)
                stepped_index['Fragment actual intensity'] = stepped_index['Fragment actual intensity'].replace(np.nan, 0)

                # theo_favor_step_index = stepped_index[stepped_index['Fragment actual intensity'] == 0]
                # theo_favor_step_index = theo_favor_step_index[theo_favor_step_index['Theoretical intensity'] != 0]                
                # theo_favor_step_index['Fragment ion'] = theo_favor_step_index['Fragment ion'].str.replace('-H2O','')
                # theo_favor_step_index['Fragment ion'] = theo_favor_step_index['Fragment ion'].str.replace('-NH3','')
                # theo_favor_step_index = theo_favor_step_index.drop_duplicates(subset='Fragment ion')

                filter_xcorr_input = stepped_index
                filter_xcorr_input = filter_xcorr_input.sort_values(by='Fragment actual intensity',ascending=False)
                filter_xcorr_input = filter_xcorr_input.drop_duplicates()
                normalization_factor = filter_xcorr_input['Fragment actual intensity'].max()
                normalization_multiple_factor = normalization/normalization_factor
                #filter_xcorr_input['Fragment normalized intensity'] = 50
                
                filter_xcorr_input['Fragment normalized intensity'] = filter_xcorr_input['Fragment actual intensity']*normalization_multiple_factor #can be modified and removed on a later version to save time, just updated so that if actual intensity >0, normalized intensity =50
                #filter_xcorr_input['Fragment normalized intensity'].values[filter_xcorr_input['Fragment normalized intensity']> 0] = 50 #within each bin the highest peak is equal to 50 for intensity
                filter_xcorr_input['Bin #'] = filter_xcorr_input['Bin #'].round(2)
                filter_xcorr_input = filter_xcorr_input.drop_duplicates()
                filter_xcorr_input = filter_xcorr_input.drop_duplicates(subset=['Fragment ion', 'Fragment theoretical Monoisotopic Mass','ion_format',
                                                                                'Theoretical intensity','Theoretical fragment z','Fragment theoretical m/z',
                                                                                'Actual m/z','Actual fragment charge','Fragment actual intensity',
                                                                                'Scan #','Step Size','Fragment normalized intensity'])
                #filter_xcorr_input = filter_xcorr_input.sort_values(by = ['Fragment theoretical m/z', 'Fragment actual intensity'], ascending = [True, False])
                #filter_xcorr_input = filter_xcorr_input[(~filter_xcorr_input['Fragment theoretical m/z'].duplicated()) | filter_xcorr_input['Fragment theoretical m/z'].isna()] #drop Fragment theoretical duplicates without dropping nans
                
                #filter_xcorr_input = filter_xcorr_input.sort_values(by = ['Actual m/z', 'Theoretical intensity'], ascending = [True, False])
                #filter_xcorr_input = filter_xcorr_input[(~filter_xcorr_input['Actual m/z'].duplicated()) | filter_xcorr_input['Actual m/z'].isna()] #drop Fragment theoretical duplicates without dropping nans
                
                #filter_xcorr_input['Fragment ion no ion'] =  filter_xcorr_input['Fragment ion']
                #filter_xcorr_input['Fragment ion no ion'] = filter_xcorr_input['Fragment ion no ion'].str.replace('-H2O','')
                #filter_xcorr_input['Fragment ion no ion'] = filter_xcorr_input['Fragment ion no ion'].str.replace('-NH3','')
                
                #filter_xcorr_input = filter_xcorr_input.sort_values(by = ['Fragment ion no ion', 'Fragment normalized intensity'], ascending = [True, False])
                #filter_xcorr_input = filter_xcorr_input[(~filter_xcorr_input['Fragment ion no ion'].duplicated()) | filter_xcorr_input['Fragment ion no ion'].isna()] #drop Fragment theoretical duplicates without dropping nans
                
                
                #max_verified_intensity_get = filter_xcorr_input[filter_xcorr_input['Theoretical intensity'] != 0] #determines the maximum actual intensity associated with an expected fragment ion
                #max_verified_intensity = max_verified_intensity_get['Fragment actual intensity'].max() #determines the maximum actual intensity associated with an expected fragment ion
                
                #filter_xcorr_input['max real intensity'] = max_verified_intensity
                
                #filter_xcorr_input = filter_xcorr_input[filter_xcorr_input['Fragment actual intensity'] <= max_verified_intensity] #remove any peaks above the threshold that might be interfering peaks
                
                
                filter_xcorr_input = filter_xcorr_input.sort_values(by='Bin #')
                #filter_xcorr_input['Theoretical intensity'].replace(0, np.nan, inplace=True)
                #filter_xcorr_input['Fragment normalized intensity'].replace(0, np.nan, inplace=True)
                
                file_path = output_folder + '\\'+ l +'xcorr_input.csv'
                with open(file_path,'w',newline='') as filec:
                        writerc = csv.writer(filec)
                        filter_xcorr_input.to_csv(filec,index=False) 
                
                #correlation = filter_xcorr_input['Fragment normalized intensity'].corr(filter_xcorr_input['Theoretical intensity'])
                
                theoretical_array = filter_xcorr_input['Theoretical intensity'].values.tolist()
                actual_array = filter_xcorr_input['Fragment normalized intensity'].values.tolist() 
                
                correlation_array = signal.correlate(theoretical_array,actual_array,method='fft')
                correlation_norm = correlation_array.max()
                correlation_max = correlation_norm
                correlation_mean = correlation_array.mean()
                correlation = correlation_max-correlation_mean
                correlation_log =  np.log10(correlation_norm)
                #correlation = filter_xcorr_input['Theoretical intensity'].corr(filter_xcorr_input['Fragment normalized intensity'])
                corr_correlation.append(correlation)
                corr_pep.append(l)
                corr_scan.append(scan_to_report)
                correlation_max_store.append(correlation_max)
                correlation_mean_store.append(correlation_mean)
                log10_corr_store.append(correlation_log)
                #correlation_array_store.append(correlation_array)
                #end of XCorr calculation
### end of fragment AMM ###

### start of results compiling ###
final_report = pd.DataFrame()
final_report['Neuropeptide'] = final_report_storage_peptide
final_report['Scan'] = final_report_storage_scan
final_report['Sequence Coverage'] = final_report_storage_seq_coverage
final_report = final_report[final_report['Sequence Coverage'] > 0]

final_report = final_report.sort_values(by='Sequence Coverage',ascending=False)
final_report_seq = final_report.drop_duplicates(subset=['Neuropeptide'])
final_report_scan = final_report_seq.drop_duplicates(subset=['Scan'])

corr_report = pd.DataFrame()
corr_report['Neuropeptide'] = corr_pep
corr_report['Scan'] = corr_scan
#corr_report['Correlation array'] = correlation_array_store
corr_report['Correlation max'] = correlation_max_store
corr_report['Correlation mean'] = correlation_mean_store
corr_report['Correlation score'] = corr_correlation
corr_report['Log Correlation'] = log10_corr_store

final_report_w_corr = final_report_scan.merge(corr_report,on=['Neuropeptide','Scan'],how='left')
final_report_w_corr = final_report_w_corr.drop_duplicates()

file_path = output_folder + '\\final_report.csv'
with open(file_path,'w',newline='') as filec:
        writerc = csv.writer(filec)
        final_report_w_corr.to_csv(filec,index=False) 
        
### end of results compiling ###


end = time.time()
elapsed_time = (end - start)
param_log_export()
print('Analysis complete. Time elapsed:',(end - start),'s')

