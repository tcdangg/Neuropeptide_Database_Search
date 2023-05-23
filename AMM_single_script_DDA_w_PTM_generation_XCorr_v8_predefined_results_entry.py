# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 11:17:13 2023

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
from datetime import datetime
import scipy
from scipy.spatial import distance
import numpy as np
#import pysptools
#from distance_metrics_mcda import distance_metrics
start = time.time()

##User input##
output_parent_directory = r"C:\Users\lawashburn\Documents\DBpep_v2\finale\Correlation_Opt\v14_double_AMM_test" #folder in which all output directories will be generated
db_path = r"C:\Users\lawashburn\Documents\DBpep_v2\finale\Reference_DB\target_df_full.csv" #database fasta path
base_file_path = r"C:\Users\lawashburn\Documents\DBpep_v2\finale\Reference_DB"

precursor_error_cutoff = 50 #ppm
fragment_error_cutoff = 0.02 #Da
precursor_charges = [2,3,4,5,6,7,8]
fragment_charges = [1,2,3,4]
h_mass = 1.00784
bin_size = 5
number_of_steps = 7
min_seq_coverage = 25
normalization=50
subsequent_matching_rounds = 5
spectra_segments = 50

amidation = True
oxidation_M_status = True
pyroglu_E_status = True
pyroglu_Q_status = True
sulfo_Y_status = True
max_modifications = 2

filter_high_int_mz = False

# raw_converter_path_input = [r"C:/Users/lawashburn/Documents/DBpep_v2/XCorr_Opt/XCorr_validation/20230207/KD_Training_Spectra/MS2_formatted/Brain1_top10_formatted.txt",
#                             r"C:/Users/lawashburn/Documents/DBpep_v2/XCorr_Opt/XCorr_validation/20230207/KD_Training_Spectra/MS2_formatted/Brain1_top20_formatted.txt",
#                             r"C:/Users/lawashburn/Documents/DBpep_v2/XCorr_Opt/XCorr_validation/20230207/KD_Training_Spectra/MS2_formatted/Brain2_top10_formatted.txt",
#                             r"C:/Users/lawashburn/Documents/DBpep_v2/XCorr_Opt/XCorr_validation/20230207/KD_Training_Spectra/MS2_formatted/Brain2_top20_formatted.txt",
#                             r"C:/Users/lawashburn/Documents/DBpep_v2/XCorr_Opt/XCorr_validation/20230207/KD_Training_Spectra/MS2_formatted/Brain3_top10_formatted.txt",
#                             r"C:/Users/lawashburn/Documents/DBpep_v2/XCorr_Opt/XCorr_validation/20230207/KD_Training_Spectra/MS2_formatted/Brain3_top20_formatted.txt",
#                             r"C:/Users/lawashburn/Documents/DBpep_v2/XCorr_Opt/XCorr_validation/20230207/KD_Training_Spectra/MS2_formatted/CoG_Unlabeled_DDA_TR1_formatted.txt",
#                             r"C:/Users/lawashburn/Documents/DBpep_v2/XCorr_Opt/XCorr_validation/20230207/KD_Training_Spectra/MS2_formatted/CoG_Unlabeled_DDA_TR2_formatted.txt",
#                             r"C:/Users/lawashburn/Documents/DBpep_v2/XCorr_Opt/XCorr_validation/20230207/KD_Training_Spectra/MS2_formatted/PO_DDA_top10_TR1_180525095121_formatted.txt",
#                             r"C:/Users/lawashburn/Documents/DBpep_v2/XCorr_Opt/XCorr_validation/20230207/KD_Training_Spectra/MS2_formatted/PO1_top20_formatted.txt",
#                             r"C:/Users/lawashburn/Documents/DBpep_v2/XCorr_Opt/XCorr_validation/20230207/KD_Training_Spectra/MS2_formatted/PO2_top10_formatted.txt",
#                             r"C:/Users/lawashburn/Documents/DBpep_v2/XCorr_Opt/XCorr_validation/20230207/KD_Training_Spectra/MS2_formatted/PO2_top20_formatted.txt",
#                             r"C:/Users/lawashburn/Documents/DBpep_v2/XCorr_Opt/XCorr_validation/20230207/KD_Training_Spectra/MS2_formatted/PO3_top20_formatted.txt",
#                             r"C:/Users/lawashburn/Documents/DBpep_v2/XCorr_Opt/XCorr_validation/20230207/KD_Training_Spectra/MS2_formatted/SG_Unlabeled_DDA_TR1_formatted.txt",
#                             r"C:/Users/lawashburn/Documents/DBpep_v2/XCorr_Opt/XCorr_validation/20230207/KD_Training_Spectra/MS2_formatted/SG_Unlabeled_DDA_TR2_formatted.txt",
#                             r"C:/Users/lawashburn/Documents/DBpep_v2/XCorr_Opt/XCorr_validation/20230207/KD_Training_Spectra/MS2_formatted/SG_Unlabeled_DDA_TR3_formatted.txt"]

raw_converter_path_input = [r"C:\Users\lawashburn\Documents\DBpep_v2\finale\Reference_DB\perfect_spectra.txt"]


high_intensity_mz = [110.0717,120.0811,129.1025,136.0758,156.0769,157.1337,298.1035,407.1847,504.2565]

##Definition storage

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

## Start of theoretical spectra generator ##
def theoretical_spectra_generator(peptide_to_check):
    ###pull theoretical masses###
    plain_peptide = re.sub("[\(\[].*?[\)\]]", "",peptide_to_check) #removes any modification for mass calculations
    res_list_for_fragment = list_of_residues(peptide_to_check)
    mass_of_residues = []
    for residue in plain_peptide:
        residue_mass = aa_masses[residue]
        mass_of_residues.append(residue_mass)
    peptide_mass = (sum(mass_of_residues)) + check_termini(peptide_to_check) + check_PTM(peptide_to_check) #calculated MH mass
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
    ion_report = pd.concat([ion_report,b_ions_report])
    ion_report = pd.concat([ion_report,y_ions_report])
    ion_report = pd.concat([ion_report,b_ions_water_adduct])
    ion_report = pd.concat([ion_report,b_ions_ammonia_adduct])
    ion_report = pd.concat([ion_report,y_ions_ammonia_adduct])
    ion_report = pd.concat([ion_report,y_ions_water_adduct])
    ion_report = ion_report.drop_duplicates()
    
    ion_report = ion_report.rename(columns={'mass':'Fragment theoretical monoisotopic mass'})
    return ion_report

def scrambled(orig):
    dest = orig[:]
    random.shuffle(dest)
    return dest

def param_log_export(output_folder,sample_name):
    finish = time.time()
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
    elapsed_time = start-finish
    elapsed_time_report = 'Time elapsed : ' + str(elapsed_time) 
    
    param_file_entries = [out_path_report,rawconverter_path_report,db_path_report,sample_ID_report,p_err_report,f_err_report,max_p_charge_report,max_f_charge_report,
                          bin_size_report,no_steps_report,normalization_report,max_modifications_report,amidation_modifications_report,oxidation_M_status_modifications_report,
                          pyroglu_E_status_modifications_report,pyroglu_Q_status_modifications_report,sulfo_Y_status_modifications_report]

    param_file_path = output_folder + '\\' + 'parameter_file.txt'
    with open(param_file_path,'a') as f:
        f.writelines('\n'.join(param_file_entries))

def db_seq_mass_compile(fasta_to_df):
    final_seq_list = []
    
    complete_db = pd.DataFrame()
    
    for b in fasta_to_df:

        aas_present = []
        aas_present_index = []
        aas_absent = []
        aas_absent_index = []
        
        seq_len = len(b)
        for c in range(0,seq_len):
            aa_ID = b[c] #look at a single residue
            if aa_ID in modded_aas: #determine if residue is susceptible to modification
                aas_present.append(aa_ID)
                aas_present_index.append(str(c))
                aas_absent.append(aa_ID)
                aas_absent_index.append(str(c))
            else:
                aas_absent.append(aa_ID)
                aas_absent_index.append(str(c))

        number_mods = len(aas_present) #determine the number of AAs able to be modified

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
            merge_seq_log['Index'] = merge_seq_log['Index'].astype(int)
            ptm_applied_seq_log = merge_seq_log.sort_values(by='Index',ascending=True)
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
                
        ###start of mass calculations ###
        fasta_monoiso_mass = []
        pg_filtered_list = []

        for sequence in final_seq_list:
            if '(Pyro-glu)' in sequence:
                count_pg = sequence.count('(Pyro-glu')
                if count_pg == 1:
                    if sequence[1] == '(':
                        pg_filtered_list.append(sequence)
                    else:
                        pass
                else:
                        pass
            else:
                pg_filtered_list.append(sequence)
        finalized_mod_list = []  

        for seq in pg_filtered_list:
            if ')(' not in seq:
                finalized_mod_list.append(seq)
            else:
                pass
        for sequence in finalized_mod_list:
            mono_mass = monoisotopic_mass_calculator(sequence)
            fasta_monoiso_mass.append(mono_mass)
        
        db = pd.DataFrame()
        db['Sequence'] = finalized_mod_list
        db['Precursor theoretical monoisotopic mass'] = fasta_monoiso_mass

        
        complete_db = pd.concat([complete_db,db])
    complete_db = complete_db.drop_duplicates() 
    return complete_db
    
        ### end of convert fasta to dataframe with monoisotopic masses ###

def raw_file_detail_extraction(raw_file_path):
    raw_file_sample_name1 = raw_converter_path.replace(base_file_path,'')
    raw_file_sample_name2 = raw_file_sample_name1.replace('_formatted','')
    raw_file_sample_name3 = raw_file_sample_name2.replace('\\','')
    sample_name = raw_file_sample_name3.replace('.txt','')
    output_folder = output_parent_directory+'\\'+sample_name
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    ### generate output folder ###
    return sample_name, output_folder


def raw_file_data_extraction(raw_file_path):
    raw_converter = pd.read_csv(raw_converter_path, sep=",",skiprows=[0], names= ["m/z","resolution","charge","intensity", "MS2", "scan_number","precursor_charge",'null','Sample','Identifier','Iteration'])
    raw_converter = raw_converter.rename(columns={'m/z':'Fragment actual m/z',
                                                  'charge': 'Fragment actual charge',
                                                  'intensity':'Fragment actual intensity',
                                                  'MS2':'Precursor actual m/z',
                                                  "precursor_charge":'Precursor actual charge',
                                                  'scan_number':'Scan'})

    raw_converter['Fragment actual charge'] = raw_converter['Fragment actual charge'].replace(to_replace=0,value=1) #assume z=0 is z=1   
    exp_precursor = raw_converter.drop_duplicates() 
    exp_precursor = exp_precursor.copy()
    exp_precursor['Precursor actual monoisotopic mass'] =  ((exp_precursor['Precursor actual m/z']) * 
                                                            (exp_precursor['Precursor actual charge']))-(h_mass*(exp_precursor['Precursor actual charge']))
    print('yes')
    return exp_precursor

# def precursor_amm(db,exp_precursor):
#     #precursor_temp_cutoff = precursor_error_cutoff #rough estimate of ppm to Da to minimize extra search space
    
#     if len(db)<1: #throws an error if database file is empty
#         raise ValueError('Database file is empty')
    
#     db_sorted = db.sort_values(by = 'Precursor theoretical monoisotopic mass')
#     exp_precursor_sorted = exp_precursor.sort_values(by = 'Precursor actual monoisotopic mass')
    
#     scan_numbers = len(set(exp_precursor_sorted['Scan'].values.tolist()))
#     scans_check = set(exp_precursor_sorted['Scan'].values.tolist())
    
#     segement_size = int(scan_numbers/spectra_segments)
    
#     merge_match = pd.DataFrame()
    
#     scans_accounted = []
#     for width in range(1,spectra_segments):
#         start_index = (width-1)*segement_size
#         end_index = width*segement_size
        
#         scans_to_include = []
#         for xxx in range(start_index,end_index):
#             scans_to_include.append(xxx)
#             scans_accounted.append(xxx)

#         exp_precursor_segment = exp_precursor_sorted[exp_precursor_sorted['Scan'].isin(scans_to_include)]
    
#         merge_match2 = pd.merge_asof(exp_precursor_segment,db_sorted, left_on='Precursor actual monoisotopic mass', 
#                                     right_on='Precursor theoretical monoisotopic mass',
#                                     tolerance = precursor_error_cutoff, allow_exact_matches=True,direction='forward') 
#         merge_match3 = pd.merge_asof(exp_precursor_segment,db_sorted, left_on='Precursor actual monoisotopic mass', 
#                                     right_on='Precursor theoretical monoisotopic mass',
#                                     tolerance = precursor_error_cutoff, allow_exact_matches=True,direction='backward') 
        
#         merge_match = pd.concat([merge_match, merge_match2, merge_match3])
    
#     scans_missing_check = []
    
#     for xxxx in scans_check:
#         if xxxx not in scans_accounted:
#             scans_missing_check.append(xxxx)
#         else:
#             pass
    
#     if len(scans_missing_check)>0:
#         exp_precursor_segment = exp_precursor_sorted[exp_precursor_sorted['Scan'].isin(scans_missing_check)]

#         merge_match2 = pd.merge_asof(exp_precursor_segment,db_sorted, left_on='Precursor actual monoisotopic mass', 
#                                     right_on='Precursor theoretical monoisotopic mass',
#                                     tolerance = precursor_error_cutoff, allow_exact_matches=True,direction='forward') 
#         merge_match3 = pd.merge_asof(exp_precursor_segment,db_sorted, left_on='Precursor actual monoisotopic mass', 
#                                     right_on='Precursor theoretical monoisotopic mass',
#                                     tolerance = precursor_error_cutoff, allow_exact_matches=True,direction='backward') 
        
#         merge_match = pd.concat([merge_match, merge_match2, merge_match3])
#     else:
#         pass
    
#     merge_match= merge_match.drop_duplicates()

#     merge_match_filtered = merge_match.dropna(subset=['Sequence','Precursor theoretical monoisotopic mass'])
#     merge_match_filtered = merge_match_filtered.copy() #silences pesky warning
#     merge_match_filtered['Precursor error (ppm)'] = ((abs((merge_match_filtered['Precursor theoretical monoisotopic mass'])-
#                                                           (merge_match_filtered['Precursor actual monoisotopic mass'])))/
#                                                       (merge_match_filtered['Precursor theoretical monoisotopic mass'])) * 1E6

#     precursor_amm_results = merge_match_filtered[merge_match_filtered['Precursor error (ppm)'] <= precursor_error_cutoff]


#     return precursor_amm_results

def precursor_amm(db,exp_precursor):
    precursor_temp_cutoff = precursor_error_cutoff*3 #rough estimate of ppm to Da to minimize extra search space
    
    if len(db)<1: #throws an error if database file is empty
        raise ValueError('Database file is empty')
    
    db_sorted = db.sort_values(by = 'Precursor theoretical monoisotopic mass') #sort database monoisotopic mass and experimental for mergeasof
    exp_precursor_sorted = exp_precursor.sort_values(by = 'Precursor actual monoisotopic mass')
    
    scan_numbers = len(set(exp_precursor_sorted['Scan'].values.tolist()))
    scans_check = set(exp_precursor_sorted['Scan'].values.tolist())
    
    segement_size = int(scan_numbers/spectra_segments) #to improve results, amm is done in fragments of spectra
    
    merge_match = pd.DataFrame()
    
    scans_accounted = []
    for width in range(1,spectra_segments):
        start_index = (width-1)*segement_size
        end_index = width*segement_size
        
        scans_to_include = []
        for xxx in range(start_index,end_index):
            scans_to_include.append(xxx)
            scans_accounted.append(xxx)

        exp_precursor_segment = exp_precursor_sorted[exp_precursor_sorted['Scan'].isin(scans_to_include)]
    
        merge_match2 = pd.merge_asof(exp_precursor_segment,db_sorted, left_on='Precursor actual monoisotopic mass', 
                                    right_on='Precursor theoretical monoisotopic mass',
                                    tolerance = precursor_error_cutoff, allow_exact_matches=True,direction='forward') 
        merge_match3 = pd.merge_asof(exp_precursor_segment,db_sorted, left_on='Precursor actual monoisotopic mass', 
                                    right_on='Precursor theoretical monoisotopic mass',
                                    tolerance = precursor_error_cutoff, allow_exact_matches=True,direction='backward') 
        
        merge_match = pd.concat([merge_match, merge_match2, merge_match3])
    
    scans_missing_check = []
    
    for xxxx in scans_check:
        if xxxx not in scans_accounted:
            scans_missing_check.append(xxxx)
        else:
            pass
    
    if len(scans_missing_check)>0:
        exp_precursor_segment = exp_precursor_sorted[exp_precursor_sorted['Scan'].isin(scans_missing_check)]

        merge_match2 = pd.merge_asof(exp_precursor_segment,db_sorted, left_on='Precursor actual monoisotopic mass', 
                                    right_on='Precursor theoretical monoisotopic mass',
                                    tolerance = precursor_error_cutoff, allow_exact_matches=True,direction='forward') 
        merge_match3 = pd.merge_asof(exp_precursor_segment,db_sorted, left_on='Precursor actual monoisotopic mass', 
                                    right_on='Precursor theoretical monoisotopic mass',
                                    tolerance = precursor_error_cutoff, allow_exact_matches=True,direction='backward') 
        
        merge_match = pd.concat([merge_match, merge_match2, merge_match3])
    else:
        pass
    
    merge_match= merge_match.drop_duplicates()

    merge_match_filtered = merge_match.dropna(subset=['Sequence','Precursor theoretical monoisotopic mass'])
    merge_match_filtered = merge_match_filtered.copy() #silences pesky warning
    merge_match_filtered['Precursor error (ppm)'] = ((abs((merge_match_filtered['Precursor theoretical monoisotopic mass'])-
                                                          (merge_match_filtered['Precursor actual monoisotopic mass'])))/
                                                      (merge_match_filtered['Precursor theoretical monoisotopic mass'])) * 1E6

    precursor_amm_results = merge_match_filtered[merge_match_filtered['Precursor error (ppm)'] <= precursor_error_cutoff]


    return precursor_amm_results

def prelim_amm_candidate_seqs(secondary_amm_prep_clean):
    candidate_sequences_raw = secondary_amm_prep_clean['Sequence'].values.tolist()   
    candidate_sequences = []
    for m in candidate_sequences_raw:
        if m not in candidate_sequences:
            candidate_sequences.append(m)
    return candidate_sequences

def seq_coverage_calc(merge_fragment_match_filtered,ion_report,scan_to_report,peptide):
    if len(merge_fragment_match_filtered)>0:
        merge_fragment_match_filtered['ion_format'] = merge_fragment_match_filtered['ion'].str.extract('(\d+)', expand=False)
        ion_report['ion_format'] = ion_report['ion'].str.extract('(\d+)', expand=False)
        merge_fragment_match_filtered_ion = merge_fragment_match_filtered.drop_duplicates(subset=['ion_format'],keep='first').reset_index(drop=True)
        ion_report_filtered_ion = ion_report.drop_duplicates(subset=['ion_format'],keep='first').reset_index(drop=True)
        seq_coverage = (len(merge_fragment_match_filtered_ion)/len(ion_report_filtered_ion))*100
        
        seq_data = {'Sequence':[peptide],
                    'Scan':[scan_to_report],
                    'Sequence coverage':[seq_coverage]}
        
        seq_coverage_rep = pd.DataFrame(seq_data)
        return seq_coverage_rep
    else:
        seq_data = {'Sequence':[peptide],
                    'Scan':[scan_to_report],
                    'Sequence coverage':[0]}
        
        seq_coverage_rep = pd.DataFrame(seq_data)
        return seq_coverage_rep


def fragment_amm(secondary_amm_prep_clean,peptide,ion_report,peptide_report_output,exp_precursor):
    
    
    sequence_coverage_rep_final = pd.DataFrame()
    
    filtered_secondary_amm = secondary_amm_prep_clean[secondary_amm_prep_clean['Sequence'] == peptide] #filter amm from before for the sequence we are looking at 
    filtered_secondary_amm = filtered_secondary_amm.copy() #silences pesky warning
    filtered_secondary_amm['Fragment actual monoisotopic mass'] = (filtered_secondary_amm['Fragment actual m/z'] * 
                                                                   filtered_secondary_amm['Fragment actual charge']) - (h_mass*filtered_secondary_amm['Fragment actual charge'])         
    
    scans_present_raw = filtered_secondary_amm['Scan'].values.tolist()
    scans_present = []
    for z in scans_present_raw:
        if z not in scans_present:
            scans_present.append(z)
    for scan_to_report in scans_present:
        
        peptide_rep_output_folder = peptide_report_output+'\\fragment_matches'
        if not os.path.exists(peptide_rep_output_folder):
            os.makedirs(peptide_rep_output_folder)
        
        scans_filtered_secondary_amm = filtered_secondary_amm[filtered_secondary_amm['Scan'] == scan_to_report] #filter to look at just one scan
        scans_filtered_secondary_amm = scans_filtered_secondary_amm.sort_values(by='Fragment actual monoisotopic mass')

        ion_report = ion_report.sort_values(by='Fragment theoretical monoisotopic mass')
        prelim_fragment_matches = pd.merge_asof(scans_filtered_secondary_amm,ion_report,left_on='Fragment actual monoisotopic mass',
                                                right_on='Fragment theoretical monoisotopic mass', tolerance=fragment_error_cutoff,allow_exact_matches=True, direction='nearest')
        merge_fragment_match_filtered = prelim_fragment_matches.dropna(subset=['ion','Fragment theoretical monoisotopic mass'])
        merge_fragment_match_filtered = merge_fragment_match_filtered.copy()
        merge_fragment_match_filtered['Fragment error (Da)'] = merge_fragment_match_filtered['Fragment actual monoisotopic mass'] - merge_fragment_match_filtered['Fragment theoretical monoisotopic mass']

        merge_fragment_match_filtered = merge_fragment_match_filtered[merge_fragment_match_filtered['Fragment error (Da)'] <= fragment_error_cutoff]
        
        if len(merge_fragment_match_filtered)>0: #only export fragment report if fragments are found
        
            output_path_rep = peptide_rep_output_folder + '\\' + peptide + '_' + str(scan_to_report) + '_fragment_report.csv'
    
            with open(output_path_rep,'w',newline='') as filec:
                    writerc = csv.writer(filec)
                    merge_fragment_match_filtered.to_csv(filec,index=False) 
        else:
            pass
        seq_coverage = seq_coverage_calc(merge_fragment_match_filtered,ion_report,scan_to_report,peptide)
        sequence_coverage_rep_final = pd.concat([sequence_coverage_rep_final,seq_coverage])
    return sequence_coverage_rep_final

        
        
        
        ###theoretical spectra for XCorr
def xcorr_calc(ion_report,exp_precursor,scan_to_report,peptide_report_output,peptide):        
    
    xcorr_rep_output_folder = peptide_report_output+'\\xcorr_data'
    if not os.path.exists(xcorr_rep_output_folder):
        os.makedirs(xcorr_rep_output_folder)
    
    actual_peak_intensities = exp_precursor[exp_precursor['Scan'] == scan_to_report]
    theoretical_mass_list_z1 = ion_report['Fragment theoretical monoisotopic mass'].values.tolist()        
    theoretical_mass_list_z_all = []        
    for zz in theoretical_mass_list_z1:
        for yy in fragment_charges:
            new_mz = (zz + (h_mass*yy))/yy
            theoretical_mass_list_z_all.append(new_mz)
    theo_mass_list = []
    theo_mass_list_intensity = []
    
    for dd in theoretical_mass_list_z_all:
        theo_mass_list.append(dd)
        theo_mass_list_intensity.append(50)
        
        theo_minus_flank1 = dd - bin_size
        theo_plus_flank1 = dd + bin_size
        
        theo_mass_list.append(theo_minus_flank1)
        theo_mass_list_intensity.append(25)
        
        theo_mass_list.append(theo_plus_flank1)
        theo_mass_list_intensity.append(25)
        
        theo_minus_flank2 = theo_minus_flank1 - bin_size
        theo_mass_list.append(theo_minus_flank2)
        theo_mass_list_intensity.append(10)
        
        theo_plus_flank2 = theo_plus_flank1 + bin_size
        theo_mass_list.append(theo_plus_flank2)
        theo_mass_list_intensity.append(10)
        
    theoretical_spectra_w_flanks = pd.DataFrame()
    theoretical_spectra_w_flanks['theoretical m/z'] = theo_mass_list
    theoretical_spectra_w_flanks['theoretical m/z'] = theoretical_spectra_w_flanks['theoretical m/z'].round(4)
    theoretical_spectra_w_flanks['theoretical intensity'] = theo_mass_list_intensity

    output_path_rep = xcorr_rep_output_folder + '\\' + peptide + '_' + str(scan_to_report) + '_theo_rep.csv'
    with open(output_path_rep,'w',newline='') as filec:
            writerc = csv.writer(filec)
            theoretical_spectra_w_flanks.to_csv(filec,index=False) 
    
    experimental_spectra = pd.DataFrame()
    experimental_spectra['experimental m/z'] = actual_peak_intensities['Fragment actual m/z']
    experimental_spectra['experimental non-normalized intensity'] = actual_peak_intensities['Fragment actual intensity']
    experimental_spectra = experimental_spectra.sort_values(by='experimental m/z')
    
    
    filtered_experimental = pd.DataFrame()
    
    if filter_high_int_mz == True:
        banned_mz = pd.DataFrame()
        banned_mz['Banned m/z'] = high_intensity_mz
        banned_mz = banned_mz.sort_values(by='Banned m/z')       
        filtered_experimental_prelim = pd.merge_asof(experimental_spectra,banned_mz,left_on='experimental m/z',
                                                right_on='Banned m/z', tolerance=0.02,allow_exact_matches=True, direction='nearest')        
        filtered_experimental_prelim['Banned m/z'] = filtered_experimental_prelim['Banned m/z'].replace(np.nan, 0)        
        filtered_experimental_prelim = filtered_experimental_prelim[filtered_experimental_prelim['Banned m/z'] == 0]
        filtered_experimental = pd.concat([filtered_experimental,filtered_experimental_prelim])
    
    else:
        filtered_experimental = pd.concat([filtered_experimental,experimental_spectra])
    filtered_experimental = filtered_experimental.replace(np.nan, 0)
    
    

    exp_max_int = filtered_experimental['experimental non-normalized intensity'].max()    
    
    filtered_experimental['experimental intensity'] = 50
    
    output_path_rep = xcorr_rep_output_folder + '\\' + peptide + '_' + str(scan_to_report) + '_exp_rep.csv'
    with open(output_path_rep,'w',newline='') as filec:
            writerc = csv.writer(filec)
            filtered_experimental.to_csv(filec,index=False) 
    
    range_eval = []    
    range_eval.append(float(theoretical_spectra_w_flanks['theoretical m/z'].min()))
    range_eval.append(float(theoretical_spectra_w_flanks['theoretical m/z'].max()))
    range_eval.append(float(filtered_experimental['experimental m/z'].max()))
    range_eval.append(float(filtered_experimental['experimental m/z'].min()))
    xcorr_start_whole = (min(range_eval)) - 10
    xcorr_end_whole = (max(range_eval)) + 10
    xcorr_start = int(xcorr_start_whole)
    xcorr_end = int(xcorr_end_whole)
    
    interval_bins = []
    
    for aa in np.arange(xcorr_start,xcorr_end,bin_size):
        interval_bins.append(float(aa))
    
    bin_tolerance = bin_size
    
    xcorr_array = pd.DataFrame()
    xcorr_array['Bin #'] = interval_bins
    #xcorr_array['Bin #'] = xcorr_array['Bin #'].round(2)
    
    xcorr_array = xcorr_array.sort_values(by='Bin #')
    theoretical_spectra_w_flanks = theoretical_spectra_w_flanks.sort_values(by='theoretical m/z')
    filtered_experimental = filtered_experimental.sort_values(by='experimental m/z')
    
    xcorr_array = pd.merge_asof(xcorr_array,filtered_experimental,left_on='Bin #',
                                            right_on='experimental m/z', tolerance=bin_tolerance,allow_exact_matches=True, direction='nearest')
    
    xcorr_array = pd.merge_asof(xcorr_array,theoretical_spectra_w_flanks,left_on='Bin #',
                                            right_on='theoretical m/z', tolerance=bin_tolerance,allow_exact_matches=True, direction='nearest')
    
    xcorr_array = xcorr_array.replace(np.nan, 0)

    theoretical_array = xcorr_array['experimental intensity'].values.tolist()
    experimental_array = xcorr_array['theoretical intensity'].values.tolist()
    correlation_score = np.dot(theoretical_array,experimental_array)

    return (correlation_score)

def perform_db_search(fasta_w_mass,exp_precursor,peptide_report_output):
    no_ban_seq = []
    fasta_w_mass_monitor = []
    fasta_w_mass_monitor2 = []
    results_table = pd.DataFrame()
    
    
    precursor_amm_results = precursor_amm(fasta_w_mass,exp_precursor)
    precursor_candidate_seqs = prelim_amm_candidate_seqs(precursor_amm_results)
    precursor_amm_results = precursor_amm_results.drop_duplicates()    
    file_path = peptide_report_output + '\\precursor_AMM_results.csv'
    with open(file_path,'w',newline='') as filec:
            writerc = csv.writer(filec)
            precursor_amm_results.to_csv(filec,index=False)         
        
    for peptide in precursor_candidate_seqs:
        ion_report = theoretical_spectra_generator(peptide)
        fragment_amm_res = fragment_amm(precursor_amm_results,peptide,ion_report,peptide_report_output,exp_precursor)
        scan = (fragment_amm_res['Scan'].values.tolist())
        for ss in scan:
            xcorr_res = xcorr_calc(ion_report,exp_precursor,ss,peptide_report_output,peptide)
            fragment_amm_res['Correlation value'] = xcorr_res
            results_table = pd.concat([results_table,fragment_amm_res])

    number_IDs = results_table['Sequence'].values.tolist()
    number_IDs = list(set(number_IDs))
        
    
    file_path = peptide_report_output + '\\final_report.csv'
    with open(file_path,'w',newline='') as filec:
            writerc = csv.writer(filec)
            results_table.to_csv(filec,index=False)    
    
    results_table = results_table.drop_duplicates()
    results_table['Scan count'] = results_table['Scan'].map(results_table.groupby('Scan')['Scan'].count())

    finalized_psms = results_table[results_table['Scan count'] == 1]
    candidate_psms = results_table[results_table['Scan count'] > 1]

    if len(candidate_psms)>0:
        scans_to_check = set(candidate_psms['Scan'].values.tolist())
        for bbb in scans_to_check:
            scan_zoom_in = candidate_psms[candidate_psms['Scan'] == bbb]
            ###filter 1, keep only max coverage
            max_coverage_scan = scan_zoom_in['Sequence coverage'].max()
            scan_filter1 = scan_zoom_in[scan_zoom_in['Sequence coverage'] == max_coverage_scan]
            if len(scan_filter1) == 1:
                finalized_psms = pd.concat([finalized_psms,scan_filter1])
            if len(scan_filter1) > 1:
                #filter 2, keep only unique peptides from finalized list
                scan_filter1 = scan_filter1.copy()
                scan_filter1['Seq count'] = scan_filter1['Sequence'].map(finalized_psms.groupby('Sequence')['Sequence'].count())
                #scan_filter2 = scan_filter1['Seq count'].isnull()
                scan_filter2 = scan_filter1[scan_filter1['Seq count'].isnull()]
                if len(scan_filter2) == 1:
                    scan_filter2 = scan_filter2.drop(columns=['Seq count'])
                    finalized_psms = pd.concat([finalized_psms,scan_filter2])
                if len(scan_filter2) > 1:
                    #filter 3, choose randomly
                    random_selection = scan_filter2.sample()
                    random_selection = random_selection.drop(columns=['Seq count'])
                    finalized_psms = pd.concat([finalized_psms,random_selection])
                else:
                    pass
            else:
                pass
    else:
        pass

    finalized_psms = finalized_psms[finalized_psms['Sequence'].notnull()]              
    number_seqs = set(finalized_psms['Sequence'].values.tolist())  
    number_scans = set(finalized_psms['Scan'].values.tolist())  

    file_path = peptide_report_output + '\\first_psm_report_check.csv'
    with open(file_path,'w',newline='') as filec:
            writerc = csv.writer(filec)
            finalized_psms.to_csv(filec,index=False)
      
    used_scans =  finalized_psms['Scan'].values.tolist()
    all_candidate_scans = results_table['Scan'].values.tolist()

    available_scans = []
    for jj in all_candidate_scans:
        if jj not in used_scans:
            if jj not in available_scans:
                available_scans.append(jj)

    for scan in available_scans:
        scan_check1 = results_table[results_table['Scan'] == scan]
        if len(scan_check1) == 1:
            finalized_psms = pd.concat([finalized_psms,scan_check1])
        if len(scan_check1) > 1:
            min_scanno = scan_check1['Scan count'].min()
            max_seq_cov = scan_check1['Sequence coverage'].max()

            seq_filtered_scan_check1 = scan_check1[scan_check1['Sequence coverage'] == max_seq_cov]
            if len(seq_filtered_scan_check1) == 1:
                finalized_psms = pd.concat([finalized_psms,seq_filtered_scan_check1])
            if len(seq_filtered_scan_check1) > 1:
                min_scanno = seq_filtered_scan_check1['Scan count'].min()
                scan_no_seq_filtered_scan_check1 = seq_filtered_scan_check1[seq_filtered_scan_check1['Scan count'] == min_scanno]
                if len(scan_no_seq_filtered_scan_check1) == 1:
                    finalized_psms = pd.concat([finalized_psms,scan_no_seq_filtered_scan_check1])
                if len(scan_no_seq_filtered_scan_check1) > 1:
                    seq_taken = scan_no_seq_filtered_scan_check1['Sequence'].values.tolist()
                    seq_scan_no_seq_filtered_scan_check1 = scan_no_seq_filtered_scan_check1[~scan_no_seq_filtered_scan_check1['Sequence'].isin(seq_taken)]
                    if len(seq_scan_no_seq_filtered_scan_check1) == 1:
                        finalized_psms = pd.concat([finalized_psms,seq_scan_no_seq_filtered_scan_check1])
                    if len(seq_scan_no_seq_filtered_scan_check1) != 1:
                        rand_seq_scan_no_seq_filtered_scan_check1 = scan_no_seq_filtered_scan_check1.sample(n=1, random_state=1)
                        finalized_psms = pd.concat([finalized_psms,rand_seq_scan_no_seq_filtered_scan_check1])
    return finalized_psms
###################Commands#########################

fasta_w_mass = pd.read_csv(db_path)
fasta_w_mass['Precursor theoretical monoisotopic mass'].round(decimals = 4)



for raw_converter_path in raw_converter_path_input:
    rounds_number = [1]
    unique_IDS_number = []
    
    details = raw_file_detail_extraction(raw_converter_path)
    sample_name = details[0]
    peptide_report_output = details[1]
    exp_precursor = raw_file_data_extraction(raw_converter_path)
    
    finalized_psms = perform_db_search(fasta_w_mass,exp_precursor,peptide_report_output)

    file_path = peptide_report_output + '\\final_psm_report.csv'
    with open(file_path,'w',newline='') as filec:
            writerc = csv.writer(filec)
            finalized_psms.to_csv(filec,index=False)
    unique_IDS = len(set(finalized_psms['Sequence'].values.tolist()))
    unique_IDS_number.append(unique_IDS)
    rounds = 3
    
    for a in range(0,rounds):
        rounds_number.append(a+2)
        all_final_IDs = finalized_psms['Sequence'].values.tolist()
        fasta_w_mass_filtered = fasta_w_mass[~fasta_w_mass['Sequence'].isin(all_final_IDs)]
        if len(fasta_w_mass_filtered)>1:
            finalized_psms_sorted = finalized_psms.sort_values(by=['Sequence coverage'],ascending=False)
            finalized_psms_no_dubs = finalized_psms_sorted.drop_duplicates(subset=['Sequence'])
            psm_final_taken_scans = finalized_psms_no_dubs['Scan'].values.tolist()
            exp_precursor_filtered = exp_precursor[~exp_precursor['Scan'].isin(psm_final_taken_scans)]
            if len(exp_precursor_filtered)>1:
                finalized_psms_round2 = perform_db_search(fasta_w_mass_filtered,exp_precursor_filtered,peptide_report_output)
                finalized_psms = pd.concat([finalized_psms,finalized_psms_round2])
        
        finalized_psms = finalized_psms.sort_values(by='Sequence coverage',ascending=False)
        finalized_psms = finalized_psms.drop_duplicates(subset='Scan')
        unique_IDS = len(set(finalized_psms['Sequence'].values.tolist()))
        unique_IDS_number.append(unique_IDS)
        
        file_path = peptide_report_output + '\\final_psm_report'+str(a+1)+'.csv'
        with open(file_path,'w',newline='') as filec:
                writerc = csv.writer(filec)
                finalized_psms.to_csv(filec,index=False)
        
    rounds_report = pd.DataFrame()
    rounds_report['Round #'] = rounds_number
    rounds_report['# Unique IDs'] = unique_IDS_number
    
    file_path = peptide_report_output + '\\rounds_report.csv'
    with open(file_path,'w',newline='') as filec:
            writerc = csv.writer(filec)
            rounds_report.to_csv(filec,index=False)

    param_log_export(peptide_report_output,sample_name)
    