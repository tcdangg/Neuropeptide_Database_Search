import pandas as pd
import random
import collections
import csv

#from AMM_single_script_DDA_w_PTM_generation_v2 import db
#database = db['Sequence']

###
database_csv = pd.read_csv(r"C:\Users\lawashburn\Desktop\ALC50_Mass_Search_Files\duplicate_removed_crustacean_database_validated_formatted20220725.csv")
unmodified_database_list = database_csv['Sequence'].to_list()
modified_db_list = []

for db in unmodified_database_list:
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
    
def shuffled_decoy_database():
    """
    ### Residues in each target sequences are randomly shuffled and generated as the decoy sequences
    
    Return:
        shuffled_decoy_list (list): Shuffle decoy database with Identify Threshold applied
    """
    shuffled_decoy_list = []
    

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
    return target_decoy_df

def decoy_database_iterations(n):
    """
    ### regenerates decoy databases based on the input number
    
    Args:
        n (int): number of decoy database iteratures
    Return:
        concatenated_decoy_db (list): list contains combined decoy databases
    """
    concatenated_decoy_db = []
    for i in range(n):
        decoy_database = shuffled_decoy_database()
        for decoy_sequence in decoy_database:
            if decoy_sequence in concatenated_decoy_db:
                if len(decoy_sequence) <= 3:
                    pass
                else: 
                    decoy_sequence = shuffled_algorithm(decoy_sequence)
            concatenated_decoy_db.append(decoy_sequence)
    return concatenated_decoy_db

first_it = decoy_database_iterations(10)
decoy_df = pd.DataFrame()
decoy_df['Sequence'] = first_it
decoy_df['Status'] = 'decoy'

target_df = pd.DataFrame()
target_df['Sequence'] = database
target_df['Status'] = 'target'

database_out = pd.concat([decoy_df,target_df])
database_out = database_out.sample(frac=1)

file_path = r"C:\Users\lawashburn\Documents\DBpep_v2\finale\20230414\target_decoy_database_iteration10.csv"
with open(file_path,'w',newline='') as filec:
        writerc = csv.writer(filec)
        database_out.to_csv(filec,index=False) 
