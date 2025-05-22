# Import necessary libraries
import pandas as pd
import sys
import matplotlib.pyplot as plt
import numpy as np
import re

# Set necessary variables
species = input('Choose the species you want to check (human/yeast/E.coli)：')
L = 0                         

# Open the database required
file = open(str(species) + '_RSCU.tsv', 'r')            
skip_row = 0
for line in file:
    if 'CODON' in line:
        break
    skip_row += 1
RSCU = pd.read_table(str(species) + '_RSCU.tsv', skiprows = skip_row)
# Change T in the database to U
for i in range(len(RSCU['CODON'].tolist())):
    RSCU.loc[i,'CODON'] = re.sub('T', 'U', RSCU['CODON'][i]) # sub T to U

seq = ("AUGGCCAAGGUUACCGAUCAUCCUGAAGAGCUUCAGUUCUUCCAGAAGGCCCAGU"
"ACUUCGAGCAGAUCCUCAACAGUCGGACUGAGUUCUUGACCCGGCUGGAACAGUA"
"AGGAGGGUGA")


# Add a dictionary called genetic_code to store all codons and their according amino acid
genetic_code = { 
# Phenylalanine
'UUU': 'Phe', 'UUC': 'Phe',
# Leucine
'UUA': 'Leu', 'UUG': 'Leu',
'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
# Isoleucine
'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile',
# Methionine (Start)
'AUG': 'Met',
# Valine
'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
# Serine
'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser',
'AGU': 'Ser', 'AGC': 'Ser',
# Proline
'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
# Threonine
'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
# Alanine
'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
# Tyrosine
'UAU': 'Tyr', 'UAC': 'Tyr',
# Histidine
'CAU': 'His', 'CAC': 'His',
# Glutamine
'CAA': 'Gln', 'CAG': 'Gln',
# Asparagine
'AAU': 'Asn', 'AAC': 'Asn',
# Lysine
'AAA': 'Lys', 'AAG': 'Lys',
# Aspartic acid
'GAU': 'Asp', 'GAC': 'Asp',
# Glutamic acid
'GAA': 'Glu', 'GAG': 'Glu',
# Cysteine
'UGU': 'Cys', 'UGC': 'Cys',
# Tryptophan
'UGG': 'Trp',
# Arginine
'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
'AGA': 'Arg', 'AGG': 'Arg',
# Glycine
'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly',
# Stop codons
'UAA': 'Stop', 'UAG': 'Stop', 'UGA': 'Stop'
}

#Check if the sequence is valid
def check_sequence(seq):
    std_nucleotide = 'AGCU' 
    # Check for correct species
    if species != 'human' and species != 'yeast' and species != 'E.coli':
        print('Wrong species, please check your input.')
        sys.exit()        
    # Check if the sequence length is too short
    seq_len = len(seq)
    if seq_len < 3:
        print("The sequence length is too short! Try again!")
        sys.exit()
    # Check if the sequence starts with start codon AUG
    # If not, print the error information and stop the program
    if seq[0:3] != 'AUG':                           
        print("Start codon not valid! Try again!") 
        sys.exit()
    # Check if all nucleotides are valid
    for i in seq:
        # If invalid, print error information (with error nucleotide & location), then stop the program
        if i not in std_nucleotide:            
            print(f'The sequence has invalid nucleotide! Try again!\nInvalid nucleotide: {i}, location: {seq.find(i)}')
            sys.exit()
    # Check if the length is divisible by 3 (exclution: stop codon exists)        
    flag = True
    # If not, check if stop codon exists                             
    if seq_len % 3 != 0:                    
        flag = False
        stop_codon = ('UAG', 'UAA', 'UGA')
        for i in range(0, seq_len, 3):
            # If stop codon exists, the sequence is valid.
            if seq[i: i+3] in stop_codon:   
                flag = True
                break
    if not flag:
        print('The sequence length cannot be divided by 3! Try again!')
        sys.exit()

# Function 1: most frequent trinucleotide
# Count the number of every codon
stop_codon = ('UAG', 'UAA', 'UGA')                              
seq_len = len(seq)
freq_dict = {}                                                  
for i in range(0, seq_len, 3):
    #Slice every three nucleotides                           
    now_seq = seq[i : i + 3]
    # If the codon is in the dictionary, add 1 to the frequency
    # If not, add it as a key and give the frequency as 1                                
    freq_dict[now_seq] = freq_dict.get(now_seq, 0) + 1
    # Calculate number of codons (for additional function)                   
    L += 1
    # If meets stop codon, break the loop                     
    if now_seq in stop_codon:                                   
        break

def most_frequent_trinucleotide(seq):                                    
    mmax_value = 0                                                  
    mmax_key = []
    for key in freq_dict:
        # If the value > mmax_value, clear the answer list and append 
        if freq_dict[key] > mmax_value:                             
            mmax_key.clear()
            mmax_key.append(key)
            mmax_value = freq_dict[key]
        # If they are the same, append the key into answer list
        elif freq_dict[key] == mmax_value:                          
            mmax_key.append(key)
    return mmax_key, mmax_value

# Function 2: most frequent amino acid
def most_frequent_amino_acid(seq):
    most_frequent_codon, max_count = most_frequent_trinucleotide(seq)
    # Translate the codon to amino acid.
    most_frequent_AA = [genetic_code.get(codon, 'Unknown') for codon in most_frequent_codon] 
    return most_frequent_AA

#Function 3: Plot the frequency of amino acids:
# Initialize a dictionary to store AA counts
freq_dict_AA = {}
#Translate all codons to amino acid counts 
for codon, cnt in freq_dict.items():
    now_aa = genetic_code[codon]
    freq_dict_AA[now_aa] = freq_dict_AA.get(now_aa, 0) + cnt   

def amino_acid_frequency_plot(seq):
    # Extract the AA count data
    amino_acids = list(freq_dict_AA.keys())                       
    frequencies = list(freq_dict_AA.values())
    # Set the figure size
    plt.figure(figsize=(12, 6))                                         
    # Decorate the bar graph using viridis transformation
    num_bars = len(amino_acids)
    colors = plt.cm.viridis(np.linspace(0, 1, num_bars))
    
    # Plot the bar graph
    plt.bar(x = amino_acids, height = frequencies, color = colors, width = 0.8)  
    # Set the graph title, labels
    plt.title('Frequency Distribution of Encoded Amino Acids')
    plt.xlabel('Amino Acid')
    plt.ylabel('Frequency')
    # Rotate xticks for a more comfortable visualization
    plt.xticks(rotation=45)                 
    plt.tight_layout()
    # Show the graph
    plt.show()

# Additional Function: Calculate CAI value of the sequence
def largest_RSCU_dict(df):
    AA_RSCU_max_dict = {}               
    # Scan through the DataFrame to calculate the largest RSCU for each AA
    for index in df.index:
        #Find the AA and RSCU on each row              
        now_AA = df.loc[index, 'Amino acid']
        now_RSCU = df.loc[index, 'RSCU']        
        AA_RSCU_max_dict[now_AA] = max(now_RSCU, AA_RSCU_max_dict.get(now_AA, 0))
        #If the AA exists in the dictionary, hoose the larger RSCU for the AA;
        #If not, store it into the dictionary
    return AA_RSCU_max_dict

def CAI_calculator(sequence, freq_dict, L):
    RSCU_max_dic = largest_RSCU_dict(RSCU)                                        
    cai = 1
    # Remove number of AUG & UGG from L
    freq_dict['AUG'] = freq_dict.get('AUG', 0)
    freq_dict['UGG'] = freq_dict.get('UGG', 0)                                         
    L = L - freq_dict['AUG'] - freq_dict['UGG'] 
    # calculate CAI        
    for keys in freq_dict:
        if keys in RSCU['CODON'].tolist():
            AA = genetic_code[keys]
            RSCU_k = RSCU[RSCU['CODON'] == keys]['RSCU'].tolist()
            RSCU_max = RSCU_max_dic[AA]
            cai *= RSCU_k[0] / RSCU_max
        else:
            # If the codon is not in the RSCU database, assume it as 0.5
            cai = cai * 0.5  
    cai = cai ** (1/L)                          
    return cai

check_sequence(seq)
mmax_key, mmax_count = most_frequent_trinucleotide(seq)
print(f"Most frequent trinucleotide: {mmax_key}, max count: {mmax_count}")
print(f"Most frequent amino acid: {most_frequent_amino_acid(seq)}")
amino_acid_frequency_plot(seq)
print(f"CAI value for {species}: {CAI_calculator(seq, freq_dict, L).round(3)}")