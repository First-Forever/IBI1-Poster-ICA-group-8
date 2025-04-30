#Import necessary libraries
import math
import pandas as pd
import re
import sys
import matplotlib.pyplot as plt

# Add a dictionary called genetic_code to store all codons and their according amino acid
genetic_code = { 
# Phenylalanine
'UUU': 'Phenylalanine', 'UUC': 'Phenylalanine',
# Leucine
'UUA': 'Leucine', 'UUG': 'Leucine',
'CUU': 'Leucine', 'CUC': 'Leucine', 'CUA': 'Leucine', 'CUG': 'Leucine',
# Isoleucine
'AUU': 'Isoleucine', 'AUC': 'Isoleucine', 'AUA': 'Isoleucine',
# Methionine (Start)
'AUG': 'Methionine',
# Valine
'GUU': 'Valine', 'GUC': 'Valine', 'GUA': 'Valine', 'GUG': 'Valine',
# Serine
'UCU': 'Serine', 'UCC': 'Serine', 'UCA': 'Serine', 'UCG': 'Serine',
'AGU': 'Serine', 'AGC': 'Serine',
# Proline
'CCU': 'Proline', 'CCC': 'Proline', 'CCA': 'Proline', 'CCG': 'Proline',
# Threonine
'ACU': 'Threonine', 'ACC': 'Threonine', 'ACA': 'Threonine', 'ACG': 'Threonine',
# Alanine
'GCU': 'Alanine', 'GCC': 'Alanine', 'GCA': 'Alanine', 'GCG': 'Alanine',
# Tyrosine
'UAU': 'Tyrosine', 'UAC': 'Tyrosine',
# Histidine
'CAU': 'Histidine', 'CAC': 'Histidine',
# Glutamine
'CAA': 'Glutamine', 'CAG': 'Glutamine',
# Asparagine
'AAU': 'Asparagine', 'AAC': 'Asparagine',
# Lysine
'AAA': 'Lysine', 'AAG': 'Lysine',
# Aspartic acid
'GAU': 'Aspartic acid', 'GAC': 'Aspartic acid',
# Glutamic acid
'GAA': 'Glutamic acid', 'GAG': 'Glutamic acid',
# Cysteine
'UGU': 'Cysteine', 'UGC': 'Cysteine',
# Tryptophan
'UGG': 'Tryptophan',
# Arginine
'CGU': 'Arginine', 'CGC': 'Arginine', 'CGA': 'Arginine', 'CGG': 'Arginine',
'AGA': 'Arginine', 'AGG': 'Arginine',
# Glycine
'GGU': 'Glycine', 'GGC': 'Glycine', 'GGA': 'Glycine', 'GGG': 'Glycine',
# Stop codons
'UAA': 'Stop', 'UAG': 'Stop', 'UGA': 'Stop'
}

#Additional Function 1: Check if the sequence is valid
def check_sequence(seq):
    std_nucleotide = 'AGCU'         #Valid nucleotides
    #Check if the sequence starts with start codon AUG
    if seq[0:3] != 'AUG':                           #If the sequence does not start with AUG, print the error information and stop the program
        print("Start codon not valid! Try again!") 
        sys.exit()
    #Check if all nucleotides are valid
    for i in seq:
        #If invalid nucleotide is detected, print error information (with error nucleotide & location), then stop the program
        if i not in std_nucleotide:            
            print(f'The sequence has invalid nucleotide! Try again!\nInvalid nucleotide: {i}, location: {seq.find(i)}')
            sys.exit()

#Additional Function 2: Find the largest RSCU for each AA
def largest_RSCU_dict(df):
    AA_RSCU_max_dict = {}               #Initialize the dictionary to store all maximum RSCU for each amino acid
    for index in df.index:
        now_AA = df.loc[index, 'Amino acid']
        now_RSCU = df.loc[index, 'RSCU']        #Find the AA and RSCU on each row
        AA_RSCU_max_dict[now_AA] = max(now_RSCU, AA_RSCU_max_dict.get(now_AA, 0))
        #If the AA exists in the dictionary, compare its value to current RSCU, choose the larger one;
        #If not, store it into the dictionary
    return AA_RSCU_max_dict

#Additional Function 3: Calculate CAI value of the sequence:
def CAI_calculator(sequence):
    RSCU_max = largest_RSCU_dict(RSCU)
    L = len(sequence)
    CAI = 1
    for i in range (0, L, 3):
        codon = sequence[i: i+3]
        AA = genetic_code[codon]
        RSCU_codon = RSCU[RSCU['CODON'] == codon]['RSCU'][0]
        RSCU_AA = RSCU_max[AA]
        CAI *= RSCU_codon / RSCU_AA
        if codon in ('UAA', 'UAG', 'UGA'):
            L = i // 3
            break
    CAI = CAI ** (1/L)
    return CAI

#Additional function 4: count the number of every codon
def count_codon(seq):
    stop_codon = ('UAG', 'UAA', 'UGA')                              #Define the sequence of stop codon
    seq_len = len(seq)
    freq_dict = {}                                                  #Use dictionary to calculate all the frequency of codons (to save time)
    for i in range(0, seq_len, 3):                           
        now_seq = seq[i : i + 3]                                #Slice every three nucleotides
        if now_seq in stop_codon:                                   #If stop codon is met, break the loop
            break
        freq_dict[now_seq] = freq_dict.get(now_seq, 0) + 1                   #If the codon is in the dictionary, add 1 to the value; if not, add it as a key and give the value as 1
    return freq_dict
        
#Additional function 5: count the number of every AA
def count_AA(seq):
    freq_dict_codon = count_codon(seq)
    freq_dict_AA = {}
    for codon in freq_dict_codon.keys():
        AA = genetic_code[codon]
        freq_codon = freq_dict_codon[codon]
        freq_dict_AA[AA] = freq_dict_AA.get(AA, 0) + freq_codon
    return freq_dict_AA

#Function 1: most frequent trinucleotide
def most_frequent_trinucleotide(seq):
    freq_dict = count_codon(seq)
    mmax_value = 0                                                  #Initialize the variables to find the most frequent sequence
    mmax_key = []
    for key in freq_dict:
        if freq_dict[key] > mmax_value:                             #If the value is bigger than mmax_value, clear the answer list and append 
            mmax_key.clear()
            mmax_key.append(key)
            mmax_value = freq_dict[key]
        elif freq_dict[key] == mmax_value:                          #If they are the same, append the key into answer list
            mmax_key.append(key)
    return mmax_key, mmax_value

#Function 2: most frequent amino acid
def most_frequent_amino_acid(seq):
    stop_codon = ('UAG', 'UAA', 'UGA')
    codons = [] # Empty list called codons.
    for i in range(0, len(seq), 3): # For loop to iterate over the mRNA_sequence in steps of 3.
        codon = seq[i:i+3] # Extract each codon.
        if len(codon) == 3: # Check if the codon is of length 3.
            codons.append(codon) # Append the codon to the codons list.
        if codon in stop_codon:                                   #If stop codon is met, break the loop
            break
    aa_list = [genetic_code[codon] for codon in codons]
    aa_count = {}
    for aa in aa_list:
        aa_count[aa] = aa_count.get(aa, 0) + 1

    max_count = max(aa_count.values())
    most_frequent_amino_acid = [k for k, v in aa_count.items() if v == max_count]
    
    return most_frequent_amino_acid, max_count # Return the result.


#Function 3: Plot the frequency of amino acids:
def amino_acid_frequency_plot(seq):
    amino_acid_freq = count_AA(seq)
    plt.figure(figsize=(10, 6))
    plt.bar(x = amino_acid_freq.keys(), height = amino_acid_freq.values())
    plt.title('Frequency Distribution of Encoded Amino Acids')
    plt.xlabel('Amino Acid')
    plt.ylabel('Frequency')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()

""" file = open('', 'r')
skip_row = 0
for line in file:
    if 'CODON' in line:
        break
    skip_row += 1
RSCU = pd.read_table('', skiprows = skip_row) """
seq = ("AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUCACUGAUGAGUAU"
    "UUUCUUCACCCUAACUAUGGGAGCGAACUCCAUCUGGACUACAAUGCUGAGGUA"
    "UUGACUGCUGAAGUGGCACUUGGACCCUACUAUGGUAAGAAUACUGCCAGACCG"
    "CGCCUUCUUGACAUUAUGGUGACUGUACAGAAUGUGAUGGGUUCUCCUGACUAC"
    "AUGGUGGAGGCCAUUAAGGAAGCUGAUAACCCUAUUGAUUACUUGAUCAGACUG"
    "ACUGGUGCACUUGGUGUUGUCAUGACUGGUGCCCGUAAGUUCUUGAAAGACGGU"
    "GGUACUCCUCGUCCUAAGGAACUGACUGGUCAGCUGCUGAAUAAGCAUAAGACC"
    "AUCAGUCGCCUUAUUGACAAUAAGUAUGGUGACUUCACUGAUGAGUUCUAA")
check_sequence(seq)
print(most_frequent_trinucleotide(seq))
print(most_frequent_amino_acid(seq))
amino_acid_frequency_plot(seq)