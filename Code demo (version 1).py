#Import necessary libraries
import pandas as pd
import sys
import matplotlib.pyplot as plt
import numpy as np

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
    #Check if the sequence length is too short
    seq_len = len(seq)
    if seq_len < 3:
        print("The sequence length is too short! Try again!")
        sys.exit()
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
    #Check if the length is divisible by 3 (exclution: stop codon exists)        
    flag = True                             #A boolean flag to decide validity
    if seq_len % 3 != 0:                    #If the length is not divisible by 3, check if stop codon exists
        flag = False
        stop_codon = ('UAG', 'UAA', 'UGA')
        for i in range(0, seq_len, 3):
            if seq[i: i+3] in stop_codon:   #If stop codon exists, the sequence is valid.
                flag = True
                break
    if not flag:
        print('The sequence length cannot be divided by 3! Try again!')
        sys.exit()

#Additional Function 2: Find the largest RSCU for each AA
def largest_RSCU_dict(df):
    AA_RSCU_max_dict = {}               #Initialize the dictionary to store all maximum RSCU for each amino acid
    #Scan through all elements in the DataFrame to calculate the largest RSCU for each AA
    for index in df.index:              
        now_AA = df.loc[index, 'Amino acid']
        now_RSCU = df.loc[index, 'RSCU']        #Find the AA and RSCU on each row
        AA_RSCU_max_dict[now_AA] = max(now_RSCU, AA_RSCU_max_dict.get(now_AA, 0))
        #If the AA exists in the dictionary, compare its value to current RSCU, choose the larger one;
        #If not, store it into the dictionary
    return AA_RSCU_max_dict

#Additional Function 3: Calculate CAI value of the sequence:
def CAI_calculator(sequence):
    RSCU_max = largest_RSCU_dict(RSCU)              #First give the dictionary that stores the maximum RSCU for each AA
    L = 0                                           #length of the peptide sequence (including 'stop')
    CAI = 1                                         #Cumulative production of CAI
    for i in range (0, len(sequence), 3):           #Scan through the sequence
        #Extract the codon and according AA
        codon = sequence[i: i+3]
        AA = genetic_code[codon]
        #Find the RSCU for each codon
        matches = RSCU[RSCU['CODON'] == codon]
        if matches.empty:                                 #Check if the codons all exist in RSCU table
            print(f"Codon {codon} in the sequence is not found in RSCU table! Try again!")
            sys.exit()
        else:
            RSCU_codon = matches.loc[0, 'RSCU']
        RSCU_AA = RSCU_max[AA]                      #Find the RSCU for the according AA
        CAI *= RSCU_codon / RSCU_AA                 #Calculate the CAI value
        L += 1
        if codon in ('UAA', 'UAG', 'UGA'):          #If stop codon exists, break the loop
            break
    CAI = CAI ** (1/L)                          #Final CAI value
    return CAI

#Additional function 4: count the number of every codon
def count_codon(seq):
    stop_codon = ('UAG', 'UAA', 'UGA')                              #Define the sequence of stop codon
    seq_len = len(seq)
    freq_dict = {}                                                  #Use dictionary to calculate all the frequency of codons (to save time)
    for i in range(0, seq_len, 3):                           
        now_seq = seq[i : i + 3]                                #Slice every three nucleotides
        freq_dict[now_seq] = freq_dict.get(now_seq, 0) + 1                   #If the codon is in the dictionary, add 1 to the value; if not, add it as a key and give the value as 1
        if now_seq in stop_codon:                                   #If stop codon is met, break the loop
            break
    return freq_dict
        
#Additional function 5: count the number of every AA
def count_AA(seq):
    freq_dict_codon = count_codon(seq)                      #Find the codon count of the sequence
    freq_dict_AA = {}                                       #Initialize a dictionary to store AA counts
    for codon in freq_dict_codon.keys():                    #Translate all codons to amino acid counts
        AA = genetic_code[codon]
        freq_codon = freq_dict_codon[codon]
        if AA == 'Stop':                                    #If stop codon is met, break the loop
            break
        freq_dict_AA[AA] = freq_dict_AA.get(AA, 0) + freq_codon #Add translated codon count into freq_dict_AA
    return freq_dict_AA

#Function 1: most frequent trinucleotide
def most_frequent_trinucleotide(seq):
    freq_dict = count_codon(seq)                                    #Count all codons
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
        if codon in stop_codon:                                   #If stop codon is met, break the loop
            break
        if len(codon) == 3: # Check if the codon is of length 3.
            codons.append(codon) # Append the codon to the codons list.
    aa_list = [genetic_code[codon] for codon in codons]           #Store amino acid sequence after translation
    aa_count = {}                                                 #Initialize a dictionary to count amino acids
    for aa in aa_list:                                            #Count all amino acids
        aa_count[aa] = aa_count.get(aa, 0) + 1

    max_count = max(aa_count.values())                            #Find the largest count number
    most_frequent_amino_acid = [k for k, v in aa_count.items() if v == max_count]   #Find all amino acids whose count number coordinate with the max_count
    
    return most_frequent_amino_acid, max_count # Return the result.


#Function 3: Plot the frequency of amino acids:
def amino_acid_frequency_plot(seq):
    amino_acid_freq = count_AA(seq)                                     #Find the AA count
    amino_acids = list(amino_acid_freq.keys())                          #Extract the AA count data
    frequencies = list(amino_acid_freq.values())

    plt.figure(figsize=(12, 6))                                         #Set the figure size
    #Decorate the bar graph using viridis transformation
    num_bars = len(amino_acids)
    colors = plt.cm.viridis(np.linspace(0, 1, num_bars))
    plt.bar(x = amino_acids, height = frequencies, color = colors, width = 0.8)  #Plot the bar graph

    #Set the graph title, labels
    plt.title('Frequency Distribution of Encoded Amino Acids')
    plt.xlabel('Amino Acid')
    plt.ylabel('Frequency')
    plt.xticks(rotation=45)                 #Rotate xticks for a more comfortable visualization
    plt.tight_layout()
    plt.show()



file = open('human_RSCU.tsv', 'r')
skip_row = 0
for line in file:
    if 'CODON' in line:
        break
    skip_row += 1
RSCU = pd.read_table('human_RSCU.tsv', skiprows = skip_row)
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
print(largest_RSCU_dict(RSCU))
print(CAI_calculator(seq))