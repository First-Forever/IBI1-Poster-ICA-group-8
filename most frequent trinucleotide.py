import math

def check_seq(seq):
    std_nucleotide = 'AGCU'         #Valid nucleotides
    for i in seq:
        if i not in std_nucleotide:         #Check if all nucleotides in the sequence is valid; if not, raise ValueError
            raise ValueError(f'The sequence has invalid nucleotide! Try again!\nInvalid nucleotide: {i}, location: {seq.find(i)}')

def most_frequent_trinucleotide(seq):
    stop_codon = ('UAG', 'UAA', 'UGA')                              #Define the sequence of stop codon
    try:                                                            #Check if the first three nucleotides is a start codon
        if seq[0:3] != 'AUG':
            raise ValueError("Start codon not valid! Try again!")   
    except ValueError as e:
        print(e)
    seq_len = len(seq)
    freq_dict = {}                                                  #Use dictionary to calculate all the frequency of codons (to save time)
    for i in range(math.ceil(seq_len/3)):                           
        now_seq = seq[i*3 : i*3 + 3]                                #Slice every three nucleotides
        if now_seq in stop_codon:                                   #If stop codon is met, break the loop
            break
        freq_dict[now_seq] = freq_dict.get(now_seq, 0) + 1                   #If the codon is in the dictionary, add 1 to the value; if not, add it as a key and give the value as 1
    mmax_value = 0                                                  #Initialize the variables to find the most frequent sequence
    mmax_key = []
    for key in freq_dict:
        if freq_dict[key] > mmax_value:                             #If the value is bigger than mmax_value, clear the answer list and append 
            mmax_key.clear()
            mmax_key.append(key)
            mmax_value = freq_dict[key]
        elif freq_dict[key] == mmax_value:                          #If they are the same, append the key into answer list
            mmax_key.append(key)
    return mmax_key, mmax_value, freq_dict                                     #Return the most sequence codons and their frequency
seq = 'AUGGGAGGGGGGUUU'
check_seq(seq)
print(most_frequent_trinucleotide(seq)[0:2])