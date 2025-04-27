# 1. Define a new function called most_frequent_AA, and make it take a single argument, mRNA_sequence.
# 2. Inside the function, create an empty list called codons.
# 3. Use a for loop to iterate over the mRNA_sequence in steps of 3, extracting each codon (3 nucleotides) and appending it to the codons list.
# 4. Use the max function to find the most common codon in the codons list. 
# 5. Create a dictionary called genetic_code that maps each codon to its corresponding amino acid or "Stop" for stop codons.
# 6. Use the genetic_code dictionary to find the corresponding amino acid for the most common codon.
# 7. Return the most frequent amino acid from the function.
# 8. Call the function with an example mRNA sequence and print the result.

def most_frequent_AA (mRNA_sequence): # Function name: most_frequent_AA, input: mRNA_sequence
    codons = [] # Empty list called codons.
    for i in range(0, len(mRNA_sequence), 3): # For loop to iterate over the mRNA_sequence in steps of 3.
        codon = mRNA_sequence[i:i+3] # Extract each codon.
        if len(codon) == 3: # Check if the codon is of length 3.
            codons.append(codon) # Append the codon to the codons list.
    most_common_codon = max(codons, key=codons.count) # Find the most common codon in the codons list.
    genetic_code = { # Dictionary called genetic_code.
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
    most_frequent_amino_acid = genetic_code.get(most_common_codon, 'Unknown') # Find the corresponding amino acid.
    return most_frequent_amino_acid # Return the result.

print(most_frequent_AA(input("Please enter the mRNA sequence:"))) # Let the user to use this function.