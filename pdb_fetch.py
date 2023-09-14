# Purpose: Fetch a PDB file from the RCSB PDB database using its ID
import requests
from Bio.Data import CodonTable
# User defines the input PDB ID
pdb_id = input("Enter a PDB ID: ")
try:
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    print(f"Fetching PDB file from {url}")
    response = requests.get(url)
    print(f"HTTP response code: {response.status_code}")
    #code 200 means the request was successful
    if response.status_code == 200:
        print("PDB file fetched successfully")
        # Save the PDB file to a variable
        pdb_data = response.text
        # extracting the sequence from the PDB file and filter out non-ATCG characters starting from breaking sequence into single lines
        lines = pdb_data.split("\n")
        sequence = ""
        for line in lines:
            if line.startswith("SEQRES"):    # Sequence is stored in lines starting with "SEQRES"
                parts = line.split() # slitting string of characters into list
                sequence += "".join(parts[4:]) # because sequence starts as 5th element in the previously splitted list: it starts like this: 'SEQRES', '1', 'A', '5', 'A', 'T', 'C' etc., so the 5th element is the first base. 
                # using an empty string ("") as the delimiter -  nothing gets put between the list elements when joining them.
        filtered_sequence = "".join([base for base in sequence if base in "ATCG"])
        # if there was a problem with fetching the PDB file:
    else:
        print("Couldn't fetch the PDB file. Are you sure the ID is correct?")
        # if there was a problem with the HTTP request:
except Exception as error:
        print(f"Something went wrong: {error}")
# AT content
def at_content(sequence):
    # return the percentage of A and T bases in the sequence
    return(sequence.count("A") + sequence.count("T")) / len(filtered_sequence)*100
print(f"AT content: {at_content(filtered_sequence)}%")
# GC content
def gc_content(sequence):
    # return the percentage of G and C bases in the sequence
    return (sequence.count("G") + sequence.count("C")) / len(filtered_sequence)*100
print(f"GC content: {gc_content(filtered_sequence)}%")

# Transcription
def transcription(sequence):
    # replace each T with U
    return sequence.replace("T", "U")
print(f"Transcription: {transcription(filtered_sequence)}")

# Reverse Complement 
def reverse_complement(sequence):
    # reverse the sequence and replace each base with its complement
    return sequence[::-1].translate(str.maketrans("ATCG", "TAGC"))
print(f"Reverse Complement: {reverse_complement(filtered_sequence)}")

# Translation
def translation(sequence):
    # getting first unambiguous DNA codon table and storing it in the variable table 
    table = CodonTable.unambiguous_dna_by_id[1]
    # creating an empty string for the protein sequence
    protein = ""
    if len(sequence) % 3 != 0: # if the sequence length is not divisible by 3, remove the extra bases
     return "The sequence length is not divisible by 3"
 # if it is divisible by 3, translate the sequence by iterating over it in steps of 3
    for i in range(0, len(sequence), 3):
        # get the codon and translate it to an amino acid i:i+3 means that we are taking 3 bases at a time
        codon = sequence[i:i+3]
        # if the codon is not in the table, return "?"
        amino_acid = table.forward_table.get(codon, '?')
        # asterisk means stop codon, so if we encounter it, stop the translation
        if amino_acid == '*':
            break
        # add the amino acid to the protein sequence string
        protein += amino_acid
    return protein
print(f"Translation: {translation(filtered_sequence)}")
    
    
        