import os
import numpy as np
import gzip, shutil
from Bio.Data.IUPACData import protein_letters_3to1

def generate_fasta(file_name):
    seq = ""
    with open(file_name, "r") as fin:
        for line in fin:
            if line.startswith('SEQRES'):
                split_seq = line.split()
                split_seq = split_seq[4:]  # remove initial data until protein codes
                for item in split_seq:
                    seq += str(amino_acid_dic.get(item.capitalize()))

    with open("../output/fasta/" + file_name + '.fasta', "w") as fout:
        fout.write(seq)
    return seq


def generate_ss(file_name, seq):
    #generate amino acid sequence
    amino_acid_seq = ""
    seq_len = len(seq)
    amino_acid_seq = amino_acid_seq.ljust(seq_len, "C")
    with open(file_name, "r") as fin, open("../output/ss/" + file_name + '.ss', "w") as fout:
        index = 1
        for line in fin:
            if line.startswith('HELIX') or line.startswith('SHEET'):
                fout.write(f"{index}\n")
                if line.startswith('HELIX'):
                    start = int(line[21:25].strip())
                    end = int(line[33:37].strip())
                    if not (start > seq_len or end > seq_len or start < 0): #TODO: why some start indexs are - and some numbers are bigger
                        fout.write(f"{start}:{end}\n")
                        fout.write(f"{seq[start-1:end-1]}\n")
                        replacement_string = ""
                        replacement_string = replacement_string.ljust(end - start, "H")
                        amino_acid_seq = amino_acid_seq[:start-1] + replacement_string + amino_acid_seq[end-1:]
                else:
                    start = int(line[22:26].strip())
                    end = int(line[33:37].strip())
                    if not (start > seq_len or end > seq_len or start < 0):
                        replacement_string = ""
                        replacement_string = replacement_string.ljust(end - start, "B")
                        fout.write(f"{start}:{end}\n")
                        fout.write(f"{seq[start-1:end-1]}\n")
                        amino_acid_seq = amino_acid_seq[:start-1] + replacement_string + amino_acid_seq[end-1:]
                fout.write(f"{line[0:6]}\n")
                index += 1

    with open("../output/AA-SS.txt", "a") as fout:
        fout.write(f"{seq}\n")
        fout.write(f"{amino_acid_seq}\n\n\n")
    calculate_frequency(seq, amino_acid_seq)



def calculate_frequency(seq, amino_acid_seq):
    prev_char = ""
    for i in range(0, len(amino_acid_seq)):
        prev_char = prev_char + amino_acid_seq[i]
        if len(prev_char) >= 2:
            row = prev_char[-2]  # get second to last character
            col = prev_char[-1]  # get last character
            col = 0 if col == "H" else 1 if col == "B" else 2
            secondary_seq.get(row)[col] += 1

        if seq[i] in amino_structure: #check for key existance TODO: check protein names other than 20
            if amino_acid_seq[i] == "H":
                amino_structure.get(seq[i])[0] += 1
            elif amino_acid_seq[i] == "B":
                amino_structure.get(seq[i])[1] += 1
            elif amino_acid_seq[i] == "C":
                amino_structure.get(seq[i])[2] += 1


def calculate_structure_probability():
    for key in amino_structure:
        total_occurrence = np.sum(amino_structure[key])
        if total_occurrence != 0:
            for i, item in enumerate(amino_structure[key]):
                amino_structure[key][i] = round(item/total_occurrence,2)


def calculate_binary_seq_probability():
    for key in secondary_seq:
        total_occurrence = np.sum(secondary_seq[key])
        if total_occurrence != 0:
            for i, item in enumerate(secondary_seq[key]):
                secondary_seq[key][i] = round(item/total_occurrence, 2)


def generate_summary(file_count, file_names, file_lengths):
    avg = np.average(file_lengths)
    aaa2standard_deviation = np.std(file_lengths)
    with open("../output/summary.txt", "w") as fout:
        fout.write(f"- Total number of processed PDBs: {file_count} \n")
        fout.write(f"- Average of processed protein lengths: {avg} \n")
        fout.write(f"- Standard deviation of processed protein lengths: {standard_deviation} \n")
        fout.write(f"- Names of processed proteins: \n")
        fout.write(f"{file_names} \n\n\n")

        # write the probability of the presence of each of the twenty amino acids in each of the three states of the
        # type II structure
        fout.write(f"    Helix Sheet coil\n")
        structure_matrix = str(amino_structure)
        structure_matrix = structure_matrix.replace(", '", "\n")
        fout.write(f"{structure_matrix[2:]}\n\n\n")

        # write the probability of binary sequences of second type structure elements
        fout.write(f"     H     B      C\n")
        secondary_seq_matrix = str(secondary_seq)
        secondary_seq_matrix = secondary_seq_matrix.replace(", '", "\n")
        fout.write(secondary_seq_matrix[2:])
    print("summary file created")



amino_acid_dic = {
        "Ala": "A",
        "Arg": "R",
        "Asn": "N",
        "Asp": "D",
        "Cys": "C",
        "Glu": "E",
        "Gln": "Q",
        "Gly": "G",
        "His": "H",
        "Ile": "I",
        "Leu": "L",
        "Lys": "K",
        "Met": "M",
        "Phe": "F",
        "Pro": "P",
        "Ser": "S",
        "Thr": "T",
        "Trp": "W",
        "Tyr": "Y",
        "Val": "V",
        "Ace": "X",
        "Dm0": "K",
}

amino_structure = {
    # The first number shows that amino frequency in helix type, the second and third numbers are beta and coil, respectively
    "A": [0, 0, 0],
    "R": [0, 0, 0],
    "N": [0, 0, 0],
    "D": [0, 0, 0],
    "C": [0, 0, 0],
    "Q": [0, 0, 0],
    "E": [0, 0, 0],
    "G": [0, 0, 0],
    "H": [0, 0, 0],
    "I": [0, 0, 0],
    "L": [0, 0, 0],
    "K": [0, 0, 0],
    "M": [0, 0, 0],
    "F": [0, 0, 0],
    "P": [0, 0, 0],
    "S": [0, 0, 0],
    "T": [0, 0, 0],
    "W": [0, 0, 0],
    "Y": [0, 0, 0],
    "V": [0, 0, 0],
    "X": [0, 0, 0],
}

secondary_seq = {
    # in array The first number indicates helix, the second and third numbers are beta and coil, respectively
    "H": [0, 0, 0],
    "B": [0, 0, 0],
    "C": [0, 0, 0],
}

#make output directories
if not os.path.exists("output"):
    os.makedirs("output")
    if not os.path.exists("output/fasta"):
        os.makedirs("output/fasta")
    if not os.path.exists("output/ss"):
        os.makedirs("output/ss")


file_names = []
file_lengths = []
def gz_extract(directory):
    extension = ".gz"
    os.chdir(directory)
    file_count = 0
    for item in os.listdir("."):  # loop through items in dataset dir
        file_count += 1
        if item.endswith(extension):  # check for ".gz" extension
            #decompress the pdf file
            gz_name = os.path.abspath(item)  # get full path of files
            file_name = (os.path.basename(gz_name)).rsplit('.', 1)[0]  # get file name for file within
            with gzip.open(gz_name, "rb") as f_in, open(file_name, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
            os.remove(gz_name)  # delete zipped file

            # process the pdf file
            sequence = generate_fasta(file_name)
            file_lengths.append(len(sequence))
            generate_ss(file_name, sequence)
            file_name = file_name.split('.', 1)[0]
            file_names.append(file_name[3:])
    return file_count


processed_file_count = gz_extract('dataset')
calculate_structure_probability()
calculate_binary_seq_probability()
generate_summary(processed_file_count, file_names, file_lengths)