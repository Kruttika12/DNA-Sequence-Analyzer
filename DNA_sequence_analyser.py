
from Bio.Seq import Seq
from collections import Counter
import csv
file_name = input("Enter the FASTA file path: ")
sequences = {}
with open(file_name, "r") as file:
    current_id = ""
    current_sequence = ""
    for line in file:
        line = line.strip()
        if line.startswith(">"):
            if current_id != "":
                sequences[current_id] = current_sequence
            current_id = line[1:].split()[0]
            current_sequence = ""
        else:
            current_sequence = current_sequence + line.upper()
    if current_id != "":
        sequences[current_id] = current_sequence
if len(sequences) == 0:
    print("No sequences found.")
    exit()
for sequence_id in sequences:
    sequence = sequences[sequence_id]
    if sequence == "":
        print("Empty sequence found:", sequence_id)
        exit()
    for base in sequence:
        if base not in "ATGCN":
            print("Invalid character found in", sequence_id)
            exit()
print("\nFASTA file loaded successfully.")
print("Number of sequences:", len(sequences))
print("\nAvailable sequences:")
for sequence_id in sequences:
    print(sequence_id)
selected_id = input("\nEnter the sequence ID: ")
if selected_id not in sequences:
    print("Invalid sequence ID.")
    exit()
sequence = sequences[selected_id]
print("\n1. SEQUENCE STATISTICS")
length = len(sequence)
counts = Counter(sequence)
A = counts["A"]
T = counts["T"]
G = counts["G"]
C = counts["C"]
N = counts["N"]
gc_content = ((G + C) / length) * 100
at_content = ((A + T) / length) * 100
print("Sequence length:", length)
print("A count:", A)
print("T count:", T)
print("G count:", G)
print("C count:", C)
print("N count:", N)
print("GC content:", round(gc_content, 2), "%")
print("AT content:", round(at_content, 2), "%")
print("\nNucleotide percentages:")
for base in "ATGCN":
    percentage = (counts[base] / length) * 100
    print(base, ":", round(percentage, 2), "%")
print("\n2. REVERSE COMPLEMENT")
complement = {
    "A": "T",
    "T": "A",
    "G": "C",
    "C": "G",
    "N": "N"
}
reverse_sequence = ""
for base in sequence[::-1]:
    reverse_sequence = reverse_sequence + complement[base]
print(reverse_sequence)
print("\n3. MOTIF DETECTION")
motif = input("Enter a DNA motif: ").upper()
if motif == "":
    print("Motif cannot be empty.")
elif not set(motif).issubset(set("ATGCN")):
    print("Invalid motif.")
else:
    positions = []
    for i in range(len(sequence) - len(motif) + 1):
        if sequence[i:i + len(motif)] == motif:
            positions.append(i + 1)
    print("Number of occurrences:", len(positions))
    print("Motif positions:", positions)
print("\n4. OPEN READING FRAME DETECTION")
stop_codons = ["TAA", "TAG", "TGA"]
reverse_complement = str(Seq(sequence).reverse_complement())
orfs = []
for strand, dna_sequence in [
    ("Forward", sequence),
    ("Reverse", reverse_complement)
]:
    for frame in range(3):
        for i in range(frame, len(dna_sequence) - 2, 3):
            codon = dna_sequence[i:i + 3]
            if codon == "ATG":
                for j in range(i + 3, len(dna_sequence) - 2, 3):
                    stop_codon = dna_sequence[j:j + 3]
                    if stop_codon in stop_codons:
                        orf_sequence = dna_sequence[i:j + 3]
                        start = i + 1
                        end = j + 3
                        if strand == "Reverse":
                            original_start = length - start + 1
                            original_end = length - end + 1
                        else:
                            original_start = start
                            original_end = end
                        orfs.append({
                            "Strand": strand,
                            "Frame": frame + 1,
                            "Start": original_start,
                            "End": original_end,
                            "Length": len(orf_sequence),
                            "Sequence": orf_sequence
                        })
                        break
print("Number of ORFs detected:", len(orfs))
for i, orf in enumerate(orfs):
    print("\nORF", i + 1)
    print("Strand:", orf["Strand"])
    print("Frame:", orf["Frame"])
    print("Start:", orf["Start"])
    print("End:", orf["End"])
    print("Length:", orf["Length"])
print("\n5. PROTEIN TRANSLATION")
for i, orf in enumerate(orfs):
    protein = str(Seq(orf["Sequence"]).translate(table=1))
    print("\nORF", i + 1)
    print("Protein sequence:")
    print(protein)
print("\n6. SEQUENCE COMPARISON")
comparison = []
for sequence_id in sequences:
    dna = sequences[sequence_id]
    base_counts = Counter(dna)
    gc = ((base_counts["G"] + base_counts["C"]) / len(dna)) * 100
    at = ((base_counts["A"] + base_counts["T"]) / len(dna)) * 100
    comparison.append({
        "Sequence_ID": sequence_id,
        "Length": len(dna),
        "A": base_counts["A"],
        "T": base_counts["T"],
        "G": base_counts["G"],
        "C": base_counts["C"],
        "GC_Content": round(gc, 2),
        "AT_Content": round(at, 2)

    })
for result in comparison:
    print(result)
with open("sequence_comparison.csv", "w", newline="") as file:
    columns = [
        "Sequence_ID",
        "Length",
        "A",
        "T",
        "G",
        "C",
        "GC_Content",
        "AT_Content"
    ]
    writer = csv.DictWriter(file, fieldnames=columns)
    writer.writeheader()
    writer.writerows(comparison)
with open("orf_results.csv", "w", newline="") as file:
    columns = [
        "Strand",
        "Frame",
        "Start",
        "End",
        "Length",
        "Sequence"
    ]
    writer = csv.DictWriter(file, fieldnames=columns)
    writer.writeheader()
    writer.writerows(orfs)
print("\nResults saved successfully.")
print("Files created:")
print("1. sequence_comparison.csv")
print("2. orf_results.csv")
