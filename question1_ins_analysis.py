import truststore
truststore.inject_into_ssl()


from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from pathlib import Path

# --------------------------------
# 1. NCBI configuration
# --------------------------------

Entrez.email = "nethmikavindya828@gmail.com"

ACCESSION = "NC_060935.1"
START = 2247427
END = 2248857

# Create output folder
output_folder = Path("question1_output")
output_folder.mkdir(exist_ok=True)

# --------------------------------
# 2. Retrieve INS DNA sequence
# --------------------------------

print("Retrieving INS gene sequence from NCBI...")

handle = Entrez.efetch(
    db="nuccore",
    id=ACCESSION,
    seq_start=START,
    seq_stop=END,
    strand=1,
    rettype="fasta",
    retmode="text"
)

record = SeqIO.read(handle, "fasta")
handle.close()

dna_sequence = str(record.seq).upper()

print("\nSequence retrieved successfully!")

# --------------------------------
# 3. Display sequence information
# --------------------------------

print("\n--- INS Gene Information ---")

print("Accession Number:", ACCESSION)
print("Selected Region:", START, "-", END)
print("Sequence Length:", len(dna_sequence), "bp")

# Check expected sequence length
if len(dna_sequence) != (END - START + 1):
    raise ValueError("Sequence length does not match expected region!")

# --------------------------------
# 4. Save FASTA sequence
# --------------------------------

fasta_path = output_folder / "INS_sequence.fasta"

SeqIO.write(record, fasta_path, "fasta")

print("\nFASTA file saved:", fasta_path)

# --------------------------------
# 5. Calculate GC content
# --------------------------------

a_count = dna_sequence.count("A")
t_count = dna_sequence.count("T")
g_count = dna_sequence.count("G")
c_count = dna_sequence.count("C")

valid_bases = a_count + t_count + g_count + c_count

gc_content = ((g_count + c_count) / valid_bases) * 100

print("\n--- GC Content Analysis ---")

print("Adenine (A):", a_count)
print("Thymine (T):", t_count)
print("Guanine (G):", g_count)
print("Cytosine (C):", c_count)

print("GC Content: {:.2f}%".format(gc_content))

# --------------------------------
# 6. Generate reverse complement
# --------------------------------

sequence_object = Seq(dna_sequence)

reverse_complement = str(
    sequence_object.reverse_complement()
)

print("\n--- Reverse Complement ---")

print("Reverse complement generated successfully!")

print("First 100 bases:")
print(reverse_complement[:100])

# Save reverse complement as FASTA
reverse_record = record[:]
reverse_record.id = "INS_reverse_complement"
reverse_record.description = "INS reverse complement sequence"

reverse_record.seq = Seq(reverse_complement)

reverse_path = output_folder / "INS_reverse_complement.fasta"

SeqIO.write(reverse_record, reverse_path, "fasta")

print("\nReverse complement saved:", reverse_path)

# --------------------------------
# 7. Save analysis results
# --------------------------------

results_path = output_folder / "INS_analysis_results.txt"

with open(results_path, "w") as file:
    file.write("INS GENE SEQUENCE ANALYSIS\n")
    file.write("==========================\n")
    file.write(f"Accession: {ACCESSION}\n")
    file.write(f"Region: {START}-{END}\n")
    file.write(f"Sequence Length: {len(dna_sequence)} bp\n")
    file.write(f"A Count: {a_count}\n")
    file.write(f"T Count: {t_count}\n")
    file.write(f"G Count: {g_count}\n")
    file.write(f"C Count: {c_count}\n")
    file.write(f"GC Content: {gc_content:.2f}%\n")
    file.write("\nReverse Complement:\n")
    file.write(reverse_complement)

print("\nAnalysis results saved:", results_path)

print("\n================================")
print("INS ANALYSIS COMPLETED SUCCESSFULLY!")
print("================================")
