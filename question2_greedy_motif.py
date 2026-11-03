
import truststore
truststore.inject_into_ssl()

from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from pathlib import Path
from collections import Counter

# ==========================================
# 1. CONFIGURATION
# ==========================================

Entrez.email = "nethmikavindya828@gmail.com"

GENES = {
    "G6PD": "NG_009015.2",
    "F8": "NG_011403.2",
    "DMD": "NG_012232.1",
    "PAH": "NG_008690.2"
}

K = 10
WINDOW_SIZE = 1000
PSEUDOCOUNT = 1

output_folder = Path("question2_genomic_output")
output_folder.mkdir(exist_ok=True)

# ==========================================
# 2. RETRIEVE GENOMIC RECORDS
#    AND EXTRACT GENE REGIONS
# ==========================================

def retrieve_gene_region(gene, accession):

    print(f"\nRetrieving {gene}: {accession}")

    handle = Entrez.efetch(
        db="nuccore",
        id=accession,
        rettype="gb",
        retmode="text"
    )

    record = SeqIO.read(handle, "genbank")
    handle.close()

    # Save the complete genomic RefSeqGene record
    SeqIO.write(
        record,
        output_folder / f"{gene}_genomic_record.fasta",
        "fasta"
    )

    # Find the annotated gene feature
    gene_feature = None

    for feature in record.features:

        if feature.type == "gene":

            gene_names = feature.qualifiers.get("gene", [])

            if gene in gene_names:
                gene_feature = feature
                break

    if gene_feature is None:
        raise ValueError(
            f"Gene feature {gene} not found in {accession}"
        )

    # Extract gene sequence in biological orientation.
    # Biopython handles reverse-strand genes automatically.
    gene_sequence = str(
        gene_feature.extract(record.seq)
    ).upper()

    if len(gene_sequence) < WINDOW_SIZE:
        raise ValueError(
            f"{gene} has fewer than {WINDOW_SIZE} bases"
        )

    # Select the first 1000 bases from the gene's 5' end
    selected_sequence = gene_sequence[:WINDOW_SIZE]

    if any(base not in "ACGT" for base in selected_sequence):
        raise ValueError(
            f"{gene} selected region contains ambiguous bases. "
            "Inspect the NCBI record before continuing."
        )

    print("Complete record length:", len(record.seq))
    print("Annotated gene length:", len(gene_sequence))
    print("Selected region length:", len(selected_sequence))
    print("Gene strand:", gene_feature.location.strand)

    # Save the selected gene region
    selected_record = SeqIO.SeqRecord(
        Seq(selected_sequence),
        id=gene,
        description=(
            f"{accession}; first {WINDOW_SIZE} bases "
            "of gene in 5-to-3 orientation"
        )
    )

    SeqIO.write(
        selected_record,
        output_folder / f"{gene}_selected_region.fasta",
        "fasta"
    )

    return selected_sequence


# ==========================================
# 3. RETRIEVE ALL FOUR GENES
# ==========================================

sequences = {}

for gene, accession in GENES.items():
    sequences[gene] = retrieve_gene_region(
        gene,
        accession
    )


# ==========================================
# 4. CREATE PROFILE WITH PSEUDOCOUNTS
# ==========================================

def create_profile(motifs):

    profile = []

    for position in range(K):

        counts = {
            nucleotide: PSEUDOCOUNT
            for nucleotide in "ACGT"
        }

        for motif in motifs:
            counts[motif[position]] += 1

        total = sum(counts.values())

        probabilities = {
            nucleotide: counts[nucleotide] / total
            for nucleotide in "ACGT"
        }

        profile.append(probabilities)

    return profile


# ==========================================
# 5. FIND MOST PROBABLE K-MER
# ==========================================

def find_best_kmer(sequence, profile):

    best_kmer = None
    best_position = None
    best_probability = -1.0

    for i in range(len(sequence) - K + 1):

        kmer = sequence[i:i + K]

        if any(base not in "ACGT" for base in kmer):
            continue

        probability = 1.0

        for position, nucleotide in enumerate(kmer):
            probability *= profile[position][nucleotide]

        if probability > best_probability:
            best_probability = probability
            best_kmer = kmer
            best_position = i

    if best_kmer is None:
        raise ValueError("No valid DNA k-mer found.")

    return best_kmer, best_position


# ==========================================
# 6. CALCULATE MOTIF SCORE
# ==========================================

def calculate_score(motifs):

    score = 0

    for position in range(K):

        column = [
            motif[position]
            for motif in motifs
        ]

        most_common_count = max(
            Counter(column).values()
        )

        score += len(motifs) - most_common_count

    return score


# ==========================================
# 7. GREEDY MOTIF SEARCH
# ==========================================

def greedy_motif_search(dna_sequences):

    best_motifs = None
    best_positions = None
    best_score = float("inf")

    first_sequence = dna_sequences[0]

    for i in range(len(first_sequence) - K + 1):

        first_kmer = first_sequence[i:i + K]

        if any(base not in "ACGT" for base in first_kmer):
            continue

        motifs = [first_kmer]
        positions = [i]

        for j in range(1, len(dna_sequences)):

            profile = create_profile(motifs)

            next_motif, next_position = find_best_kmer(
                dna_sequences[j],
                profile
            )

            motifs.append(next_motif)
            positions.append(next_position)

        score = calculate_score(motifs)

        if score < best_score:

            best_score = score
            best_motifs = motifs.copy()
            best_positions = positions.copy()

    if best_motifs is None:
        raise ValueError("Greedy search could not find valid motifs.")

    return best_motifs, best_positions, best_score


# ==========================================
# 8. RUN THE ALGORITHM
# ==========================================

print("\n======================================")
print("STARTING GREEDY MOTIF SEARCH")
print("======================================")

gene_names = list(GENES.keys())

dna_sequences = [
    sequences[gene]
    for gene in gene_names
]

best_motifs, positions, best_score = greedy_motif_search(
    dna_sequences
)


# ==========================================
# 9. DISPLAY RESULTS
# ==========================================

print("\n========== FINAL MOTIF RESULTS ==========")

for gene, motif, position in zip(
    gene_names,
    best_motifs,
    positions
):

    print(
        f"{gene}: {motif} "
        f"(position {position + 1})"
    )

print("\nMotif length:", K)
print("Pseudocount:", PSEUDOCOUNT)
print("Best motif score:", best_score)


# ==========================================
# 10. SAVE RESULTS
# ==========================================

results_path = output_folder / "motif_results.txt"

with open(results_path, "w") as file:

    file.write("GENOMIC DNA GREEDY MOTIF SEARCH\n")
    file.write("================================\n\n")

    file.write("NCBI genomic accessions:\n")

    for gene, accession in GENES.items():
        file.write(f"{gene}: {accession}\n")

    file.write(f"\nSelected region per gene: {WINDOW_SIZE} bp\n")
    file.write("Region: first 1000 bases in gene 5'-to-3' orientation\n")
    file.write(f"Motif length: {K}\n")
    file.write(f"Pseudocount: {PSEUDOCOUNT}\n\n")

    for gene, motif, position in zip(
        gene_names,
        best_motifs,
        positions
    ):

        file.write(
            f"{gene}: {motif} "
            f"(position {position + 1})\n"
        )

    file.write(f"\nBest motif score: {best_score}\n")

print("\nResults saved:", results_path)
print("\nGENOMIC GREEDY MOTIF SEARCH COMPLETED!")
