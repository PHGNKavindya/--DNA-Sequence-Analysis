
import truststore
truststore.inject_into_ssl()

from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from pathlib import Path

# ==========================================
# 1. CONFIGURATION
# ==========================================

Entrez.email = "nethmikavindya828@gmail.com"

BRCA_GENES = {
    "BRCA1": "NG_005905.2",
    "BRCA2": "NG_012772.3"
}

MATCH = 1
MISMATCH = -1
GAP = -1

output_folder = Path("question3_output")
output_folder.mkdir(exist_ok=True)


# ==========================================
# 2. RETRIEVE ANNOTATED CDS
# ==========================================

def retrieve_cds(gene, accession):

    print(f"\nRetrieving {gene}: {accession}")

    handle = Entrez.efetch(
        db="nuccore",
        id=accession,
        rettype="gb",
        retmode="text"
    )

    record = SeqIO.read(handle, "genbank")
    handle.close()

    cds_features = [
        feature for feature in record.features
        if feature.type == "CDS"
    ]

    if not cds_features:
        raise ValueError(
            f"No annotated CDS found for {gene}"
        )

    # Select the longest annotated CDS feature
    cds_feature = max(
        cds_features,
        key=lambda feature: len(feature.location)
    )

    cds_sequence = str(
        cds_feature.extract(record.seq)
    ).upper()

    if any(base not in "ACGT" for base in cds_sequence):
        raise ValueError(
            f"Ambiguous nucleotide found in {gene} CDS"
        )

    print("CDS length:", len(cds_sequence), "bp")

    cds_record = SeqIO.SeqRecord(
        Seq(cds_sequence),
        id=gene,
        description=f"Longest annotated CDS from {accession}"
    )

    SeqIO.write(
        cds_record,
        output_folder / f"{gene}_CDS.fasta",
        "fasta"
    )

    return cds_sequence


# ==========================================
# 3. NEEDLEMAN-WUNSCH SCORE ROW
#    Uses linear memory
# ==========================================

def nw_score_row(seq1, seq2):

    previous = [
        j * GAP for j in range(len(seq2) + 1)
    ]

    for i in range(1, len(seq1) + 1):

        current = [i * GAP]

        for j in range(1, len(seq2) + 1):

            if seq1[i - 1] == seq2[j - 1]:
                substitution = MATCH
            else:
                substitution = MISMATCH

            diagonal = (
                previous[j - 1] + substitution
            )

            deletion = previous[j] + GAP
            insertion = current[j - 1] + GAP

            current.append(
                max(diagonal, deletion, insertion)
            )

        previous = current

    return previous


# ==========================================
# 4. NEEDLEMAN-WUNSCH BASE ALIGNMENT
#    For small sequences only
# ==========================================

def small_alignment(seq1, seq2):

    rows = len(seq1) + 1
    cols = len(seq2) + 1

    dp = [
        [0] * cols
        for _ in range(rows)
    ]

    for i in range(1, rows):
        dp[i][0] = i * GAP

    for j in range(1, cols):
        dp[0][j] = j * GAP

    for i in range(1, rows):

        for j in range(1, cols):

            substitution = (
                MATCH
                if seq1[i - 1] == seq2[j - 1]
                else MISMATCH
            )

            diagonal = (
                dp[i - 1][j - 1] + substitution
            )

            deletion = dp[i - 1][j] + GAP
            insertion = dp[i][j - 1] + GAP

            dp[i][j] = max(
                diagonal,
                deletion,
                insertion
            )

    aligned1 = []
    aligned2 = []

    i = len(seq1)
    j = len(seq2)

    while i > 0 or j > 0:

        if i > 0 and j > 0:

            substitution = (
                MATCH
                if seq1[i - 1] == seq2[j - 1]
                else MISMATCH
            )

            if dp[i][j] == (
                dp[i - 1][j - 1] + substitution
            ):
                aligned1.append(seq1[i - 1])
                aligned2.append(seq2[j - 1])
                i -= 1
                j -= 1
                continue

        if i > 0 and dp[i][j] == dp[i - 1][j] + GAP:

            aligned1.append(seq1[i - 1])
            aligned2.append("-")
            i -= 1

        else:

            aligned1.append("-")
            aligned2.append(seq2[j - 1])
            j -= 1

    return (
        "".join(reversed(aligned1)),
        "".join(reversed(aligned2))
    )


# ==========================================
# 5. HIRSCHBERG DIVIDE-AND-CONQUER
# ==========================================

def hirschberg(seq1, seq2):

    if len(seq1) == 0:
        return "-" * len(seq2), seq2

    if len(seq2) == 0:
        return seq1, "-" * len(seq1)

    if len(seq1) == 1 or len(seq2) == 1:
        return small_alignment(seq1, seq2)

    middle = len(seq1) // 2

    left_scores = nw_score_row(
        seq1[:middle],
        seq2
    )

    right_scores = nw_score_row(
        seq1[middle:][::-1],
        seq2[::-1]
    )

    split = max(
        range(len(seq2) + 1),
        key=lambda j: (
            left_scores[j] + right_scores[len(seq2) - j]
        )
    )

    left_alignment = hirschberg(
        seq1[:middle],
        seq2[:split]
    )

    right_alignment = hirschberg(
        seq1[middle:],
        seq2[split:]
    )

    return (
        left_alignment[0] + right_alignment[0],
        left_alignment[1] + right_alignment[1]
    )


# ==========================================
# 6. CALCULATE ALIGNMENT SCORE
# ==========================================

def calculate_alignment_score(aligned1, aligned2):

    score = 0

    matches = 0
    mismatches = 0
    gaps = 0

    for base1, base2 in zip(aligned1, aligned2):

        if base1 == "-" or base2 == "-":

            score += GAP
            gaps += 1

        elif base1 == base2:

            score += MATCH
            matches += 1

        else:

            score += MISMATCH
            mismatches += 1

    return score, matches, mismatches, gaps


# ==========================================
# 7. RETRIEVE BRCA CDS SEQUENCES
# ==========================================

sequences = {}

for gene, accession in BRCA_GENES.items():

    sequences[gene] = retrieve_cds(
        gene,
        accession
    )


# ==========================================
# 8. PERFORM HIRSCHBERG ALIGNMENT
# ==========================================

print("\nStarting Hirschberg alignment...")

aligned_brca1, aligned_brca2 = hirschberg(
    sequences["BRCA1"],
    sequences["BRCA2"]
)


# ==========================================
# 9. CALCULATE RESULTS
# ==========================================

alignment_score, matches, mismatches, gaps = (
    calculate_alignment_score(
        aligned_brca1,
        aligned_brca2
    )
)

alignment_length = len(aligned_brca1)

similarity_percentage = (
    matches / alignment_length * 100
)

print("\n========== ALIGNMENT RESULTS ==========")

print("BRCA1 CDS length:", len(sequences["BRCA1"]))
print("BRCA2 CDS length:", len(sequences["BRCA2"]))

print("Alignment length:", alignment_length)
print("Matching pairs:", matches)
print("Mismatching pairs:", mismatches)
print("Gap positions:", gaps)

print("Alignment score:", alignment_score)

print(
    "Identity percentage: {:.2f}%".format(
        similarity_percentage
    )
)


# ==========================================
# 10. SAVE ALIGNMENT
# ==========================================

alignment_path = output_folder / "BRCA_alignment.txt"

with open(alignment_path, "w") as file:

    file.write("BRCA1 - BRCA2 HIRSCHBERG ALIGNMENT\n")
    file.write("=================================\n\n")

    file.write(f"BRCA1 accession: {BRCA_GENES['BRCA1']}\n")
    file.write(f"BRCA2 accession: {BRCA_GENES['BRCA2']}\n")

    file.write(f"Match score: {MATCH}\n")
    file.write(f"Mismatch score: {MISMATCH}\n")
    file.write(f"Gap penalty: {GAP}\n\n")

    file.write("BRCA1 CDS:\n")
    file.write(sequences["BRCA1"] + "\n\n")

    file.write("BRCA2 CDS:\n")
    file.write(sequences["BRCA2"] + "\n\n")

    file.write("ALIGNED BRCA1:\n")
    file.write(aligned_brca1 + "\n\n")

    file.write("ALIGNED BRCA2:\n")
    file.write(aligned_brca2 + "\n\n")

    file.write(f"Alignment length: {alignment_length}\n")
    file.write(f"Matching pairs: {matches}\n")
    file.write(f"Mismatching pairs: {mismatches}\n")
    file.write(f"Gap positions: {gaps}\n")
    file.write(f"Alignment score: {alignment_score}\n")
    file.write(
        f"Identity percentage: {similarity_percentage:.2f}%\n"
    )

print("\nAlignment saved:", alignment_path)
print("\nHIRSCHBERG ALIGNMENT COMPLETED!")
