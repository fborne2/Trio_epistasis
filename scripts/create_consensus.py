#Create consensus sequence from multiple alignments
import os
import re
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.Data import CodonTable
from collections import Counter

# Directory containing FASTA files
input_folder = "fasta_alignments"  # <-- Replace with your folder name
genetic_code = CodonTable.unambiguous_dna_by_id[1]

# Create output folder
output_folder = os.path.join(input_folder, "consensus")
os.makedirs(output_folder, exist_ok=True)

def process_fasta_file(file_path):
    try:
        alignment = AlignIO.read(file_path, "fasta")
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return

    ingroup_records = [rec for rec in alignment if rec.id.startswith("Dmel")]
    outgroup_records = [rec for rec in alignment if rec.id.startswith("Dsim")]

    if len(outgroup_records) != 1:
        print(f"Skipping {file_path}: expected 1 Dsim sequence, found {len(outgroup_records)}")
        return

    if not ingroup_records:
        print(f"Skipping {file_path}: no ingroup sequences found")
        return

    # Use filename as base name and remove suffix "_trimmed_with_anc_NT" if present
    fname = os.path.splitext(os.path.basename(file_path))[0]
    base_name = re.sub(r'_trimmed_with_anc_NT$', '', fname)

    outgroup_seq = outgroup_records[0].seq
    ingroup_seqs = [rec.seq for rec in ingroup_records]
    all_seqs = ingroup_seqs + [outgroup_seq]
    seq_length = len(outgroup_seq)

    # Replace 'N' with most common base at each position
    cleaned_seqs = []
    for seq in all_seqs:
        seq_list = list(seq)
        for i in range(seq_length):
            if seq_list[i] == "N":
                column = [s[i] for s in all_seqs if s[i] not in ("N", "-")]
                if column:
                    most_common = Counter(column).most_common(1)[0][0]
                    seq_list[i] = most_common
                else:
                    seq_list[i] = "N"
        cleaned_seqs.append(Seq("".join(seq_list)))

    ingroup_seqs = cleaned_seqs[:-1]
    outgroup_seq = cleaned_seqs[-1]

    poly_consensus = ""
    div_consensus = ""

    for i in range(0, seq_length, 3):
        out_codon = str(outgroup_seq[i:i+3])
        ingroup_codons = [str(seq[i:i+3]) for seq in ingroup_seqs]

        # Skip incomplete codons
        if len(out_codon) < 3 or any(len(c) < 3 for c in ingroup_codons):
            continue

        # Condition 1: all ingroup codons are NNN or all are ---
        if all(c == "NNN" for c in ingroup_codons) or all(c == "---" for c in ingroup_codons):
            poly_consensus += "-"
            div_consensus += "-"
            continue

        # Condition 2: outgroup codon is ---
        if out_codon == "---":
            poly_consensus += "-"
            div_consensus += "-"
            continue

        # Condition 3: any ingroup codon is ---
        if any(c == "---" for c in ingroup_codons):
            poly_consensus += "-"
            div_consensus += "-"
            continue

        # Condition 4: any codon contains 'N'
        if any('N' in c for c in ingroup_codons + [out_codon]):
            poly_consensus += "-"
            div_consensus += "-"
            continue

        ingroup_unique = set(ingroup_codons)

        # Polymorphism consensus
        if len(ingroup_unique) == 1:
            poly_consensus += "*"
        else:
            try:
                aa_set = set(Seq(c).translate(genetic_code) for c in ingroup_unique)
            except:
                aa_set = {"X"}
            poly_consensus += "S" if len(aa_set) == 1 else "R"

        # Divergence consensus
        if out_codon in ingroup_unique:
            div_consensus += "*"
        elif len(ingroup_unique) > 1:
            div_consensus += "-"
        else:
            a_codon = next(iter(ingroup_unique))
            try:
                a_aa = Seq(a_codon).translate(genetic_code)
                b_aa = Seq(out_codon).translate(genetic_code)
                div_consensus += "S" if a_aa == b_aa else "R"
            except:
                div_consensus += "-"

    # Write output files to consensus folder
    base_out = os.path.join(output_folder, base_name)
    with open(base_out + "_polymorphism_consensus.txt", "w") as f:
        f.write(poly_consensus + "\n")
    with open(base_out + "_divergence_consensus.txt", "w") as f:
        f.write(div_consensus + "\n")

    print(f"✓ Processed {os.path.basename(file_path)} → {base_name}_*.txt in 'consensus' folder")

# Process all .fasta files in the input folder
for filename in os.listdir(input_folder):
    if filename.endswith(".fasta"):
        process_fasta_file(os.path.join(input_folder, filename))

