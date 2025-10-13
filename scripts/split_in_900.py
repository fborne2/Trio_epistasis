#Split consensus sequences in chunks of 900 codons
import os
import math
import re
import shutil

input_folder = "consensus_output"
output_folder = "consensus_split"
max_len = 900

# Create output folder
os.makedirs(output_folder, exist_ok=True)

def read_sequence(filepath):
    with open(filepath, 'r') as f:
        return ''.join(line.strip() for line in f)

def write_parts(seq, basepath, prefix):
    """Split sequence into parts and write FASTA files with headers"""
    total_len = len(seq)
    num_parts = math.ceil(total_len / max_len)
    for i in range(num_parts):
        start = i * max_len
        end = min((i + 1) * max_len, total_len)
        part_seq = seq[start:end]
        outname = f"{basepath}_{i+1}_{prefix}_consensus.txt"
        with open(outname, 'w') as out:
            out.write(f">{prefix}\n")
            out.write(part_seq + "\n")
    return num_parts

# Collect transcript IDs from files
files = os.listdir(input_folder)
transcripts = set()

for f in files:
    m = re.match(r"(.+?)_(divergence|polymorphism)_consensus\.txt$", f)
    if m:
        transcripts.add(m.group(1))

# Process each transcript
for tid in transcripts:
    div_file = os.path.join(input_folder, f"{tid}_divergence_consensus.txt")
    poly_file = os.path.join(input_folder, f"{tid}_polymorphism_consensus.txt")

    if not (os.path.isfile(div_file) and os.path.isfile(poly_file)):
        continue

    div_seq = read_sequence(div_file)
    poly_seq = read_sequence(poly_file)

    # Sanity check
    if len(div_seq) != len(poly_seq):
        print(f"⚠️ Warning: sequence length mismatch for {tid} ({len(div_seq)} vs {len(poly_seq)})")
        continue

    base_out = os.path.join(output_folder, tid)

    # Case 1: sequences short enough → copy original files with headers
    if len(div_seq) <= max_len and len(poly_seq) <= max_len:
        for seq, prefix in [(div_seq, "divergence"), (poly_seq, "polymorphism")]:
            outname = f"{base_out}_{prefix}_consensus.txt"
            with open(outname, "w") as out:
                out.write(f">{prefix}\n")
                out.write(seq + "\n")
        print(f"✅ Copied {tid} without splitting")
        continue

    # Case 2: sequences too long → split into chunks
    num_parts = write_parts(div_seq, base_out, "divergence")
    write_parts(poly_seq, base_out, "polymorphism")
    print(f"✂️ Split {tid} into {num_parts} parts")

