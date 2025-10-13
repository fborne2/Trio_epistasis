#Filter gaps in alignments
import os
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

# Input and output directories
input_dir = "alignments"       # folder containing your OG*.fasta alignments
output_dir = "alignments_nogaps"

os.makedirs(output_dir, exist_ok=True)

# Loop over all fasta files
for filename in os.listdir(input_dir):
    if filename.endswith(".fasta"):
        infile = os.path.join(input_dir, filename)
        outfile = os.path.join(output_dir, filename)

        # Read alignment
        alignment = AlignIO.read(infile, "fasta")

        # Identify columns with no gaps
        keep_columns = []
        for i in range(alignment.get_alignment_length()):
            column = alignment[:, i]
            if "-" not in column:
                keep_columns.append(i)

        # Build filtered alignment
        filtered_records = []
        for record in alignment:
            new_seq = "".join(record.seq[i] for i in keep_columns)
            record.seq = record.seq.__class__(new_seq)
            filtered_records.append(record)

        filtered_alignment = MultipleSeqAlignment(filtered_records)

        # Write new alignment
        AlignIO.write(filtered_alignment, outfile, "fasta")

        print(f"Processed {filename}: {alignment.get_alignment_length()} → {filtered_alignment.get_alignment_length()} bp")
