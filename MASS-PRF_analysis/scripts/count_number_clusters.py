import ast

output_file = "all_OG_clusters_experimentally_tractable.txt"

total_clusters = 0

with open(output_file) as f:
    for line in f:
        parts = line.strip().split(" ", 1)
        if len(parts) < 2:
            continue
        clusters = ast.literal_eval(parts[1])  # safely convert string to list
        total_clusters += len(clusters)

print("Total number of clusters:", total_clusters)
