##Create list of sequences that contains at least one cluster of adaptive substitutions
import os
import glob
import pandas as pd

# Input file pattern
input_pattern = "table_OG*_MASS-PRF_BIC*.txt"
output_file = "OG_clusters_filtered.txt"

def find_clusters(positions, max_distance=20):
    """
    Find clusters of consecutive positions where each position is within max_distance
    from the previous one.
    Returns all clusters.
    """
    positions = sorted(positions)
    clusters = []
    cluster = []

    for pos in positions:
        if not cluster:
            cluster = [pos]
        elif pos - cluster[-1] <= max_distance:
            cluster.append(pos)
        else:
            clusters.append(cluster)
            cluster = [pos]
    if cluster:
        clusters.append(cluster)

    return clusters

results = []

for fname in glob.glob(input_pattern):
    og = os.path.basename(fname).split("_")[1]  # OG ID from filename, e.g., OG0001408

    # Load table
    df = pd.read_csv(fname, sep=r'\s+', engine='python', dtype=str)
    df.columns = df.columns.str.strip()

    # Convert numeric columns
    df["Position"] = pd.to_numeric(df["Position"], errors="coerce")
    df["Lower_CI_Gamma"] = pd.to_numeric(df["Lower_CI_Gamma"].str.replace("*","", regex=False), errors="coerce")

    # Clean status column
    status_col = "DivergenceMutationStatus"
    df[status_col] = df[status_col].astype(str).str.strip()

    # Filter positions of interest
    subset = df[(df["Lower_CI_Gamma"] >= 4) & (df[status_col] == "R")]

    if subset.empty:
        continue

    positions = subset["Position"].dropna().astype(int).tolist()
    clusters = find_clusters(positions, max_distance=20)

    # Keep only clusters of exactly 2 or 3
    valid_clusters = [c for c in clusters if len(c) in (2,3)]

    if valid_clusters:
        results.append((og, valid_clusters))

# Save results
with open(output_file, "w") as out:
    for og, clusters in results:
        out.write(f"{og}\t{clusters}\n")

print(f"✅ Done! Found {len(results)} OGs with valid clusters. Results saved in {output_file}")

