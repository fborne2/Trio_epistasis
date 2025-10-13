import os
import glob
import pandas as pd

input_pattern = "table_OG*_MASS-PRF_BIC*.txt"
output_file = "OG_clusters_experimentally_tractable.txt"

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

for fname in sorted(glob.glob(input_pattern)):
    df = pd.read_csv(fname, sep=r'\s+', engine='python', dtype=str)
    df.columns = df.columns.str.strip()
    
    # Convert numeric columns safely
    df["Position"] = pd.to_numeric(df["Position"], errors="coerce")
    df["Lower_CI_Gamma"] = pd.to_numeric(df["Lower_CI_Gamma"].str.replace("*","", regex=False), errors="coerce")
    
    # Clean DivergenceMutationStatus column
    status_col = "DivergenceMutationStatus"
    df[status_col] = df[status_col].astype(str).str.strip()
    
    # Filter positions where Lower_CI_Gamma >= 4 and status is R
    subset = df[(df["Lower_CI_Gamma"] >= 4) & (df[status_col] == "R")]
    positions = subset["Position"].dropna().astype(int).tolist()
    
    if not positions:
        continue
    
    clusters = find_clusters(positions, max_distance=20)
    
    # Keep only clusters of exactly 2 or 3 positions
    valid_clusters = [c for c in clusters if len(c) in (2,3)]
    
    if valid_clusters:
        # Clean OG ID: remove 'table_' prefix and '_MASS-PRF_BIC' suffix
        og_file = os.path.basename(fname).replace("table_", "").replace("_MASS-PRF_BIC","").replace(".txt","")
        results.append((og_file, valid_clusters))

# Save results
with open(output_file, "w") as out:
    for og_file, clusters in results:
        out.write(f"{og_file} {clusters}\n")

print(f"✅ Done! Found {len(results)} files with valid clusters. Results saved in {output_file}")

