import os
import pandas as pd
from collections import Counter
import shutil

def classify_file(filepath):
    df = pd.read_csv(filepath, sep="\t")

    # Force Lower_CI_Gamma to numeric (non-numeric like '*' → NaN)
    df["Lower_CI_Gamma"] = pd.to_numeric(df["Lower_CI_Gamma"], errors="coerce")

    positions = df["Position"].tolist()
    gamma_vals = df["Lower_CI_Gamma"].tolist()

    # Keep only valid numeric gamma values
    sites = [pos for pos, g in zip(positions, gamma_vals) if pd.notna(g) and g >= 4]

    if len(sites) == 0:
        return "All < 4"
    elif len(sites) == 1:
        return "One ≥ 4"
    else:
        sites_sorted = sorted(sites)
        close_pairs = any(
            (sites_sorted[i+1] - sites_sorted[i]) < 20
            for i in range(len(sites_sorted)-1)
        )
        if close_pairs:
            return "Several ≥ 4 (some <20 apart)"
        else:
            return "Several ≥ 4 (all ≥20 apart)"

# ---- Run on one folder ----
folder = "tables_only"
results = {}

# Folder to move files into
target_folder = os.path.join(folder, "several_ge4_some_lt20")
os.makedirs(target_folder, exist_ok=True)

for fname in os.listdir(folder):
    if fname.endswith(".txt"):
        fpath = os.path.join(folder, fname)
        try:
            category = classify_file(fpath)
            results[fname] = category

            # Move files of the specific category
            if category == "Several ≥ 4 (some <20 apart)":
                shutil.move(fpath, os.path.join(target_folder, fname))

        except Exception as e:
            results[fname] = f"Error: {e}"

# ---- Counts for all categories ----
category_counts = Counter(results.values())

print("Category breakdown:")
for cat, n in category_counts.items():
    print(f"{cat}: {n}")

print(f"\nTotal files: {sum(category_counts.values())}")
print(f"\nMoved files of category 'Several ≥ 4 (some <20 apart)' to {target_folder}")

