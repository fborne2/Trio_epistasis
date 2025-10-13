import os
import pandas as pd
from collections import defaultdict, Counter
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

# === CONFIG ===
folder = "."
# priority list: first match wins
priority_order = [
    "Several ≥ 4 (some <20 apart)",
    "Several ≥ 4 (all ≥20 apart)",
    "One ≥ 4",
    "All < 4",
    "Error"   # treat classification exceptions as lowest priority
]
# OPTIONAL: move entire OGs of a chosen assigned category
move_entire_og = False
move_target_cat = "Several ≥ 4 (some <20 apart)"
# ================

# collect file list first (so moves don't affect the iteration)
files = [f for f in os.listdir(folder) if f.startswith("table_") and f.endswith(".txt")]

# maps
file_category = {}
og_files = defaultdict(list)  # OG -> list of filenames

# classify each file (exceptions map to "Error")
for fname in files:
    fpath = os.path.join(folder, fname)
    try:
        cat = classify_file(fpath)
    except Exception as e:
        cat = "Error"
    file_category[fname] = cat

    # Correct OG ID extraction: after 'table_'
    ogid = fname.split("_")[1]  # e.g., 'OG0001404'
    og_files[ogid].append(fname)

# assign one category per OG using the priority rule
og_assigned = {}
for ogid, fnames in og_files.items():
    cats_present = {file_category[f] for f in fnames}
    assigned = None
    for p in priority_order:
        if p in cats_present:
            assigned = p
            break
    if assigned is None:
        # fallback: any uncaught category (rare)
        assigned = sorted(cats_present)[0]
    og_assigned[ogid] = assigned

# counts
counts = Counter(og_assigned.values())

print("OG category breakdown (priority applied):")
for p in priority_order:
    print(f"{p}: {counts.get(p, 0)}")

# print any other (unexpected) categories
others = [c for c in counts.keys() if c not in priority_order]
for c in others:
    print(f"{c}: {counts[c]}")

print(f"\nTotal unique OGs: {len(og_assigned)}")

# OPTIONAL: move entire OGs if requested
if move_entire_og:
    target_folder = os.path.join(folder, move_target_cat.replace(" ", "_"))
    os.makedirs(target_folder, exist_ok=True)
    for ogid, assigned in og_assigned.items():
        if assigned == move_target_cat:
            for fname in og_files[ogid]:
                src = os.path.join(folder, fname)
                dst = os.path.join(target_folder, fname)
                if os.path.exists(src):
                    shutil.move(src, dst)
    print(f"\nMoved OGs assigned to '{move_target_cat}' into {target_folder}")

