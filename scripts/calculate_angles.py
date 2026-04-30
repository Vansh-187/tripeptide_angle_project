import os
import sys
import gzip
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from tqdm import tqdm

# -----------------------------
# INPUTS
# -----------------------------
contexts_dir = sys.argv[1]   
output_file = sys.argv[2]   
target_aa = sys.argv[3]      

pdb_dir = "pdbs"  


size_map = {
    "K":"Large","R":"Bulky","D":"Intermediate","E":"Large",
    "G":"Tiny","A":"Tiny","V":"Small","L":"Intermediate","I":"Intermediate",
    "M":"Large","P":"Small","S":"Small","T":"Small","N":"Intermediate",
    "Q":"Large","C":"Small","F":"Bulky","Y":"Bulky","W":"Bulky","H":"Large"
}

def signed_angle_3d(v1, v2, axis):
    v1 = v1 / np.linalg.norm(v1)
    v2 = v2 / np.linalg.norm(v2)
    axis = axis / np.linalg.norm(axis)

    cross = np.cross(v1, v2)
    dot = np.dot(v1, v2)

    angle = np.arctan2(np.dot(cross, axis), dot)
    return np.degrees(angle)

def load_structure(pdb_id):
    file_path = os.path.join(pdb_dir, f"{pdb_id}.pdb.gz")
    if not os.path.exists(file_path):
        return None
    parser = PDBParser(QUIET=True)
    with gzip.open(file_path, "rt") as f:
        structure = parser.get_structure(pdb_id, f)
    return structure

def get_residue(structure, chain_id, resnum):
    try:
        for model in structure:
            chain = model[chain_id]
            for res in chain:
                if res.id[1] == resnum:
                    return res
    except:
        return None
    return None

def get_ca(res):
    if res and "CA" in res:
        return res["CA"].get_coord()
    return None

def get_centroid(res):
    coords = []
    for atom in res:
        if atom.get_name() not in ["N", "CA", "C", "O"]:
            coords.append(atom.get_coord())
    if not coords:
        return None
    return np.mean(coords, axis=0)


results = []

files = [f for f in os.listdir(contexts_dir) if f.endswith(".tsv")]

for file in tqdm(files, desc="Processing contexts"):

    path = os.path.join(contexts_dir, file)

    # skip empty files
    if os.path.getsize(path) == 0:
        continue

    try:
        df = pd.read_csv(path, sep="\t", header=None)
    except:
        continue

    if df.empty:
        continue

    pdb_id = df.iloc[0, 12]
    structure = load_structure(pdb_id)

    if structure is None:
        continue

    # process tripeptides (3 rows per context)
    for i in range(0, len(df), 3):

        try:
            prev = df.iloc[i]
            center = df.iloc[i+1]
            next_ = df.iloc[i+2]
        except:
            continue

        # filter center = ARG
        if center[0] != target_aa:
            continue

        # filter HHH
        if center[11] != "HHH":
            continue

        chain = center[1]

        r_prev = int(prev[2])
        r_cent = int(center[2])
        r_next = int(next_[2])

        # get residues
        res_prev = get_residue(structure, chain, r_prev)
        res_cent = get_residue(structure, chain, r_cent)
        res_next = get_residue(structure, chain, r_next)

        if not res_prev or not res_cent or not res_next:
            continue

        # get coordinates
        ca_prev = get_ca(res_prev)
        ca_cent = get_ca(res_cent)
        ca_next = get_ca(res_next)

        cen_prev = get_centroid(res_prev)
        cen_cent = get_centroid(res_cent)

        if any(x is None for x in [ca_prev, ca_cent, ca_next, cen_prev, cen_cent]):
            continue


        
        v1 = cen_prev - ca_prev      
        v2 = cen_cent - ca_cent     

       
        axis = np.cross(
            ca_cent - ca_prev,
            ca_next - ca_cent
        )

        if np.linalg.norm(axis) == 0:
            continue

        try:
            angle = signed_angle_3d(v1, v2, axis)
        except:
            continue

        # classify based on LEFT residue
        prev_aa1 = prev[9]
        size_class = size_map.get(prev_aa1, "Unknown")

        if size_class == "Unknown":
            continue

        results.append([
            pdb_id,
            prev_aa1,
            size_class,
            angle
        ])

out_df = pd.DataFrame(results, columns=[
    "pdb", "left_aa", "size_class", "angle"
])

os.makedirs(os.path.dirname(output_file), exist_ok=True)
out_df.to_csv(output_file, sep="\t", index=False)

print("Saved:", output_file)
print("Total angles:", len(out_df))