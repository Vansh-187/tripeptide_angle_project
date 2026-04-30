configfile: "config.yaml"

import glob

AA = config["target_aa"]
PDB_DIR = os.environ.get("PDB_DIR", "pdbs")

PDBS = [f.split("/")[-1].replace(".ss.out", "") for f in glob.glob("stride_out/*.ss.out")]

rule all:
    input:
        "final/angles.tsv"


rule calculate_angles:
    input:
        expand("contexts/context_for_{aa}_in_{pdb}.tsv", pdb=PDBS, aa=[AA])
    output:
        "final/angles.tsv"
    params:
        aa=AA,
        pdb_dir=PDB_DIR
    shell:
        """
        python scripts/calculate_angles.py contexts {output} {params.aa} {params.pdb_dir}
        """