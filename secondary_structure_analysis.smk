import os
import glob

STRIDE_BIN = os.environ.get("STRIDE_BIN", "stride")
PDBS = [os.path.basename(f).replace(".pdb.gz", "") for f in glob.glob("pdbs/*.pdb.gz")]

rule all:
    input:
        expand("stride_out/{pdb}.ss.out", pdb=PDBS)

rule unzip_pdb:
	input: 
		pdb = "pdbs/{pdb}.pdb.gz"
	output: 
		unzipped = temp("unzipped_pdbs/{pdb}.pdb")
	shell: 
		" zcat {input.pdb} > {output.unzipped} "
rule run_stride:
	input: 
		unzipped = "unzipped_pdbs/{pdb}.pdb"
	output:
		ss_out = "stride_out/{pdb}.ss.out"
	params:
	shell:
		f"{STRIDE_BIN} {{input.unzipped}} > {{output.ss_out}}"

