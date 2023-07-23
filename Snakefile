# Snakefile
from datetime import datetime

current_date = datetime.now().strftime("%Y-%m-%d")

data_path = "data"
seqs_path = "{data_path}/sequences"
feat_path = "{data_path}/features"

# Define input files
archaea_assembly_summary = "{data_path}/archaea_assembly_summary_{current_date}.txt"
bacteria_assembly_summary = "{data_path}/bacteria_assembly_summary_{current_date}.txt"

rule download_assembly_summaries:
    shell:
        "python3 py/dl_genome_meta.py {input}"

rule download_features:
	shell:
		"""
		mkdir {feat_path}
		cd {feat_path}

		python3 ../../py/download.py ../archaea_assembly_features_urls.txt 10
		python3 ../../py/download.py ../bacteria_assembly_features_urls.txt 10

		cd ../../
		"""

rule download_sequences:
	shell:
		"""
		mkdir {seqs_path}
		cd data/sequences

		python3 ../../py/download.py ../archaea_assembly_sequences_urls.txt 10
		python3 ../../py/download.py ../bacteria_assembly_sequences_urls.txt 10

		cd ../../
		"""
