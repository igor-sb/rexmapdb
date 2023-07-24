import os
import fire
import logging
import multiprocessing

from genome_assemblies.gff import (
    load_16s_records_from_gff,
    add_fasta_metadata_and_sequences,
    write_features_to_fasta,
)
from genome_assemblies.summary import (
    filter_missing_files,
    read_filtered_assembly_summary,
    add_features_sequences_paths,
)

from utils import partition_dataframe

logging.basicConfig(level=logging.INFO)
LOG = logging.getLogger(__name__)

def extract_16s_from_assembly_chunk(assembly_chunk_df, output_base_path):
	for assembly_index in range(len(assembly_chunk_df)):
		assembly_record = assembly_chunk_df.iloc[assembly_index]
		features_filename = assembly_record['features_path']
		sequences_filename = assembly_record['sequences_path']
		output_fasta_filename = os.path.join(
			output_base_path,
			assembly_record['#assembly_accession'] + '.fasta',
		)
		if assembly_index % 100 == 0:
			LOG.info(f' id: {assembly_index}, record: {assembly_record["#assembly_accession"]}')
		try:
			features_16s_df = load_16s_records_from_gff(features_filename)
			features_16s_df = add_fasta_metadata_and_sequences(features_16s_df, sequences_filename)
		except EOFError:
			continue
		if features_16s_df.empty:
			continue
		write_features_to_fasta(features_16s_df, output_fasta_filename)

def main(
		assembly_summary_filename,
		output_base_path='data/sequences_16s/',
		features_base_path='data/features/',
		sequences_base_path='data/sequences/',
		nprocs = 10,
	):
	LOG.info(f'Extract 16S sequences from full genomes ({nprocs} processes)')
	if not os.path.exists(output_base_path):
		os.makedirs(output_base_path)
		LOG.info(f' created output {output_base_path}.')
	assembly_df = read_filtered_assembly_summary(assembly_summary_filename)
	assembly_df = add_features_sequences_paths(assembly_df, features_base_path, sequences_base_path)
	assembly_df	= filter_missing_files(assembly_df).reset_index()
	LOG.info(f' total: {len(assembly_df)} assemblies')

	assembly_df_chunked = partition_dataframe(assembly_df, n_chunks=nprocs)
	LOG.info(' chunk sizes: %s', [len(_) for _ in assembly_df_chunked])
	procs = []
	for assembly_chunk_df in assembly_df_chunked:
		proc = multiprocessing.Process(
			target=extract_16s_from_assembly_chunk,
			args=(assembly_chunk_df, output_base_path),
		)
		procs.append(proc)
	
	LOG.info(f' starting {nprocs} processes')
	for proc in procs:
		proc.start()
	
	LOG.info(f' joining {nprocs} processes')
	for proc in procs:
		proc.join()

if __name__ == '__main__':
	fire.Fire(main)
