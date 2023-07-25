import re
import os
import pandas as pd
from Bio import SeqIO

def file_exists(filename):
	return os.path.exists(filename) and os.path.getsize(filename) > 0

def filter_missing_files(assembly_df):
	assembly_df['features_exist'] = assembly_df['features_path'].apply(file_exists)
	assembly_df['sequences_exist'] = assembly_df['sequences_path'].apply(file_exists)
	return (
		assembly_df
		.query('features_exist == True & sequences_exist == True')
		.drop(columns=['features_exist', 'sequences_exist'])
	)

def read_filtered_assembly_summary(filename):
	columns = {'#assembly_accession': str, 'ftp_path': str}
	return pd.read_csv(
		filename,
		sep='\t',
		usecols=list(columns.keys()),
		dtype=columns
	)

def build_paths(ftp_path, base_paths):
	path_prefix = os.path.basename(ftp_path)
	features_paths = os.path.join(base_paths['features'], path_prefix + '_feature_table.txt.gz')
	sequences_paths = os.path.join(base_paths['sequences'], path_prefix + '_genomic.fna.gz')
	return pd.Series(
		[features_paths, sequences_paths],
		index=['features_path', 'sequences_path']
	)

def add_features_sequences_paths(assembly_df, features_base_path, sequences_base_path):
	assembly_df[['features_path', 'sequences_path']] = assembly_df['ftp_path'].apply(
		build_paths,
		base_paths = {'features': features_base_path, 'sequences': sequences_base_path}
	)
	return assembly_df

def find_processed_assembly_ids(
		assembly_df,
		fasta_filename,
		fasta_delim=';',
		overwrite=False
	):
	processed_assembly_ids = set([])
	last_assid = ''
	start_index = 0
	if os.path.isfile(fasta_filename) and not overwrite:
		with open(fasta_filename, 'r') as fasta_out:
			n_processed_ids = 0
			for record in SeqIO.parse(fasta_out, 'fasta'):
				last_assid = record.id.strip().split(fasta_delim)[0]
				processed_assembly_ids |= set([last_assid])
				n_processed_ids += 1
		try:
			last_checked_index = assembly_df['#assembly_accession'] == last_assid
			start_index = assembly_df.index[last_checked_index].tolist()[0]                    
		except:
			start_index = 0
	return processed_assembly_ids

def load_assembly_summary_table(filename):
	columns = {
		'#assembly_accession': str,
		'refseq_category': str,
		'taxid': str,
		'organism_name': str,
		'infraspecific_name': str,
		'assembly_level': str,
		'seq_rel_date': str,
		'ftp_path': str
	}
	return pd.read_csv(
		filename,
		delimiter='\t',
		skiprows=1,
		usecols=list(columns.keys()), 
		dtype=columns
	)

def cleanup_infraspecific_names(infraspecific_name):
	return (
		str(infraspecific_name)
		.replace('strain=', '')
		.replace('nan', '')
		.replace('substr. ', '')
	)

def strain_name_from_organism_name(organism_name):
	""" Return genus species strain (if any) name from full genome species name text. 
	Sometimes, we have subspecies name exact same as species, so we remove that. """
	strain_name = (
		str(organism_name)
		.replace('substr. ', '')
		.replace('str. ', '')
		.replace('subsp. ', '')
		.replace('\'', '')
	)
	strain_name = (
		re.sub(r'\b(\w+)( \1\b)+', r'\1', strain_name)
		.replace('[', '')
		.replace(']', '')
	)
	return re.sub(r'\b(\w+)( \1\b)+', r'\1', strain_name)

def is_infraspecific_not_in_strain_name(strain_name, infraspecific_name):
	return (
		infraspecific_name not in strain_name and \
			infraspecific_name.replace(' ', '') not in strain_name and \
			'(' not in strain_name
	)

def add_detailed_strain_designation(row):
	if is_infraspecific_not_in_strain_name(
		row['strain_name'],
		row['infraspecific_name']
	):
		return row['strain_name'] + ' ' + row['infraspecific_name']
	else:
		return row['strain_name']	
