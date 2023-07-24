import os
import pandas as pd
from genome_assemblies.fasta import add_sequences

def load_gff_feature_table(filename):
	columns = {
		'assembly': str,
		'genomic_accession':str,
		'start': int,
		'end': int, 
		'strand': str,
		'name': str,
		'attributes': str
	}
	return (
		pd.read_csv(
			filename,
			sep='\t',
			compression='gzip',
			usecols=list(columns.keys()),
			dtype=columns,
			keep_default_na=False,
		)					
	)

def load_16s_records_from_gff(features_filename, length_thresholds={'min': 500, 'max': 2000}):
	match_str = '16S ribosomal RNA(?!.*methyl).*'
	# match_str = '16S ribosomal RNA.*'
	df = load_gff_feature_table(features_filename)
	features_df = df[df['name'].str.contains(match_str)].copy()
	features_df['length'] = features_df['end'] - features_df['start'] + 1
	features_df = features_df.query(
		f'length >= {length_thresholds["min"]} & length <= {length_thresholds["max"]}'
	)
	return features_df.drop(columns=['attributes'])

def create_fasta_metadata(df, delim=';'):
	return (
		'>'
		+ str(df['assembly'])
		+ delim + str(df['genomic_accession'])
		+ delim + 'loc:' + str(df['start']) + ',' + str(df['end']-1)
		+ delim + 'strand:' + df['strand']
		+ delim + 'length' + str(df['length'])
	)

def add_fasta_metadata_and_sequences(features_16s_df, fasta_filename):
	features_16s_df['metadata'] = features_16s_df.apply(create_fasta_metadata, axis=1)
	features_16s_df['sequence'] = ''
	features_16s_df = add_sequences(features_16s_df, fasta_filename)
	return features_16s_df

def write_features_to_fasta(features_16s_df, fasta_filename, overwrite=False):
	if not os.path.exists(fasta_filename) or overwrite:
		fasta_records = features_16s_df.apply(
			lambda row: row['metadata']+'\n'+row['sequence']+'\n',
			axis=1,
		).to_list()
		with open(fasta_filename, 'w') as fasta_file:
			for fasta_record in fasta_records:
				fasta_file.write(fasta_record)
