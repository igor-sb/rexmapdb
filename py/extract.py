#!~/anaconda/bin/python3
# -*- coding: utf-8 -*-
"""
Use downloaded GFF files (~/data/ncbi_genomes_bacteria/features) and FASTA
sequences (~/data/ncbi_genomes_bacteria/sequences) to extract annotated 16S
sequences and store them in a single FASTA file.


Created on Fri Jul 21 08:15:42 2017

@author: igor
"""

import gzip, os, sys
import argparse
import pandas as pd
from Bio import SeqIO
from _include import log

script_title = 'Extract 16S sequences from full genomes.'

def add_processed_assembly_ids(
		assembly_df,
		fasta_filename,
		fasta_delim=';',
		overwrite=False
	):
	processed_assembly_ids = set([])
	last_assid = ''
	start_index = 0
	if os.path.isfile(fasta_filename) and not overwrite:
		log('* output '+fasta_filename+' already exists: ', end='')
		log('appending to file...', end='')
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
		log('OK. Added % 5d sequences.' % (n_processed_ids))
	return (processed_assembly_ids, start_index)

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

def build_file_path(ftp_path, base_path, type):
	assert type in ['features', 'sequences'], f'Wrong file type: {type}'
	if type == 'features':
		suffix = '_feature_table.txt.gz'
	elif type == 'sequences':
		suffix = '_genomic.fna.gz'
	return os.path.join(
		base_path,
 		os.path.basename(ftp_path) + suffix,
	)

# ass_sub = '/Users/igor/cloud/research/microbiome/genomes/data/bacteria/bacteria_assembly_summary_sub.txt'
# ass_sub = '/Users/igor/cloud/research/microbiome/genomes/data/archaea/archaea_assembly_summary_sub_2017-08-24.txt'

def main(
		ass_sub,
		fa_out_name,
		ass_f_dir='features/',
		ass_s_dir='sequences/', 
		delim=';',
		match_str='16S ribosomal RNA.*',
		len_min=500,
		len_max=2000,
		overwrite=False
		# fa_out_exist_action='append'
	):
	
	"""
	Load accession summary table, and extract 16S sequences into a
	single-line FASTA. Arguments:
		- ass_sub: (filtered) assembly table, loaded as a pd.dataFrame
		- fa_out_name: name of the output FASTA file
		- ass_f_dir: directory with previously downloaded feature.gz files
		- ass_s_dir: directory with previously downloaded sequence.fa.gz files
		- delim: FASTA meta-data line delimiter to store extracted seq info
		- match_str: regex to match for name column in the GFF file, when looking for 16S seq
		- len_min, len_max: keep only 16S sequences within this range (inclusive both)
		- fa_out_exist_action: allowed values are either 'append' or 'overwrite', if overwrite
			it will overwrite the fa_out_name file if it already exists. if 'append' it will
			load the FASTA file, check what's already processed and skip those. this is used
			to resume the sdcript in case of early termination / errors.
		
	"""
	
	# ass -> assembly_df
	assembly_df = pd.read_csv(
		ass_sub,
		sep='\t',
		usecols=['#assembly_accession', 'ftp_path'],
		dtype={'#assembly_accession': str, 'ftp_path': str}
	)
	
	log('Extract 16S sequences from full genomes')
	log('* matching feature names: ' + str(match_str))
	
	#
	# If FASTA file already exists, extract all already processed entries
	# then keep them in a list so that we won't have to reopen gz files again.
	#
	
	# change fa_out_exist_action to overwrite bool
	processed_assembly_ids = set([]) # ass_exist
	start_index = 0 # i0
	write_mode = 'w'
	if not overwrite:
		processed_assembly_ids, start_index = add_processed_assembly_ids(assembly_df, fa_out_name)
		write_mode = 'a'
		assembly_df = assembly_df[~assembly_df['#assembly_accession'].isin(processed_assembly_ids)]

	# Check for missing or zero filesize features/sequence files
	assembly_df['features_path'] = assembly_df['ftp_path'].apply(build_file_path, ass_f_dir, 'features')
	assembly_df['sequences_path'] = assembly_df['ftp_path'].apply(build_file_path, ass_s_dir, 'sequences')
	
	log('* extracting sequences with lengths: '+str(len_min)+' <= len <= '+str(len_max))
	with open(fa_out_name, write_mode) as fa_out:
		log('* writing file ' + fa_out_name + '...')
		n = {
			'wrong_length': 0,
			'missing_assembly_ids': 0,
			'empty_assembly_ids': 0,
			'existing_assembly_ids': 0,
		}
		n_wrong_len = 0
		n_missing_assids = 0 # Keep track of missing records
		n_empty_assids = 0 # Keep track of the number of empty records (either GFF.GZ or FNA.GZ)
		n_existing_assids = 0 # Number of records already in the FASTA file
		
		for i in range(start_index, len(assembly_df)): 
			
			log('\r* % 5d out of % 5d assemblies (miss: % 4d, empty % 4d, % 5d seqs diff len, exist: % 4d).' % (i, len(ass),
					n_missing_assids, n_empty_assids, n_wrong_len, n_existing_assids), end='')
			
			# ass_id  -> assembly_id
			assembly_sub_df = assembly_df.iloc[i]
			assembly_id = assembly_sub_df['#assembly_accession']

			path = {
				'features': build_file_path(features_path, assembly_sub_df, 'features'),
				'sequences': build_file_path(sequences_path, assembly_sub_df, 'sequences'),
			}
			# ass_f_path -> path['features'], ass_s_path -> path['sequences']
			
			if not os.path.isfile(path['features']) or not os.path.isfile(path['sequences']):
				n_missing_assids += 1
				continue
			if os.path.getsize(path['features']) == 0 or os.path.getsize(path['sequences']) == 0:
				n_empty_assids += 1
				continue			
			try:
				ass_f = load_gff_feature_table(ass_f_path)
			except:
				n_missing_assids += 1
				continue
	
			ass_f_sub = ass_f.loc[(ass_f['name'].str.match(match_str))]
			
			# Now if we have full length 16Ss, iterate over all sequence entries
			# in FASTA sequence file to find them.            
			if not ass_f_sub.empty:
				# We found some 16S hits that are not partial (i.e. are full)
				# Now read use these coordinates to extract sequences.
				try:
					with gzip.open(ass_s_path, 'rt') as handle:
						for record in SeqIO.parse(handle, 'fasta'):
							# Check if record.id is in ass_f_sub['genomic_accession']
							# If yes, then extract the sequence based on start and end.
							ass_f_sub_sub = ass_f_sub[ass_f_sub['genomic_accession'] == record.id]
							if not ass_f_sub_sub.empty:
								for j in range(len(ass_f_sub_sub)):
									start = ass_f_sub_sub.iloc[j]['start']
									end = ass_f_sub_sub.iloc[j]['end'] + 1
									length = end - start
									# Check length
									if length < len_min or length > len_max:
										n_wrong_len = n_wrong_len + 1
									else:
										if ass_f_sub_sub.iloc[j]['strand'] == '-':
											seq = str(record.seq[start:end].reverse_complement())
										else:
											seq = str(record.seq[start:end])
										fa_out.write('>' + ass_id + delim + record.id + delim + 'loc:' + \
													 str(start) + ',' + str(end-1) + delim + 'strand:' + \
													 ass_f_sub_sub.iloc[j]['strand'] + delim + 'length:' + \
													 str(length) + '\n' + \
													 seq + '\n')
				except:
					# Error reading sequence data, skip it.
					n_missing_assids += 1
					continue
				
		# Last progress bar after we're done.        
		log('\r* % 5d out of % 5d assemblies (miss: % 4d, empty % 4d, % 5d seqs diff len, exist: % 4d).' % (i, len(ass),
			n_missing_assids, n_empty_assids, n_wrong_len, n_existing_assids))

			

def parse_input():
	parser = argparse.ArgumentParser(
		description=script_title,
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
	)
	parser.add_argument('--input-asssum', required=True, 
						help='Input assembly summary file (either bacteria or archaea).')
	parser.add_argument('--input-features-dir', required=True, 
						help='Folder with assembly features (typically data/features/).')
	parser.add_argument('--input-sequences-dir', required=True, 
						help='Folder with assembly sequences (typically data/sequences/).')
	parser.add_argument('--output-fasta', required=True, 
						help='Output FASTA with 16S sequences.')
	parser.add_argument('--action', required=False, default='append',
						help='Whether to "append" or "overwrite" output FASTA. Used for resuming from errors.')
	args = parser.parse_args()
	return args


#%% run this is script is run directly and not imported
if __name__ == '__main__':
#    ass_f_dir = '/Users/igor/data/ncbi_genomes/archaea/features/'
#    ass_s_dir = '/Users/igor/data/ncbi_genomes/archaea/sequences/'
#    out_fa = '/Users/igor/cloud/research/microbiome/genomes/data/archaea/16s_from_genomes_2017-08-24.fasta'
#    in_file = '/Users/igor/cloud/research/microbiome/genomes/data/archaea/archaea_assembly_summary_sub_2017-08-24.txt'
#    ass_f_dir = '/Users/igor/cloud/research/microbiome/genomes/data/bacteria/features/'
#    ass_s_dir = '/Users/igor/cloud/research/microbiome/genomes/data/bacteria/sequences/'
#    out_fa = '/Users/igor/cloud/research/microbiome/genomes/data/bacteria/16s_from_genomes_2017-07-20.fasta'
#    in_file = '/Users/igor/cloud/research/microbiome/genomes/data/bacteria/bacteria_assembly_summary_sub_2017-07-20.txt'
	# in_file = sys.argv[1]
	# ass_f_dir = sys.argv[2]
	# ass_s_dir = sys.argv[3]
	# out_fa = sys.argv[4]
	# if len(sys.argv) == 6:
	#     action = sys.argv[5]
	# else:
	#     action = 'append'
	log(script_title)
	args = parse_input()
	
	in_file = args.input_asssum
	ass_f_dir = args.input_features_dir
	ass_s_dir = args.input_sequences_dir
	out_fa = args.output_fasta
	action = args.action
	
	ass_sub = in_file
	fa_out_name = out_fa
	main(in_file, out_fa, ass_f_dir=ass_f_dir, 
		 ass_s_dir=ass_s_dir, fa_out_exist_action=action)
