import gzip
from Bio import SeqIO

def extract_sequence_from_fasta_record(fasta_record, row):
	sequence = fasta_record.seq[row['start']:row['end']]
	if row['strand'] == '-':
		sequence = sequence.reverse_complement()
	return str(sequence)

def add_sequences(df, sequences_filename):
	with gzip.open(sequences_filename, 'rt') as sequences_fasta:
		for fasta_record in SeqIO.parse(sequences_fasta, 'fasta'):
			fasta_rows = df['genomic_accession'] == fasta_record.id
			df.loc[fasta_rows, 'sequence'] = df.loc[fasta_rows].apply(
				lambda row: extract_sequence_from_fasta_record(fasta_record, row),
				axis=1,
			)
	return df
