import fire
import pandas as pd
from genome_assemblies.summary import (
    load_assembly_summary_table,
    cleanup_infraspecific_names,
    strain_name_from_organism_name,
    add_detailed_strain_designation,
)

def filter_assembly_summary(assembly_summary_filename):
	assembly_levels = ['Complete Genome', 'Chromosome', 'Scaffold', 'Contig']
	refseq_levels = ['reference genome', 'representative genome', 'na']
	
	as_df = load_assembly_summary_table(assembly_summary_filename)	
	as_df['infraspecific_name'] = as_df['infraspecific_name'].apply(cleanup_infraspecific_names)
	as_df['strain_name'] = as_df['organism_name'].apply(strain_name_from_organism_name)
	as_df['assembly_level'] = pd.Categorical(
		as_df['assembly_level'],
		categories=assembly_levels,
		ordered=True,
	)
	as_df['refseq_category'] = pd.Categorical(
		as_df['refseq_category'],
		categories=refseq_levels,
		ordered=True,
	)
	as_df['strain_name'] = as_df.apply(add_detailed_strain_designation, axis=1)
	as_df = as_df.sort_values(
		by=['assembly_level', 'refseq_category', 'seq_rel_date'],
		ascending=[False, False, False]
	)
	return (
		as_df
		.groupby('strain_name')
		.first()
		.reset_index(drop=True)
	)

if __name__ == '__main__':
	fire.Fire(filter_assembly_summary)
