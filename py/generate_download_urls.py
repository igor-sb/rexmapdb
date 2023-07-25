	# Write filtered list to file (we will use this in other scripts)
	print('- Writing clean assembly summary table...', end='')
	ass_filter_out = os.path.join(out_path, k+'_assembly_summary_filter_'+date+'.txt')
	ass_df_filter.to_csv(ass_filter_out, sep='\t', index=False)
	print('OK.')
	
	# Save a separate list of feature and sequences for download
	fea_out = os.path.join(out_path, k+'_assembly_features_urls.txt')
	seq_out = os.path.join(out_path, k+'_assembly_sequences_urls.txt')

	with open(fea_out, 'w') as f_out, open(seq_out, 'w') as s_out:
		for i, row in enumerate(ass_df_filter.itertuples()):
			f_out.write(os.path.join(row[8], os.path.basename(row[8])+'_feature_table.txt.gz\n'))
			s_out.write(os.path.join(row[8], os.path.basename(row[8])+'_genomic.fna.gz\n'))
			if i % 100 == 0 or i == len(ass_df_filter)-1:
				print('\r- Processed % 5d out of % 5d assembly records.' % (i+1, len(ass_df_filter)), end='')
	print('\n- Saved feature and sequence URLs.')
