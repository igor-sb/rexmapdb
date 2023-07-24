import math

def partition_dataframe(df, n_chunks):
	size = math.ceil(len(df) / n_chunks)
	partitioned_dfs = []
	for start in range(0, len(df), size):
		end = min(start + size, len(df))
		partitioned_dfs.append(df.iloc[start:end].copy())
	return partitioned_dfs

