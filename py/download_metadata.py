import os
import fire
from datetime import datetime
import logging
from six.moves import urllib

logging.basicConfig(level=logging.INFO)
LOG = logging.getLogger(__name__)

def download_kingdom_assembly_summary(kingdom, output_filename):
	assembly_summary_base_url = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/'
	assembly_summary_full_url = f'{assembly_summary_base_url}/{kingdom}/assembly_summary.txt'

	output_path = os.path.dirname(output_filename)
	if not os.path.isdir(output_path):
		os.makedirs(output_path)
	
	urllib.request.urlretrieve(assembly_summary_full_url, output_filename)
	
def download_assembly_summaries(output_path, kingdoms=['archaea', 'bacteria'], timestamp=True):
	suffix = '_'+datetime.now().strftime('%Y-%m-%d') if timestamp else ''
	for kingdom in kingdoms:
		output_filename = os.path.join(output_path, f'{kingdom}_assembly_summary{suffix}.txt')
		download_kingdom_assembly_summary(kingdom, output_filename)

if __name__ == '__main__':
	fire.Fire(download_assembly_summaries)
