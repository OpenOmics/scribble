import argparse

parser = argparse.ArgumentParser(description='create protein coding GTF removing transcript versions')
parser.add_argument('-i', dest='input_gtf', help='GTF file to be filtered', required=True)
parser.add_argument('-o', dest='output_gtf', help='Filtered GTF file name', required=True)
args = parser.parse_args()

with open(args.input_gtf, 'r') as input_file:
    with open(args.output_gtf, 'w') as output_file:
        for line in input_file:
            if line.startswith('#'):
                output_file.write(line)
            else:
                if 'gene_type "protein_coding"' in line:
                    gene_id = line.strip().split('gene_id "')[1].split('";')[0]
                    line = line.replace(gene_id, gene_id.split('.')[0])
                    if 'transcript_id "' in line:
                        transcript_id = line.strip().split('transcript_id "')[1].split('";')[0]
                        line = line.replace(transcript_id, transcript_id.split('.')[0] if '.' in transcript_id else transcript_id)
                    output_file.write(line)
