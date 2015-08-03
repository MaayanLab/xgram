def main():
	print('here')

	# extract some data of protein type and zscore 
	extract_prot_type_z_score('TF',2)

def extract_prot_type_z_score( inst_type, z_cutoff ):
	import json_scripts
	import d3_clustergram
	import numpy as np

	# load ccle data with zscore normalization
	ccle = json_scripts.load_to_dict('CCLE/nsclc_allzc.json')
	# convert zscored data to nparray 
	ccle['data_z'] = np.asarray(ccle['data_z'], dtype = float)

	# load gs_list to filter for kinase, tf etc
	# gs_list = json_scripts.load_to_dict('enz_and_tf_lists_gmts/categories_gs_list.json')
	gs_list = json_scripts.load_to_dict('gene_classes_harmonogram.json')

	# generate node lists 
	nodes = {}
	# get all genes 
	nodes['row'] = filter_genes(ccle, gs_list, inst_type, z_cutoff)
	# get all cell lines from CCLE 
	nodes['col'] = ccle['cell_lines']

	# minimum number intersect
	min_num_int = 3

	# only make clustergram if there are genes remaining after filtering 
	# the minimum needed to produce a clustergram 
	if len(nodes['row']) > 1:

		print('there are ' + str(len(nodes['row'])) + ' genes being clustered' )

		# Generate data_mat: used to filter data from the original ccle for a subset of genes 
		# takes inputs: node lists for rows and columns, and primary data that will be used to make the matrix
		# the last two arguments are the names of the rows and columns in the original data
		data_mat = d3_clustergram.generate_data_mat_array( nodes, ccle, 'gene', 'cell_lines', 'data_z' )

		# # cluster rows and columns 
		# clust_order = d3_clustergram.cluster_row_and_column( nodes, data_mat, 'cosine', min_num_int )

	# open file to write tsv
	fw = open('CCLE/example_tsv','w')

	# write network to file 
	###########################
	# write column names to first row 
	# write first tab
	fw.write('\t')
	# loop over column names 
	for inst_col in nodes['col']:
		fw.write(inst_col+'\t')

	# write new line
	fw.write('\n')

	# get the number of rows 
	num_rows = len(nodes['row']) - 1

	num_cols = len(nodes['col']) - 1

	# loop through rows 
	for i in range(len(nodes['row'])):

		# get inst_row - the inst row name 
		inst_row = nodes['row'][i]

		# write row name 
		fw.write(inst_row+'\t')

		# get row data 
		row_data = data_mat[i,:].tolist()

		# write row data 
		for j in range(len(row_data)):

			# get inst data point 
			inst_point = row_data[j]

			fw.write(str(inst_point))

			if j < num_cols:
				fw.write('\t')

		print(i)

		if i < num_rows:
			# write new line 
			fw.write('\n')

	fw.close()

# gather all genes and cell lines 
def filter_genes(ccle, gs_list, inst_type, z_cutoff):
	import numpy as np 

	# get all gene names 
	all_genes = ccle['gene']

	# define the list of all genes as all the genes in ccle
	gs_list['all'] = all_genes 

	# filter for zscore here 
	#########################
	filtered_genes = []

	# loop through the genes and check if they are part of the gene class of interest 
	for inst_gene in all_genes:

		# check if gene of the type of interest 
		if inst_gene in gs_list[inst_type]:
			
			# get index of gene 
			inst_index =  all_genes.index(inst_gene)

			# get all zscores 
			inst_zscores = ccle['data_z'][inst_index,:]

			# get maximum zscore
			inst_max_z = np.amax(np.absolute(inst_zscores)) 

			# check if gene has zscore above cutoff 
			if inst_max_z > z_cutoff:
				filtered_genes.append(inst_gene)

	return filtered_genes

# run main
main()