def main():

	# load data into json 
	mat = load_chikungunya_data()

	clean_chik_data(mat)

def clean_chik_data(mat):
	import json_scripts
	import numpy as np
	import scipy
	from collections import defaultdict


	# initailize cleaned chik 
	chik = {}
	chik['nodes'] = {}
	# transfer columns 
	chik['nodes']['col'] = mat['nodes']['col']
	chik['nodes']['row'] = []
	chik['mat'] = []

	# loop through the rows and check which terms appear more than once 
	# count the occurence of terms in list 
	d = defaultdict(int)
	for word in mat['nodes']['row']:
		d[word] += 1

	# get unique terms 
	all_rows = d.keys()
	all_rows = sorted(all_rows)
	print('there are '+str(len(mat['nodes']['row']))+' non-redundant proteins - in the original data')
	print('there are '+str(len(all_rows))+ ' unique proteins')

	# loop through the rows 
	for inst_row in all_rows:

		# append row to chik['nodes']
		chik['nodes']['row'].append(inst_row)

		if d[inst_row]==1:
			# get the original index of the row 
			inst_x = mat['nodes']['row'].index(inst_row)
			inst_values = mat['mat'][inst_x,:] 

			# add data to matrix 
			if len(chik['mat']) == 0:
				chik['mat'] = inst_values
			else:
				chik['mat'] = np.vstack((chik['mat'],inst_values))

		else:
			tmp_x = [i for i, x in enumerate(mat['nodes']['row']) if x == inst_row]

			print('there are ' + str(len(tmp_x)) + ' copies of ' + inst_row)

			# initialize values 
			num_cols = mat['mat'].shape[1]
			inst_values = scipy.zeros([  num_cols ])

			# make inst_values a sum of all vectors 
			for i in range(len(tmp_x)):
				tmp_values = mat['mat'][tmp_x[i],:]

				# sum over the values in inst_values 
				for j in range(len(tmp_values)):

					inst_values[j] = np.nansum([inst_values[j],tmp_values[j]])

			# add values to mat 
			if len(chik['mat']) == 0:
				chik['mat'] = inst_values
			else:
				chik['mat'] = np.vstack((chik['mat'],inst_values))

	# print(chik['mat'])
	print('there are '+str(len(chik['nodes']['row']))+' proteins in chik nodes' )
	print(chik['mat'].shape)

	# convert matrix to list and save to json
	chik['mat'] = chik['mat'].tolist()

	json_scripts.save_to_json( chik,'chik_log2.json','no-indent')	

def load_chikungunya_data():
	import json_scripts
	import numpy as np

	# initialize network data 
	mat = {}
	mat['nodes'] = {}
	mat['nodes']['row'] = []
	mat['nodes']['col'] = []

	# open file
	filename = 'virus_chikungunya/log2_heavy_light_matrix.txt'
	f = open(filename, 'r')
	lines = f.readlines()
	f.close()

	for i in range(len(lines)):
		inst_line = lines[i].strip().split('\t')

		# get column labels from first line 
		if i == 0:

			# loop through the columns
			for j in range(len(inst_line)):

				# skip first column label (not a label)
				if j > 0:

					# process the column names 
					inst_label = inst_line[j].replace('|','').replace('Skeletal','').replace('Muscle','').replace('cell','').replace('expt','')

					mat['nodes']['col'].append(inst_label)

		else:

			# ignore proteins with no name 
			if inst_line[0] != '-':

				# get inst_name
				inst_name = inst_line[0]

				# correct sept protein renaming (excel)
				if 'Sep' in inst_name:

					inst_num = inst_name.split('-')[0]

					inst_name = 'SEPT'+inst_num

				# get the row labels 
				mat['nodes']['row'].append(inst_name)

				# save values 
				###################
				inst_values = inst_line[1:]

				# transfer values to floats
				inst_values = [float(x) for x in inst_values]

				# transfer to numpy array 
				inst_values = np.asarray(inst_values)

				# transfer values to matrix 
				if i == 1:
					mat['mat'] = inst_values
				else:
					mat['mat'] = np.vstack((mat['mat'],inst_values))


	# print('col')
	# print(len(mat['nodes']['col']))
	# print('row')
	# print(len(mat['nodes']['row']))
	# print('shape matrix')
	# print(mat['mat'].shape)

	# return the matrix - which has repeats for protein measurements 
	return mat 



# run main
main()