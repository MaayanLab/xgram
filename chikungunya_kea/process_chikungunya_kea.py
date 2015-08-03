def main():

	# load data into json 
	load_chikungunya_data()

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

					mat['nodes']['col'].append(inst_line[j])

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


	print('col')
	print(len(mat['nodes']['col']))
	print('row')
	print(len(mat['nodes']['row']))
	print('shape matrix')
	print(mat['mat'].shape)

	# convert matrix to list and save to json
	mat['mat'] = mat['mat'].tolist()

	json_scripts.save_to_json(mat,'chik_log2.json','no-indent')

# run main
main()