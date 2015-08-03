def main():

	# load data into json 
	load_chikungunya_data()

def load_chikungunya_data():

	# initialize network data 
	mat = {}
	mat['nodes'] = {}
	mat['nodes']['row'] = []
	mat['nodes']['col'] = []

	# open file
	filename = 'virus_chikungunya/log2_heavy_light_matrix.txt'
	f = open(filename, 'r')

	lines = f.readlines()

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


	print('\n')
	print(mat['nodes']['col'])
	print('\n')
	print(mat['nodes']['row'])


		# print(inst_line)

	f.close()


# run main
main()