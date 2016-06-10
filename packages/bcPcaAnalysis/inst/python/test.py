import json
import urllib2


def split_edge_list(edge_list,split_char='\t'):
	new_edge_list = [[],[]]

	for edge in edge_list:
		genes = edge.split(split_char)
		#print genes
		new_edge_list[0].append(genes[0])
		new_edge_list[1].append(genes[1])

	return new_edge_list


def file_to_list(my_file):
	'''Converts each line in the given filename (str) to an item on a list, returns the list'''
	my_file = open(my_file)
	returned_list = []
	for line in my_file:
		line = line.strip()
		returned_list.append(line)
	my_file.close()
	return returned_list



def convert_go_array_to_edgewise(go_array,universe_edge_list):
	'''Given a PPI list corresponding to the PPIs in the universe (universe_edge_list)
	converts a go_array object e.g. ['GO:0004386', 'YAL019W, 'YBL023C'] to interactions
	within that go term e.g. ['GO:0004386', 'YAL019W YBL023C'] and returns that object
	'''

	go_dict = {}
	for item in go_array:
		go_dict[item[0]] = item[1:]

	inverse_go_dict = {}
	for go_term in go_dict:
		for gene in go_dict[go_term][1:]:
			if not inverse_go_dict.get(gene):
				inverse_go_dict[gene] = []
			inverse_go_dict[gene].append(go_term)

	inverse_pairwise_go_dict = {}
	for i in range(len(universe_edge_list[1])):
		gene1 = universe_edge_list[0][i]
		gene2 = universe_edge_list[1][i]
		ppi = gene1 + ' ' + gene2

		terms_gene1 = inverse_go_dict.get(gene1)
		terms_gene2 = inverse_go_dict.get(gene2)
		if terms_gene1 and terms_gene2:
			common_go_terms = list(set(terms_gene1).intersection(set(terms_gene2)))
			inverse_pairwise_go_dict[ppi] = common_go_terms



	pairwise_go_dict = {}
	for ppi in inverse_pairwise_go_dict:
		for term in inverse_pairwise_go_dict[ppi]:
			if not pairwise_go_dict.get(term):
				pairwise_go_dict[term] = []
			pairwise_go_dict[term].append(ppi)




	pairwise_go_array = []

	for item in pairwise_go_dict:
		pairwise_go_dict[item].insert(0,item)
		pairwise_go_array.append(pairwise_go_dict[item])

	return pairwise_go_array


def parse_go_terms_slim(go_slim_file,all_gene_list,mapping_index=1):
	'''go_slim_file is the name (str) of a tab-separated go slim file with the following structure
	YAL055W	PEX22	S000000051	P	peroxisome organization	GO:0007031	ORF|Verified
	
	all_gene_list and a list of genes in the 'gene universe'which is one SGD
	ID per entry

	If SGD Common is given in all_gene_list, set mapping_index to 1
	If SGD ORF is given in all_gene_list, set mapping_index to 0


	returns a two-item list
		1) a list of lists, wherein each list is a GO ID
	followed by the SGD COMMON genes it contains e.g. ['GO:0004386', 'DHH1', 'PRP43']
		2) a dictionary where every key is a go term and every
	value is a description of that go term e.g. {''GO:0004386':helicase activity}
	'''
	go_slim_array_precursor = {}
	go_slim_dict = {}

	go_slim_file = open(go_slim_file)
	for line in go_slim_file:
		line = line.strip().split('\t')
		go_slim_dict[line[-2]] = line[-3]
		if 'ORF' in line[-1] and (line[mapping_index] in all_gene_list or line[0] in all_gene_list):
			if not go_slim_array_precursor.get(line[-2]):
				go_slim_array_precursor[line[-2]] = []
			if len(line[1]) > 2 and line[mapping_index] in all_gene_list:
				go_slim_array_precursor[line[-2]].append(line[mapping_index])
			else:
				go_slim_array_precursor[line[-2]].append(line[0])

	go_slim_array = []
	go_slim_file.close()

	for item in go_slim_array_precursor:
		go_slim_array_precursor[item].insert(0,item)
		go_slim_array.append(go_slim_array_precursor[item])
	
	return [go_slim_array,go_slim_dict]


def parse_go_terms_funcassociate(go_funcassociate_file,all_gene_list):
	'''go_funcassociate_file is the name of a GO file downloaded from funcassociate with the following format:
	#Series of comment headers
	Then the GO term, its definition, and list of SGD ORFs corresponding to that
	term, for example
	GO:0000001	mitochondrion inheritance	YAL048C YDL006W

	returns a two-item list
		1) a list of lists, wherein the first item of each list is a GO ID
	followed by the SGD COMMON genes it contains e.g. ['GO:0004386', 'YAL019W, 'YBL023C']
		2) a dictionary where every key is a go term and every
	value is a description of that go term e.g. {''GO:0004386':helicase activity}
	'''
	
	go_funcassociate_array_precursor = {}
	go_funcassociate_dict = {}

	go_funcassociate_file = open(go_funcassociate_file)
	line = go_funcassociate_file.readline().strip()
	while line:
		if not line.startswith('#'):
			line_list = line.strip().split('\t')
			go_funcassociate_dict[line_list[0]] = line_list[1]
			genes = line_list[2].split(' ')
			if not go_funcassociate_array_precursor.get(line_list[0]):
				go_funcassociate_array_precursor[line_list[0]] = []
			for gene in genes:
				go_funcassociate_array_precursor[line_list[0]].append(gene)
		line = go_funcassociate_file.readline().strip()
	go_funcassociate_file.close()

	go_funcassociate_array = []

	for item in go_funcassociate_array_precursor:
		go_funcassociate_array_precursor[item].insert(0,item)
		go_funcassociate_array.append(go_funcassociate_array_precursor[item])

	return [go_funcassociate_array,go_funcassociate_dict]
	

def funcassociate_server_submission(given_gene_list,\
	association_array,\
	attribute_dict,\
	order_mode="ordered",\
	funcassociate_url='http://llama.mshri.on.ca/cgi/funcassociate/serv',\
	p_value_cutoff=0.05):
	'''given_gene_list is a list of terms e.g. ['DHH1','PRP43']

	association_array is a list of lists, wherein each list is a GO ID
	followed by the terms it contains e.g. ['GO:0004386', 'DHH1', 'PRP43']

	attribute_dict is a dictionary where every key is a go term and every
	value is a description of that go term

	order mode can be set to either "ordered" or "unordered" as a string

	funcassociate_url gives the API server

	returns the server response as a JSON string
	'''

	url = 'http://llama.mshri.on.ca/cgi/funcassociate/serv'

	values =  {
			  "id": "dc8bed03b5a9787447db75d7b948784f",
			  "method": "functionate",
			  "params": [
			              {
			                "query": given_gene_list,
			                "mode":order_mode,
							"associations":association_array,
							"attrib_dict":attribute_dict,
							"cutoff":p_value_cutoff
			              }
			            ],
			  "jsonrpc": "2.0"
			}
	
	response = urllib2.Request(url, json.dumps(values))
	response.add_header('Content-Type', 'application/json')


	response = urllib2.urlopen(response)
	response = response.read()
	return(response)


def funcassociate(gene_file,\
	universe_file,\
	go_association_file,\
	order_mode="ordered",\
	go_type='slim',\
	nametype='common',\
	mode='nodewise'):
	'''gene_file is a filename (str) of a file containg one item per line of the genes of interest

	universe_file is a filename (str) of a file containg one item per line of all possible genes of
	interest in the experiment

	go_association_file is a filename(str) of a file which maps each term in universe_file to a GO term and description

	order_mode is a string, either "ordered" or "unordered"

	go_type is a string, either "slim" or "default_funcassociate", indicating how the GO mapping file should be processed
		with the given data, slim files use COMMON names and default_funcassociate uses SGD ORF names, gene lists should
		be processed accordinly'''



	#Check for validity
	if nametype == 'common' and go_type == 'default_funcassociate':
		print('Combination of name type and GO type not supported')
		exit()
	elif nametype == 'common':
		mapping_index_slim = 1
	elif nametype == 'ORF':
		mapping_index_slim = 0
	else:
		print('Name type not supported')
		exit()

	if not go_type in ['slim','default_funcassociate']:
		print('Go type not supported')
		exit()

	gene_list = file_to_list(gene_file)
	universe_gene_list = file_to_list(universe_file)

	if mode == 'edgewise':
		edge_list = split_edge_list(gene_list)
		collapsed_edge_list = edge_list[0] + edge_list[1]
		gene_list = collapsed_edge_list

		universe_edge_list = split_edge_list(universe_gene_list)
		collapsed_universe_edge_list = universe_edge_list[0] + universe_edge_list[1]
		universe_gene_list = collapsed_universe_edge_list


	if go_type == 'slim':
		go_conversion = parse_go_terms_slim(go_association_file,universe_gene_list,mapping_index_slim)
	elif go_type=='default_funcassociate':
		go_conversion = parse_go_terms_funcassociate(go_association_file,universe_gene_list)



	go_array = go_conversion[0]
	#print go_array
	go_dict = go_conversion[1]

	if mode == 'edgewise':
		go_array = convert_go_array_to_edgewise(go_array,universe_edge_list)
		new_gene_list = []
		print edge_list
		for i in range(len(edge_list[0])):
			new_gene_list.append(edge_list[0][i] + ' ' + edge_list[1][i])
		gene_list = new_gene_list
		print gene_list

	return funcassociate_server_submission(gene_list,go_array,go_dict,order_mode)

if __name__ == '__main__':
	
	#print funcassociate('gene_space_enhanced.txt',
	#	'gene_space_total.txt',
	#	'go_slim_mapping.tab',
	#	go_type='slim',
	#	nametype='common')

	print funcassociate('enhanced_inters.tsv',
		'ppi_universe_orfs.tsv',
		'funcassociate_go_associations.txt',
		go_type='default_funcassociate',
		nametype='ORF',
		mode='edgewise')