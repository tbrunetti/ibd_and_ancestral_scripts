import argparse
import pandas
import numpy

# A collection of ancestry functions

'''
func: unrelated_inds(ibdMat, pihat)
input:
output:
'''
def unrelated_inds(ibdMat, pihat):
	unrelated_ids = []
	ibd = pandas.read_table(ibdMat, delim_whitespace=True)
	unique_inds = list(set(list(ibd['IID1'])))
	for ind in unique_inds:
		get_ind = ibd.loc[((ibd['IID1'] == ind))]
		if max(list(get_ind['PI_HAT'])) <= pihat:
			unrelated_ids.append(str(ind))

	print len(unrelated_ids)


'''
func: nearestFounder(pedigree)
input:
output:
'''
def nearestFounder(pedigree):
	import networkx as nx
	
	ancestryGraph = nx.Graph()	
	pedigreeInfo = pandas.read_table(pedigree, dtype=str)
	
	# TODO: 
	# add an exception for any individual id that is not unique in pedigreeInfo
	
	trueFoundersInfo = pedigreeInfo.loc[(pedigreeInfo['father'] == "0") & (pedigreeInfo['mother'] == "0")]
	nonFoundersInfo = pedigreeInfo.loc[(pedigreeInfo['father'] != "0") & (pedigreeInfo['mother'] != "0")]
	fatherFounderInfo = pedigreeInfo.loc[(pedigreeInfo['father'] == "0") & (pedigreeInfo['mother'] != "0")]
	motherFounderInfo = pedigreeInfo.loc[(pedigreeInfo['father'] != "0") & (pedigreeInfo['mother'] == "0")]

	founders = list(trueFoundersInfo['ind'])
	semiFounders = list(fatherFounderInfo['ind']) + list(motherFounderInfo['ind'])
	nonFounders = list(nonFoundersInfo['ind'])

	# TODO:
	# add check to make sure founders, semiFounders and nonFounders lists do no overlap by any ids;

	ancestryGraph.add_nodes_from(founders, status="founder") # adds founder individual nodes
	ancestryGraph.add_nodes_from(semiFounders, status="semiFounder")
	ancestryGraph.add_nodes_from(nonFounders, status="nonFounder") # adds non-founder individual nodes

	# TODO:
	# get tuples of all connections for non-founder nodes only and those with only one founder parent(true founders will have no mother or father connections)

	edges_father = list(zip(nonFoundersInfo.ind, nonFoundersInfo.father))
	edges_mother = list(zip(nonFoundersInfo.ind, nonFoundersInfo.mother))
	edges_semi = list(zip(fatherFounderInfo.ind, fatherFounderInfo.mother)) + list(zip(motherFounderInfo.ind, motherFounderInfo.father))

	ancestryGraph.add_edges_from(edges_father)
	ancestryGraph.add_edges_from(edges_mother)
	ancestryGraph.add_edges_from(edges_semi)




if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Set of methods that help get ancestry information from IBD estimates")
	parser.add_argument('-method', required=True, help="OPTIONS: unrelated, findFounder")
	parser.add_argument('--ibdMat', help="Path to IBD matrix from PLINK --genome IBD estimates")
	parser.add_argument('--maxPIHAT', default=0.05, type=float, help="Maximum pi-hat value in order to be considered unrelated to individual (range 0-1.0, inclusive)")
	parser.add_argument('--ped', type=str, help="tab-delimited pedigree information of all samples with following headers: ind, father, mother; founder individuals shoud have mother and/or father set to 0")

	args = parser.parse_args()

	if args.method.lower() == "unrelated":
		unrelated_inds(ibdMat=args.ibdMat, pihat=args.maxPIHAT)
	elif args.method.lower() == "findFounder":
		nearestFounder(pedigree=args.ped)
	else:
		try:
			raise SystemExit
		finally:
			print("Method: {}, does not exist!".format(str(args.method.lower())))

