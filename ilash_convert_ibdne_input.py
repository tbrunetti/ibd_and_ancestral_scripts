import pandas
import numpy as np
import argparse
import sys

# A collection of functions to be used with iLash

'''
func: getROH()
input:
output:
'''
def getROH():
	pass

'''
func: sampleSelection()
input:
output:
'''
def sampleSelection(inputFile, select, output):
	print('Running sampleSelection...\n')


'''
func: iLash2ibdNe()
input:
output:
'''
def iLash2ibdNe(inputFile, output):
	print('Running iLash2ibdNe...\n')


'''
func: formatSubset()
input:
output:
'''
def formatSubset(inputFile, select, output):
	print('Running both sample selection and formatting methods...\n')

	#ilashFile = 'Tangier_cleaned_final_added_ped_duohmm_chr22.match.keep'
	#maxInd = 'max_independent_set_removed_samples_with_9000s.txt'
	#output = 'Tangier_cleaned_final_added_ped_duohmm_chr22.match.keep.ibd'

	ilashResults = pandas.read_table(inputFile, header=None,
		names=['FID1', 'IID1', 'FID2', 'IID2', 'chrm', 'start', 'end', 'snpStart', 'snpEnd', 'depth', 'JC_similarity'],
		dtype={'start':np.int64, 'end':np.int64, 'depth':np.float64, 'JC_similarity':np.float64})

	independentInds = list()
	with open(select, 'r') as getInd:
		for line in getInd:
			independentInds.append(line.strip())


	ilashResults['IID1_name'], ilashResults['IID1_haplotype'] = ilashResults['IID1'].str.split('_', 1).str
	ilashResults['IID2_name'], ilashResults['IID2_haplotype'] = ilashResults['IID2'].str.split('_', 1).str

	try:
		ilashResults['IID1_haplotype'] = ilashResults['IID1_haplotype'].astype(int)
		ilashResults['IID2_haplotype'] = ilashResults['IID2_haplotype'].astype(int)
	except TypeError:
		print("Wrong type for IID1 and/or IID2 haplotypes")

	print(ilashResults)
	print(independentInds)
	ilashResults['firstID_check'] = ilashResults['IID1_name'].isin(independentInds)
	ilashResults['secondID_check'] = ilashResults['IID2_name'].isin(independentInds)

	indSetResults = ilashResults.loc[(ilashResults['firstID_check'] & ilashResults['secondID_check'])]
	indSetResults['fake_LOD'] = 100
	indSetResults[['IID1_name', 'IID1_haplotype', 'IID2_name', 'IID2_haplotype', 'chrm', 'start', 'end', 'fake_LOD']].to_csv('{}.ibd'.format(output), index=False, sep="\t", header=False)


'''
func: get_shared_haplotypes()
input:
output:
'''
def get_shared_haplotypes(inputFile):

	pairwise = pandas.read_table(inputFile, names=['IID1', 'haplo1', 'IID2', 'haplo2', 'chrm', 'startPos', 'endPos', 'IBDseqScore'],
		dtype={'IID1':str, 'haplo1':str, 'IID2':str, 'haplo2':str, 'chrm':str, 'startPos':np.int64, 'endPos':np.int64})

	pairwise['bp_length']=pairwise['endPos']-pairwise['startPos']
	pairwise['region']=pairwise['startPos'].astype(str)+"_"+pairwise['endPos'].astype(str)
	

	sameHaps = pairwise.loc[pairwise['haplo1'] == pairwise['haplo2']]

	pandas.value_counts(sameHaps['region'].values, sort=True)

	sameHaps.to_csv("test.tsv", index=False, sep="\t")





if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Methods for using iLash of IBD')
	
	## TO DO:  add arguments
	parser.add_argument('-inputFile', required=True, dest='inputFile', type=str, help='iLash output (.match file prefix) OR output of format_subset if using method sharedHaps')
	parser.add_argument('-method', required=True, dest='method', type=str, help='format, subset, format_subset, sharedHaps, roh')
	parser.add_argument('--keep', default=None, dest='keep', type=str, help='file path to sample ids want to keep/use in ilash dataset (one sample ID per line)')
	parser.add_argument('--output', default=None, dest='output', type=str, help='prefix to use for output file')
	args = parser.parse_args()
	

	## TO DO: pass arguments and set up logic for sampleSelection() and iLash2ibdNe()
	if args.method.lower() == 'format':
		iLash2ibdNe(inputFile = args.inputFile, output = args.output)
	elif args.method.lower() == 'subset':
		sampleSelection(inputFile = args.inputFile, select = args.keep, output = args.output)
	elif args.method.lower() == 'format_subset':
		formatSubset(inputFile = args.inputFile, select = args.keep, output = args.output)
	elif args.method.lower() == 'sharedHaps':
		get_shared_haplotypes(inputFile = args.inputFile)
	elif args.method.lower() == 'roh':
		getROH()
	else:
		try:
			raise SystemExit
		finally:
			print("Method: {}, does not exist!".format(str(args.method.lowerd())))