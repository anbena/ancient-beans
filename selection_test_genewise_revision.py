#! /usr/bin/env python

import sys
from operator import itemgetter 

out = open('selscan_hint_null', 'w')

out.write('gene\tchr\tstart\tend\tlength\tinvariants\tmissingdata\tsnp_total\ttopoB\ttopoB_pos\ttopoA\ttopoA_pos\ttopoC\ttopoC_pos\ttopoE\ttopoH\ttopoX\n')

snp_total = 0
snp_modern = 0


snp_ancient1modern = 0


snp_ancient_all = 0


topoE = 0 # WA<-->CHI,Aa,Am,OUT
topoH = 0 # OUT<-->CHI,Aa,Am,WA
topoX = 0 # variability intra A

input_file = 'shuffled.tsv'
#input_file = 'snp.genes1kb.tsv'

excluded = 0
invariant_all = 0
invariant_in = 0
positions_DR = []
positions_DA = []
positions_DAA = []

gene = 'null'
for line in open(input_file, 'r').readlines():
	line = line.split('\n')[0].split('\t')
	if line == '' or line == '\n':
		break
	if gene == 'null' or gene == line[8]:
		gene = line[8]
		start = line[3]
		end = line[4]
		chrom = line[0]
		length = int(line[4])-int(line[3])
		pos = line[10]

		#get the genotypes in each group of samples while removing missing data coded as 'N'
		alleles_ancient1 = [i for i in itemgetter(*[12,13,15,19,21,24,28,29,30])(line) if i != 'N'] #ancient1: ATR e ATM 1500-500 BP
		alleles_ancient2 = [i for i in itemgetter(*[14,16,17,18,23,25])(line) if i != 'N'] #ancient2: ANS 3500-1500 BP
		#alleles_ancient3 = [i for i in itemgetter(*[26])(line) if i != 'N'] #ancient3: SJ28 6000 BP to be confirmed
		alleles_ancient = alleles_ancient1+alleles_ancient2
		
		alleles_dom_and = [i for i in line[34:36] if i != 'N'] #chilean
		alleles_outgroup = [i for i in line[33] if i != 'N'] #Hintoni
		alleles_wild_and = [i for i in line[31] if i != 'N'] #Wild from Andes
		#alleles_outgroup = [i for i in line[32] if i != 'N'] #Wild from Chiapas
		alleles_dom_meso = [i for i in line[36:len(line)] if i != 'N'] #Domesticated from the two major mesoamerican groups

		if len(set([i for i in line[12:len(line)] if i != 'N'])) == 1:
			invariant_all +=1
		elif len(set([i for i in alleles_dom_and+alleles_ancient+alleles_wild_and+alleles_outgroup])) == 1:
			invariant_in +=1
		elif len(alleles_ancient) == 0 or len(alleles_dom_and) == 0 or len(alleles_outgroup) == 0 or len(alleles_wild_and) == 0:
			excluded += 1
		else:
			snp_total += 1
			
			if len([i for i in set(alleles_dom_and) if i in alleles_ancient+alleles_wild_and+alleles_outgroup]) == 0:
				snp_modern += 1
				positions_DR.append(pos)


			if len([i for i in set(alleles_dom_and+alleles_ancient) if i in alleles_wild_and+alleles_outgroup]) == 0:
				snp_ancient1modern += 1
				positions_DA.append(pos)


			if len([i for i in set(alleles_ancient) if i in alleles_dom_and+alleles_wild_and+alleles_outgroup]) == 0:
				snp_ancient_all += 1
				positions_DAA.append(pos)


			if len([i for i in set(alleles_ancient+alleles_dom_and+alleles_outgroup) if i in alleles_wild_and]) == 0:
				topoE += 1
			if len([i for i in set(alleles_ancient+alleles_dom_and+alleles_wild_and) if i in alleles_outgroup]) == 0:
				topoH += 1
			if len(set(alleles_ancient+alleles_dom_and)) > 1 and len([i for i in set(alleles_ancient+alleles_dom_and) if i in alleles_outgroup+alleles_wild_and]) > 0:
				topoX += 1
			#else:
			#	print(alleles_ancient, alleles_dom_and, alleles_outgroup, alleles_wild_and)

	else:
		out.write(gene+'\t'+chrom+'\t'+start+'\t'+end+'\t'+str(length)+'\t'+str(invariant_in)+'\t'+str(excluded)+'\t'+str(snp_total)+'\t'+str(snp_modern)+'\t'+(',').join(positions_DR)+'\t'+str(snp_ancient1modern)+'\t'+(',').join(positions_DA)+'\t'+str(snp_ancient_all)+'\t'+(',').join(positions_DAA)+'\t'+str(topoE)+'\t'+str(topoH)+'\t'+str(topoX)+'\n')
		
		snp_total = 0
		snp_modern = 0

		snp_ancient1modern = 0

		snp_ancient_all = 0


		topoE = 0 # WA<-->CHI,Aa,Am,OUT
		topoH = 0 # OUT<-->CHI,Aa,Am,WA
		topoX = 0 # variability intra A

		excluded = 0
		invariant_all = 0
		invariant_in = 0
		positions_DR = []
		positions_DA = []
		positions_DAA = []

		gene = line[8]
		start = line[3]
		end = line[4]
		chrom = line[0]
		length = int(line[4])-int(line[3])
		pos = line[10]

		#get the genotypes in each group of samples while removing missing data coded as 'N'
		alleles_ancient1 = [i for i in itemgetter(*[12,13,15,19,21,24,28,29,30])(line) if i != 'N'] #ancient1: ATR e ATM 1500-500 BP
		alleles_ancient2 = [i for i in itemgetter(*[14,16,17,18,23,25])(line) if i != 'N'] #ancient2: ANS 3500-1500 BP
		#alleles_ancient3 = [i for i in itemgetter(*[26])(line) if i != 'N'] #ancient3: SJ28 6000 BP to be confirmed
		alleles_ancient = alleles_ancient1+alleles_ancient2
		
		alleles_dom_and = [i for i in line[34:36] if i != 'N'] #chilean
		#alleles_outgroup = [i for i in line[33] if i != 'N'] #Hintoni
		alleles_wild_and = [i for i in line[31] if i != 'N'] #Wild from Andes
		alleles_outgroup = [i for i in line[32] if i != 'N'] #Wild from Chiapas
		alleles_dom_meso = [i for i in line[36:len(line)] if i != 'N'] #Domesticated from the two major mesoamerican groups

		if len(set([i for i in line[12:len(line)] if i != 'N'])) == 1:
			invariant_all +=1
		elif len(set([i for i in alleles_dom_and+alleles_ancient+alleles_wild_and+alleles_outgroup])) == 1:
			invariant_in +=1
		elif len(alleles_ancient) == 0 or len(alleles_dom_and) == 0 or len(alleles_outgroup) == 0 or len(alleles_wild_and) == 0:
			excluded += 1
		else:
			snp_total += 1
			
			if len([i for i in set(alleles_dom_and) if i in alleles_ancient+alleles_wild_and+alleles_outgroup]) == 0:
				snp_modern += 1
				positions_DR.append(pos)


			if len([i for i in set(alleles_dom_and+alleles_ancient) if i in alleles_wild_and+alleles_outgroup]) == 0:
				snp_ancient1modern += 1
				positions_DA.append(pos)


			if len([i for i in set(alleles_ancient) if i in alleles_dom_and+alleles_wild_and+alleles_outgroup]) == 0:
				snp_ancient_all += 1
				positions_DAA.append(pos)


			if len([i for i in set(alleles_ancient+alleles_dom_and+alleles_outgroup) if i in alleles_wild_and]) == 0:
				topoE += 1
			if len([i for i in set(alleles_ancient+alleles_dom_and+alleles_wild_and) if i in alleles_outgroup]) == 0:
				topoH += 1
			if len(set(alleles_ancient+alleles_dom_and)) > 1 and len([i for i in set(alleles_ancient+alleles_dom_and) if i in alleles_outgroup+alleles_wild_and]) > 0:
				topoX += 1
			#else:
			#	print(alleles_ancient, alleles_dom_and, alleles_outgroup, alleles_wild_and)


out.write(gene+'\t'+chrom+'\t'+start+'\t'+end+'\t'+str(length)+'\t'+str(invariant_in)+'\t'+str(excluded)+'\t'+str(snp_total)+'\t'+str(snp_modern)+'\t'+(',').join(positions_DR)+'\t'+str(snp_ancient1modern)+'\t'+(',').join(positions_DA)+'\t'+str(snp_ancient_all)+'\t'+(',').join(positions_DAA)+'\t'+str(topoE)+'\t'+str(topoH)+'\t'+str(topoX)+'\n')
out.close()

		











			

