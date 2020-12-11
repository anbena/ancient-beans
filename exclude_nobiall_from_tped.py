data = open('haplo.allchr.neutral_filtered.tped', 'r')
out = open('haplo.allchr.neutral_filtered.biall.tped', 'w')
count_tri = 0

for i in data.readlines():
	line = i.strip().split('\t')
	geno = line[4:len(line)]
	haplo = []
	for j in geno:
		haplo.append(j[0])
	haplo_nomiss = [x for x in haplo if x != '0']
	alleles = set(haplo_nomiss)
	if len(alleles) < 3:
		out.write(i)
	else:
		count_tri += 1
data.close()
out.close()

print 'SNPs with more than 2 alleles: ', count_tri

