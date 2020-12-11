#!/usr/bin/awk -f
#Usage: zcat test.haplo.gz | sed '1d' | ./pairdist.awk 
BEGIN{
	nodam=1;#0:with damage mut; 1: without damage mut
	offs=3#number of col before first indiv
	badmut=0;
}
{
if (NR==1){
        j=0;
        for (i=offs+1;i<=(NF-1);i++){
                for (k=(i+1);k<=NF;k++){
                        j++;
                        if (nodam==0){
				if ((i==(NF-1))&&(k==NF)){
                                	printf "USABLE_%s_%s\tDIFF_%s_%s\n",i-offs,k-offs,i-offs,k-offs;
				}else{
					printf "USABLE_%s_%s\tDIFF_%s_%s\t",i-offs,k-offs,i-offs,k-offs;	
				}
                        }else{
				if ((i==(NF-1))&&(k==NF)){
                                	printf "USABLE_%s_%s\tPOSTM_%s_%s\tDIFF_%s_%s\n",i-offs,k-offs,i-offs,k-offs,i-offs,k-offs;
				}else{
					printf "USABLE_%s_%s\tPOSTM_%s_%s\tDIFF_%s_%s\t",i-offs,k-offs,i-offs,k-offs,i-offs,k-offs;
				}
                        }
                }
        }
        j=0;
        for (i=offs+1;i<=(NF-1);i++){
                for (k=(i+1);k<=NF;k++){
                        j++;
                        usable[j]=0;#the site is not missing in each comparison 
                        diff[j]=0;#nuber of different base
                }
        }
}
#NB: added a check if a dna positon is compatible with a postmortem mutation (C->T;G->A) only in ancient samples (column 4 & 5)
if (nodam==0){
	j=0;
	for (i=offs+1;i<=(NF-1);i++){
		for (k=(i+1);k<=NF;k++){
			j++;
			if (($i!="N")&&($k!="N")){
				usable[j]++;
				if ($i!=$k){
					diff[j]++;
					#print NR;
					#printf "combinazione %s %s\n",i,k;
					#printf "alleli %s %s\n",$i,$k
				}
			}
		}
	}
}else{
	bad=0;
	t1=0;
	t2=0;
	t3=0;
	t4=0;
	for (i=offs+1;i<=19+offs;i++){#controllo solo i campioni antichi 1-19
		if ($i=="A"){
			t1=1;
		}
		if ($i=="G"){
			t2=1;
		}
		if ($i=="C"){
                        t3=1;
                }
		if ($i=="T"){
                        t4=1;
                }
	}
	if (((t1+t2)==2)||((t3+t4)==2)){
		bad=1;
		badmut++
	}
	if (bad==0){
		j=0;
		for (i=offs+1;i<=(NF-1);i++){
			for (k=(i+1);k<=NF;k++){
				j++;
				if (($i!="N")&&($k!="N")){
					usable[j]++;
					if ($i!=$k){diff[j]++;}
				}
			}
		}
	}
}
}
END{
	j=0;
	for (i=offs+1;i<=(NF-1);i++){
		for (k=(i+1);k<=NF;k++){
			j++;
			if (nodam==0){
				if ((i==(NF-1))&&(k==NF)){
					printf "%s\t%s\n",usable[j],diff[j];
				}else{
					printf "%s\t%s\t",usable[j],diff[j];
				}
			}else{
				if ((i==(NF-1))&&(k==NF)){
					printf "%s\t%s\t%s\n",usable[j],badmut,diff[j];
				}else{
					printf "%s\t%s\t%s\t",usable[j],badmut,diff[j];
				}
			}
		}
	}
}
