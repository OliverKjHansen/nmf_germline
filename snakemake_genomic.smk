import pickle
import os
os.environ['OPENBLAS_NUM_THREADS'] = '10'
import glob

configfile: 'config.yaml'

###IMPORTENT###
#file_path = config["file_path"]
# all the 
chrom = config["chromosomes"]
length = config["chromosome_length"]
datasets = config["datasets"]
blacklist = config["blacklist"]
exons = config["exons"]

#genome_bed = config["genome_bed"]
#topmed = config["topmed_data"] # when i want to start filtering 
#coverage_files = config["coverage"] # use this to make bedfiles, maybe later
window_sizes = config["window_size_mb"]
kmer_indels = config["kmer_indels"]
NumberWithDepth = config["NumberWithDepth"]
two_bit = config["twobitgenome"]
allelefrequency = config["allelefrequency"]
# signatures = config["NumberOfSignatures"]
# pattern_type = config["pattern_type"]
# size_partition = config["size_partition"]
methylation_data = config["methylation_data"]
replication_time = config["replication_time"]
size_partition = config["size_partition"]
complex_structure = config["complex_structure"]
recombination = config["recombination"]

def making_windows(chromosomes, length, window_sizes):
	regions = []
	for window_size in window_sizes:
		size = window_size*1000000 # for an indel analysis only windows with the size of 10^6 is doable
		for chrom in chromosomes:
			for pos in range(size, int(length[chrom]), size):
				title = chrom+"_"+str(int((pos-size)/1000000))+"mb_to_"+str(int(pos/1000000))+"mb"
				regions.append(title)
	return regions

# Check if the file exists
if os.path.exists("regions.pkl"):
    # Load the list from the file
    with open("regions.pkl", "rb") as file:
        regions = pickle.load(file)
else:
	os.makedirs("1mb_windows/regions/")
	regions = making_windows(chrom, length, window_sizes)
	for region in regions:
		print(region)
		chrom = region.split("_")[0].split("m")[0]
		start = str(int(region.split("_")[1].split("m")[0])*1000000)
		end = str(int(region.split("_")[3].split("m")[0])*1000000)
		filename = f"1mb_windows/regions/{region}.bed"
		with open(filename, "w") as file:
			file.write(f"{chrom}\t{start}\t{end}\n")
	with open("regions.pkl", "wb") as file:
		pickle.dump(regions, file)
rule all: 
	input:
		expand(["{window_sizes}mb_windows/filtered_regions_{fraction}p/{region}.bed",
		"{window_sizes}mb_windows/methylation_{fraction}p/{region}_methylation.bed",
		"{window_sizes}mb_windows/replication_timing_{fraction}p/{region}_replicationtime.bed",
		"{window_sizes}mb_windows/recombination/{region}_recombination_{fraction}p.bed",
		"{window_sizes}mb_windows/complex_structures/{region}_{complex_structure}_{fraction}p.bed"
		], window_sizes = window_sizes, fraction = NumberWithDepth, region = regions, complex_structure = complex_structure)

### DOING METHYLATION DATA for plot
rule genomic_features: #im not sure this works
	input:
		accept_regions = "{window_sizes}mb_windows/filtered_regions_{fraction}p/{region}.bed",
		methylation = methylation_data,
		rep_time = replication_time,
		recombination_map = recombination #change, make one or two rules that create this file(filter and ) output?
	conda: "envs/bedtools.yaml"
	resources:
		threads=1,
		time=60,
		mem_mb=5000 # should work with 1000mb
	output:
		#intersected_meth = "{window_sizes}mb_windows/methylation/intersected",
		ss_meth = "{window_sizes}mb_windows/methylation_{fraction}p/{region}_methylation.bed",
		reptime = "{window_sizes}mb_windows/replication_timing_{fraction}p/{region}_replicationtime.bed",
		recomb =  "{window_sizes}mb_windows/recombination/{region}_recombination_{fraction}p.bed"
	shell:"""
	check=`cat {input.accept_regions}  | wc -l`
	if [[ $check -gt 0 ]]
	then 
		bedtools intersect -wa -a {input.methylation} -b {input.accept_regions} | awk -v OFS='\t' '{{print "{wildcards.region}",$2,$3,$5,$10,$11,$12,$13,$14}}' > {output.ss_meth}
		bedtools intersect -wa -a {input.rep_time} -b {input.accept_regions} | awk -v OFS='\t' '{{print "{wildcards.region}",$0}}' > {output.reptime}
		bedtools intersect -a {input.recombination_map} -b {input.accept_regions} | awk -v OFS='\t' '{{print "{wildcards.region}",$0}}' > {output.recomb}

	else
		printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' {wildcards.region} NA NA NA NA NA NA NA NA > {output.ss_meth}
		printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' {wildcards.region} NA NA NA NA NA > {output.reptime}
		touch {output.recomb}
	fi
	"""

##replication_time
# rule replication_timing: #im not sure this works
# 	input:
# 		accept_regions = "{window_sizes}mb_windows/filtered_regions_{fraction}p/{region}.bed",
# 		rep_time = replication_time
# 	resources:
# 		threads=1,
# 		time=60,
# 		mem_mb=1000
# 	output:
# 		#intersected_meth = "{window_sizes}mb_windows/methylation/intersected",
# 		reptime = "{window_sizes}mb_windows/replication_timing_{fraction}p/{region}_replicationtime.bed"
# 	shell:"""
# 	check=`cat {input.accept_regions}  | wc -l`
# 	if [[ $check -gt 0 ]]
# 	then 
# 		bedtools intersect -wa -a {input.rep_time} -b {input.accept_regions} | awk -v OFS='\t' '{{print "{wildcards.region}",$0}}' > {output.reptime}
# 	else
# 		printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' {wildcards.region} NA NA NA NA NA > {output.reptime}
# 	fi
# 	"""

####Recombination rate######
# rule recombination_rate:
# 	input:
# 		accepted_regions = "{window_sizes}mb_windows/filtered_regions_{fraction}p/{region}.bed",
# 		recombination_map = recombination #change, make one or two rules that create this file(filter and ) output?
# 	output:
# 		recomb =  "{window_sizes}mb_windows/recombination/{region}_recombination_{fraction}p.bed"
# 	shell:"""
# 	check=`cat {input.accepted_regions} | wc -l`
# 	if [[ $check -gt 0 ]]
# 	then
# 		bedtools intersect -a {input.recombination_map} -b {input.accepted_regions} | awk -v OFS='\t' '{{print "{wildcards.region}",$0}}' > {output.recomb}
# 	else
# 		touch {output.recomb}
# 	fi
# 	"""
####COMPLEX STRUCTURES#####
rule creating_complex_structure:
	input:
		accepted_regions = "{window_sizes}mb_windows/filtered_regions_{fraction}p/{region}.bed",
		complex_structure = "files/nonBdna/{complex_structure}.bed" #change, make one or two rules that create this file(filter and ) output?
	conda: "envs/bedtools.yaml"
	output:
		cs =  "{window_sizes}mb_windows/complex_structures/{region}_{complex_structure}_{fraction}p.bed"
	shell:"""
	check=`cat {input.accepted_regions} | wc -l`
	if [[ $check -gt 0 ]]
	then
		bedtools intersect -a {input.complex_structure} -b {input.accepted_regions} | awk -v OFS='\t' '{{print $0}}' > {output.cs}
	else
		touch {output.cs}
	fi
	"""
