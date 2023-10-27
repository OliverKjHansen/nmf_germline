import pickle
import os

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
signatures = config["NumberOfSignatures"]
pattern_type = config["pattern_type"]
# methylation_data = config["methylation_data"]
# replication_time = config["replication_time"]
# size_partition = config["size_partition"]
# complex_structure = config["complex_structure"]
# recombination = config["recombination"]

#this makes the wildcards of the regions i want to investigate
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
	regions = making_windows(chrom, length, window_sizes)
	with open("regions.pkl", "wb") as file:
		pickle.dump(regions, file)

def creating_breakpoints(kmer):
	before = int(kmer)/2
	after = before-1
	return [before, after]

# def making_windows2(chromosome, length, window_size):
# 	regions = []
# 	size = window_size*1000000
# 	for pos in range(size, int(length[chrom]), size):
# 		title = str(int((pos-size)/1000000))+"mb_to_"+str(int(pos/1000000))+"mb"
# 			regions.append(title)
# 	return regions
rule all:
	input:
		expand(["files/{datasets}/coverage_files/{chrom}.BRAVO_TOPMed_coverage_hg38.txt.gz",
		"files/{datasets}/derived_files/accepted_coverage/{chrom}x10_{fraction}p.bed", # This file contains a bedfile of all the regions that passes the restriction i have put on 80% of the individuals needs to have a coverage of 10x
		"files/{datasets}/derived_files/accepted_coverage/all_coverage_x10_{fraction}p.bed",
		"files/{datasets}/vcf_files/{chrom}.BRAVO_TOPMed_Freeze_8.vcf.gz",
		"files/{datasets}/derived_files/vcf_indels/{chrom}_indel_{freq}.vcf.gz",
		"files/{datasets}/derived_files/vcf_indels/all_indels_{freq}.vcf.gz",
		"{window_sizes}mb_windows/regions/{region}.bed",
		"{window_sizes}mb_windows/filtered_regions/{region}_{fraction}p.bed",
		"{window_sizes}mb_windows/background_{kmer}mer/background_{region}_{kmer}mer_{fraction}p.bed",
		"{window_sizes}mb_windows/variants/indels_{region}_{freq}_{fraction}p.bed", # can be removed # no it shouldnt 
		"{window_sizes}mb_windows/variants/ins_{region}_{freq}_{fraction}p.bed",
		"{window_sizes}mb_windows/variants/del_{region}_{freq}_{fraction}p.bed",
		"{window_sizes}mb_windows/indels_{kmer}mer/frequency_{freq}_at_{fraction}p/ins_counts_{region}_{kmer}mer.bed",
		"{window_sizes}mb_windows/indels_{kmer}mer/frequency_{freq}_at_{fraction}p/del_counts_{region}_{kmer}mer.bed",
		"{window_sizes}mb_windows/indels_{kmer}mer/combined/frequency_{freq}_at_{fraction}p/ins_counts_{kmer}mer.bed",
		"{window_sizes}mb_windows/indels_{kmer}mer/combined/frequency_{freq}_at_{fraction}p/del_counts_{kmer}mer.bed", 
		"{window_sizes}mb_windows/background_{kmer}mer/combined/background_{kmer}mer_{fraction}p.bed",
		"{window_sizes}mb_windows/indels_{kmer}mer/combined/frequency_{freq}_at_{fraction}p/insertions_dataframe_{kmer}mer.rds",
		"{window_sizes}mb_windows/indels_{kmer}mer/combined/frequency_{freq}_at_{fraction}p/deletions_dataframe_{kmer}mer.rds",
		"{window_sizes}mb_windows/models/frequency_{freq}_at_{fraction}p/{types}_{kmer}mer/{types}_{kmer}mer_{signatures}.rds"
		], datasets = datasets, chrom = chrom, fraction = NumberWithDepth, freq = allelefrequency,
		region = regions, window_sizes = window_sizes, kmer = kmer_indels, types = pattern_type, signatures = signatures) #region = regions, window_sizes = window_sizes, kmer = kmer_indels,chrom = chrom, fraction = NumberWithDepth, freq = allelefrequency, size_partition = size_partition, complex_structure = complex_structure)

rule coverage_regions:
	input:
		seq_zipped = "files/{datasets}/coverage_files/{chrom}.BRAVO_TOPMed_coverage_hg38.txt.gz"
	params: 
		procent = lambda wildcards: float(int(wildcards.fraction)/100)
	resources:
		threads=4,
		time=120,
		mem_mb=5000
	output:
		bedfile = "files/{datasets}/derived_files/accepted_coverage/{chrom}x10_{fraction}p.bed" # This file contains a bedfile of all the regions that passes the restriction i have put on 80% of the individuals needs to have a coverage of 10x
	shell:"""
	temp_unzipped=$(mktemp -u)
    gunzip -c {input.seq_zipped} > $temp_unzipped
	python scripts/countingregions.py $temp_unzipped {params.procent} > {output.bedfile}
	gzip $temp_unzipped
	"""
rule vcf_indel:
	input:
		raw_vcf = "files/{datasets}/vcf_files/{chrom}.BRAVO_TOPMed_Freeze_8.vcf.gz"
	conda: "envs/bcftools.yaml"
	output:
		filtered = "files/{datasets}/derived_files/vcf_indels/{chrom}_indel_{freq}.vcf.gz",
	shell:"""
	tabix -f -p vcf {input.raw_vcf}
	bcftools filter -O z -o {output.filtered} -i 'AF<{wildcards.freq} && VRT=2' {input.raw_vcf}
	"""

rule aggregate_chromosomes:
	input:
		individual_coverage = expand("files/{datasets}/derived_files/accepted_coverage/{chrom}x10_{fraction}p.bed", datasets=datasets, chrom=chrom, fraction=NumberWithDepth),
		individual_vcf =  expand("files/{datasets}/derived_files/vcf_indels/{chrom}_indel_{freq}.vcf.gz", datasets=datasets, chrom=chrom, freq = allelefrequency)
	output:
		summary_coverage = expand("files/{datasets}/derived_files/accepted_coverage/all_coverage_x10_{fraction}p.bed", datasets=datasets, fraction=NumberWithDepth),
		summary_vcf= expand("files/{datasets}/derived_files/vcf_indels/all_indels_{freq}.vcf.gz", datasets=datasets, freq = allelefrequency)
	shell:"""
	cat {input.individual_coverage} >> {output.summary_coverage}
	cat {input.individual_vcf} >> {output.summary_vcf}
	"""

# a rule which makes the MegaBases bedfile 
# Creating regions which are to be investigated
rule mega_bases:
	params: 
		chrom = lambda wildcards: wildcards.region.split("_")[0].split("m")[0],
		start = lambda wildcards: str(int(wildcards.region.split("_")[1].split("m")[0])*1000000),
		end = lambda wildcards: str(int(wildcards.region.split("_")[3].split("m")[0])*1000000)
	resources:
		threads=1,
		time=1,
		mem_mb=100
	output:
		bedfiles = "{window_sizes}mb_windows/regions/{region}.bed"
	shell:"""
	printf '%s\t%s\t%s\n' {params.chrom} {params.start} {params.end} > {output.bedfiles}
	"""
# the fitlers for regions are blacklist(add ref), average coverage, and i want to add exome as well
# last couple of lines is accepting the regins, in which more then 50% of the bases are not filtered away
rule filtering_regions:
	input:
		regions = "{window_sizes}mb_windows/regions/{region}.bed",
		coverage_accepted = expand("files/{datasets}/derived_files/accepted_coverage/all_coverage_x10_{fraction}p.bed", datasets=datasets, fraction=NumberWithDepth),
		blacklist = blacklist,
		exons = exons
	conda: "envs/bedtools.yaml"
	resources:
		threads=1,
		time=60,
		mem_mb=5000
	output:
		tmp_cov = temporary("{window_sizes}mb_windows/tmp/tmp_coverage_{region}_{fraction}p.bed"),
		tmp_blacklist = temporary("{window_sizes}mb_windows/tmp/blacklist_{region}_{fraction}p.bed"),
		tmp_exons = temporary("{window_sizes}mb_windows/tmp/exons_{region}_{fraction}p.bed"),
		filtered_regions = "{window_sizes}mb_windows/filtered_regions/{region}_{fraction}p.bed"
	shell:"""
	bedtools intersect -a {input.regions} -b {input.coverage_accepted} > {output.tmp_cov} #make this temp
	bedtools intersect -v -a {output.tmp_cov} -b {input.blacklist} > {output.tmp_blacklist}
	bedtools intersect -v -a {output.tmp_blacklist} -b {input.exons} > {output.tmp_exons}
	tmp=`bedtools intersect -wo -a {input.regions} -b {output.tmp_exons}| awk '{{s+=$7}} END {{print s}}'`
	num=$(expr {window_sizes} \* 1000000 / 2)
	if [[ $tmp -ge $num ]]
	then 
		cp {output.tmp_exons} {output.filtered_regions}
	else
		touch {output.filtered_regions}
	fi
	"""

rule indel_background_counter: #im not sure this works
	input:
		filtered_regions = "{window_sizes}mb_windows/filtered_regions/{region}_{fraction}p.bed",
		genome = two_bit
	conda: "envs/kmer_counter.yaml"
	params:
		before_break  = lambda wildcards: int(creating_breakpoints(wildcards.kmer)[0]),
		after_break = lambda wildcards: int(creating_breakpoints(wildcards.kmer)[1])
	output:
		background = temporary("{window_sizes}mb_windows/tmp/background_{region}_{kmer}mer_{fraction}p.bed"),
		ss_background = "{window_sizes}mb_windows/background_{kmer}mer/background_{region}_{kmer}mer_{fraction}p.bed"
	shell:"""
	check=`cat {input.filtered_regions} | wc -l`
	if [[ $check -gt 0 ]]
	then 
		kmer_counter background --bed {input.filtered_regions} --before_after {params.before_break} {params.after_break} --reverse_complement_method both {input.genome} > {output.background}
		awk -v OFS='\t' '{{print "{wildcards.region}",$1,$2}}' {output.background} > {output.ss_background}
	else
		touch {output.background}
		touch {output.ss_background}
	fi
	"""

rule creating_indel_variants:
	input:
		filtered_regions = "{window_sizes}mb_windows/filtered_regions/{region}_{fraction}p.bed",
		vcf_file = expand("files/{datasets}/derived_files/vcf_indels/all_indels_{freq}.vcf.gz", datasets=datasets, freq = allelefrequency)
	conda: "envs/bedtools.yaml"
	output:
		variants =  "{window_sizes}mb_windows/variants/indels_{region}_{freq}_{fraction}p.bed",
		ins_variants = "{window_sizes}mb_windows/variants/ins_{region}_{freq}_{fraction}p.bed",
		del_variants = "{window_sizes}mb_windows/variants/del_{region}_{freq}_{fraction}p.bed"
	shell:"""
	check=`cat {input.filtered_regions} | wc -l`
	if [[ $check -gt 0 ]]
	then
		bedtools intersect -a {input.vcf_file} -b {input.filtered_regions} | awk -v OFS='\t' '{{print $1,$2,$4,$5}}' > {output.variants}
	else
		touch {output.variants}
	fi
	python scripts/creating_deletions.py {output.variants} > {output.del_variants}
	python scripts/creating_insertions.py {output.variants} > {output.ins_variants}
	"""

rule indel_variant_counter:
	input:
		ins_variants = "{window_sizes}mb_windows/variants/ins_{region}_{freq}_{fraction}p.bed",
		del_variants = "{window_sizes}mb_windows/variants/del_{region}_{freq}_{fraction}p.bed",
		genome = two_bit
	conda: "envs/kmer_counter.yaml"
	params:
		radius  = lambda wildcards: int(int(wildcards.kmer)/2)
	output:
		kmer_count_ins = temporary("{window_sizes}mb_windows/tmp/ndels_{kmer}mer/frequency_{freq}_at_{fraction}p/ins_counts_{region}_{kmer}mer.bed"),
		kmer_count_del = temporary("{window_sizes}mb_windows/tmp/indels_{kmer}mer/frequency_{freq}_at_{fraction}p/del_counts_{region}_{kmer}mer.bed"),
		ss_ins = "{window_sizes}mb_windows/indels_{kmer}mer/frequency_{freq}_at_{fraction}p/ins_counts_{region}_{kmer}mer.bed",
		ss_del = "{window_sizes}mb_windows/indels_{kmer}mer/frequency_{freq}_at_{fraction}p/del_counts_{region}_{kmer}mer.bed"
	shell:"""
	check=`cat {input.ins_variants} | wc -l`
	if [[ $check -gt 0 ]]
	then
		kmer_counter indel -r {params.radius} --sample {input.genome} {input.ins_variants} ins > {output.kmer_count_ins}
		kmer_counter indel -r {params.radius} --sample {input.genome} {input.del_variants} del_start > {output.kmer_count_del}
		awk -v OFS='\t' '{{print "{wildcards.region}",$1,$2}}' {output.kmer_count_ins} > {output.ss_ins} 
		awk -v OFS='\t' '{{print "{wildcards.region}",$1,$2}}' {output.kmer_count_del} > {output.ss_del}
	else
		touch {output.kmer_count_del}
		touch {output.kmer_count_ins}
		touch {output.ss_del}
		touch {output.ss_ins} 
	fi
	"""
##Make a check for the directories
rule aggregate_indels_regions:
	input:
		insertions = expand("{window_sizes}mb_windows/indels_{{kmer}}mer/frequency_{freq}_at_{fraction}p/ins_counts_{region}_{{kmer}}mer.bed", fraction = NumberWithDepth, freq = allelefrequency, region = regions, window_sizes = window_sizes, kmer = kmer_indels),
		deletions = expand("{window_sizes}mb_windows/indels_{{kmer}}mer/frequency_{freq}_at_{fraction}p/del_counts_{region}_{{kmer}}mer.bed", fraction = NumberWithDepth, freq = allelefrequency, region = regions, window_sizes = window_sizes, kmer = kmer_indels),
		background = expand("{window_sizes}mb_windows/background_{{kmer}}mer/background_{region}_{{kmer}}mer_{fraction}p.bed", fraction = NumberWithDepth, freq = allelefrequency, region = regions, window_sizes = window_sizes, kmer = kmer_indels)
	output:
		summary_insertions = expand("{window_sizes}mb_windows/indels_{{kmer}}mer/combined/frequency_{freq}_at_{fraction}p/ins_counts_{{kmer}}mer.bed", fraction = NumberWithDepth, freq = allelefrequency, window_sizes = window_sizes, kmer = kmer_indels),
		summary_deletions = expand("{window_sizes}mb_windows/indels_{{kmer}}mer/combined/frequency_{freq}_at_{fraction}p/del_counts_{{kmer}}mer.bed", fraction = NumberWithDepth, freq = allelefrequency, window_sizes = window_sizes, kmer = kmer_indels),
		summary_background = expand("{window_sizes}mb_windows/background_{{kmer}}mer/combined/background_{{kmer}}mer_{fraction}p.bed", fraction = NumberWithDepth, freq = allelefrequency, window_sizes = window_sizes, kmer = kmer_indels)
	params:
	shell:"""
	cat {input.insertions} >> {output.summary_insertions}
	cat {input.deletions} >> {output.summary_deletions}
	cat {input.background} >> {output.summary_background}
	"""
###Now let do some nmf###

rule prepare_for_nmf:
	input:
		summary_insertions = "{window_sizes}mb_windows/indels_{kmer}mer/combined/frequency_{freq}_at_{fraction}p/ins_counts_{kmer}mer.bed",
		summary_deletions = "{window_sizes}mb_windows/indels_{kmer}mer/combined/frequency_{freq}_at_{fraction}p/del_counts_{kmer}mer.bed",
		summary_background = "{window_sizes}mb_windows/background_{kmer}mer/combined/background_{kmer}mer_{fraction}p.bed"
	conda: "envs/callr.yaml"
	output:
		insertions_dataframe = "{window_sizes}mb_windows/indels_{kmer}mer/combined/frequency_{freq}_at_{fraction}p/insertions_dataframe_{kmer}mer.rds",
		deletions_dataframe = "{window_sizes}mb_windows/indels_{kmer}mer/combined/frequency_{freq}_at_{fraction}p/deletions_dataframe_{kmer}mer.rds"
	shell:"""
	Rscript scripts/creating_dataframes.R {input.summary_background} {input.summary_insertions} {input.summary_deletions} {output.deletions_dataframe} {output.insertions_dataframe}
	"""
rule modelselection:
	input:
		count_data = "{window_sizes}mb_windows/indels_{kmer}mer/combined/frequency_{freq}_at_{fraction}p/{types}_dataframe_{kmer}mer.rds"
	conda: "envs/nmf.yaml"
	resources:
		threads=2,
		time=480,
		mem_mb=4000
	output:
		model = "{window_sizes}mb_windows/models/frequency_{freq}_at_{fraction}p/{types}_{kmer}mer/{types}_{kmer}mer_{signatures}.rds"
	shell:"""
    Rscript scripts/opportunity_modelselection.R {wildcards.signatures} {input.count_data} {output.model}
    """
### Types is not implemented across all wildcards
sss
rule plotting:
	input:
	conda: 
	resources:
	output:
	shell:"""
    
    """

# ### DOING METHYLATION DATA for plot

# rule methylation_slicing: #im not sure this works
# 	input:
# 		accept_regions = "{window_sizes}mb_windows/accepted_regions/accepted_{region}_{fraction}p.bed",
# 		methylation = methylation_data
# 	resources:
# 		threads=1,
# 		time=60,
# 		mem_mb=1000
# 	output:
# 		#intersected_meth = "{window_sizes}mb_windows/methylation/intersected",
# 		ss_meth = "{window_sizes}mb_windows/methylation_{fraction}p/{region}_methylation.bed"
# 	shell:"""
# 	check=`cat {input.accept_regions}  | wc -l`
# 	if [[ $check -gt 0 ]]
# 	then 
# 		bedtools intersect -wa -a {input.methylation} -b {input.accept_regions} | awk -v OFS='\t' '{{print "{wildcards.region}",$2,$3,$5,$10,$11,$12,$13,$14}}' > {output.ss_meth}
# 	else
# 		printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' {wildcards.region} NA NA NA NA NA NA NA NA > {output.ss_meth}
# 	fi
# 	"""

# ##replication_time
# rule replication_timing: #im not sure this works
# 	input:
# 		accept_regions = "{window_sizes}mb_windows/accepted_regions/accepted_{region}_{fraction}p.bed",
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

# # DOING SIZE VARIANTS

# rule size_indel_variants:
# 	input:
# 		ins_variants = "{window_sizes}mb_windows/variants/ins_{region}_{freq}_{fraction}p.bed",
# 		del_variants = "{window_sizes}mb_windows/variants/del_{region}_{freq}_{fraction}p.bed"
# 	params:
# 		#size = lambda wildcards: ",".join(wildcards.size_partition.split(","))
# 	output:
# 		insertion = "{window_sizes}mb_windows/size_difference_{freq}_{fraction}p/variant/ins_{region}_{size_partition}.bed",
# 		deletion = "{window_sizes}mb_windows/size_difference_{freq}_{fraction}p/variant/del_{region}_{size_partition}.bed"
# 	shell:"""
# 	check=`cat {input.ins_variants} | wc -l`
# 	if [[ $check -gt 0 ]]
# 	then
# 		python splitting_indel_size.py {input.ins_variants} {wildcards.size_partition} > {output.insertion} 
# 		python splitting_indel_size.py {input.del_variants} {wildcards.size_partition} > {output.deletion}
		
# 	else
# 		touch {output.insertion}
# 		touch {output.deletion}
# 	fi
# 	"""

# rule indel_size_counter:
# 	input:
# 		insertion_size = "{window_sizes}mb_windows/size_difference_{freq}_{fraction}p/variant/ins_{region}_{size_partition}.bed",
# 		deletion_size = "{window_sizes}mb_windows/size_difference_{freq}_{fraction}p/variant/del_{region}_{size_partition}.bed",
# 		genome = two_bit
# 	params:
# 		radius  = lambda wildcards: int(int(wildcards.kmer)/2)
# 	output:
# 		insertion_sized = "{window_sizes}mb_windows/size_difference_{freq}_{fraction}p/tmp/ins_{region}_{size_partition}_{kmer}mer_tmp.bed",
# 		deletion_sized = "{window_sizes}mb_windows/size_difference_{freq}_{fraction}p/tmp/del_{region}_{size_partition}_{kmer}mer_tmp.bed",
# 		insertion_size_ss = "{window_sizes}mb_windows/size_difference_{freq}_{fraction}p/final/ins_{region}_{size_partition}_{kmer}mer_final.bed",
# 		deletion_size_ss = "{window_sizes}mb_windows/size_difference_{freq}_{fraction}p/final/del_{region}_{size_partition}_{kmer}mer_final.bed"
# 	shell:"""
# 	check=`cat {input.insertion_size} | wc -l`
# 	if [[ $check -gt 0 ]]
# 	then
# 		kmer_counter indel -r {params.radius} --sample {input.genome} {input.insertion_size} ins > {output.insertion_sized}
# 		kmer_counter indel -r {params.radius} --sample {input.genome} {input.deletion_size} del_start > {output.deletion_sized}
# 		awk -v OFS='\t' '{{print "{wildcards.region}","{wildcards.size_partition}",$1,$2}}' {output.insertion_sized} > {output.insertion_size_ss} 
# 		awk -v OFS='\t' '{{print "{wildcards.region}","{wildcards.size_partition}",$1,$2}}' {output.deletion_sized} > {output.deletion_size_ss}

# 	else
# 		touch {output.insertion_sized}
# 		touch {output.deletion_sized}
# 		touch {output.deletion_size_ss}
# 		touch {output.insertion_size_ss}
# 	fi
# 	"""

# ####COMPLEX STRUCTURES#####

# rule creating_complex_structure:
# 	input:
# 		accepted_regions = "{window_sizes}mb_windows/accepted_regions/accepted_{region}_{fraction}p.bed",
# 		complex_structure = "nonBdna/{complex_structure}.bed" #change, make one or two rules that create this file(filter and ) output?
# 	output:
# 		cs =  "{window_sizes}mb_windows/complex_structures/{region}_{complex_structure}_{fraction}p.bed"
# 	shell:"""
# 	check=`cat {input.accepted_regions} | wc -l`
# 	if [[ $check -gt 0 ]]
# 	then
# 		bedtools intersect -a {input.complex_structure} -b {input.accepted_regions} | awk -v OFS='\t' '{{print $0}}' > {output.cs}
# 	else
# 		touch {output.cs}
# 	fi
# 	"""

# ####Recombination rate######

# rule recombination_rate:
# 	input:
# 		accepted_regions = "{window_sizes}mb_windows/accepted_regions/accepted_{region}_{fraction}p.bed",
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