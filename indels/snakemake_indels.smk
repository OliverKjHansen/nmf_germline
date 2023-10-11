# snakemake --profile slurm
#sort -k1,1 -k2,2n
configfile: 'config.yaml'

chrom = config["chromosomes"]
length = config["chromosome_length"]

genome_bed = config["genome_bed"]
blacklist = config["blacklist"]
topmed = config["topmed_data"] # when i want to start filtering 
#coverage_files = config["coverage"] # use this to make bedfiles, maybe later
window_sizes = config["window_size_mb"]
kmer_indels = config["kmer_indels"]
NumberWithDepth = config["NumberWithDepth"]
two_bit = config["twobitgenome"]
allelefrequency = config["allelefrequency"]
methylation_data = config["methylation_data"]
replication_time = config["replication_time"]
size_partition = config["size_partition"]
complex_structure = config["complex_structure"]
recombination = config["recombination"]

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

regions = making_windows(chrom, length, window_sizes)

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
		expand(["coverage_bedfiles/{chrom}x10_{fraction}p.bed",
		"coverage_bedfiles/all_chromosomesx10_{fraction}p.bed",
		"raw_vcf/{chrom}_indel_{freq}.vcf.gz",
		"raw_vcf/indel_{freq}.vcf.gz",
		"{window_sizes}mb_windows/clean/{region}.bed",
		"{window_sizes}mb_windows/intersected/intersected_{region}_{fraction}p.bed", # remove
		"{window_sizes}mb_windows/final/final_{region}_{fraction}p.bed", # remove
		"{window_sizes}mb_windows/accepted_regions/accepted_{region}_{fraction}p.bed", # remove
		"{window_sizes}mb_windows/background_{kmer}mer/background_{region}_{kmer}mer_{fraction}p.bed", # remove
		"{window_sizes}mb_windows/variants/indels_{region}_{freq}_{fraction}p.bed", # remove
		"{window_sizes}mb_windows/variants/ins_{region}_{freq}_{fraction}p.bed",
		"{window_sizes}mb_windows/variants/del_{region}_{freq}_{fraction}p.bed",
		"{window_sizes}mb_windows/indels_{kmer}mer/frequency_{freq}_at_{fraction}p/ins_counts_{region}_{kmer}mer.bed", # remove
		"{window_sizes}mb_windows/indels_{kmer}mer/frequency_{freq}_at_{fraction}p/del_counts_{region}_{kmer}mer.bed", # remove
		"{window_sizes}mb_windows/indels_{kmer}mer/frequency_{freq}_at_{fraction}p/final/ins_counts_{region}_{kmer}mer.bed",
		"{window_sizes}mb_windows/indels_{kmer}mer/frequency_{freq}_at_{fraction}p/final/del_counts_{region}_{kmer}mer.bed",
		"{window_sizes}mb_windows/background_{kmer}mer/final/background_{region}_{kmer}mer_{fraction}p.bed",
		"{window_sizes}mb_windows/methylation_{fraction}p/{region}_methylation.bed",
		"{window_sizes}mb_windows/replication_timing_{fraction}p/{region}_replicationtime.bed", 
		"{window_sizes}mb_windows/size_difference_{freq}_{fraction}p/variant/ins_{region}_{size_partition}.bed",
		"{window_sizes}mb_windows/size_difference_{freq}_{fraction}p/variant/del_{region}_{size_partition}.bed",
		"{window_sizes}mb_windows/size_difference_{freq}_{fraction}p/tmp/ins_{region}_{size_partition}_{kmer}mer_tmp.bed",
		"{window_sizes}mb_windows/size_difference_{freq}_{fraction}p/tmp/del_{region}_{size_partition}_{kmer}mer_tmp.bed",
		"{window_sizes}mb_windows/size_difference_{freq}_{fraction}p/final/ins_{region}_{size_partition}_{kmer}mer_final.bed",
		"{window_sizes}mb_windows/size_difference_{freq}_{fraction}p/final/del_{region}_{size_partition}_{kmer}mer_final.bed",
		"{window_sizes}mb_windows/complex_structures/{region}_{complex_structure}_{fraction}p.bed",
		"{window_sizes}mb_windows/recombination/{region}_recombination_{fraction}p.bed"
		], region = regions, window_sizes = window_sizes, kmer = kmer_indels,
		chrom = chrom, fraction = NumberWithDepth, freq = allelefrequency, size_partition = size_partition,
		complex_structure = complex_structure)


rule coverage_regions:
	input:
		seq_summary_zipped = "{dataset}/coverage_files/{chrom}.BRAVO_TOPMed_coverage_hg38.txt.gz"
		seq_summary_unzipped = "{dataset}/coverage_files/{chrom}.BRAVO_TOPMed_coverage_hg38.txt.gz"
	params: lambda wildcards: float(int(wildcards.fraction)/100)
	output:
		bedfile = "coverage_bedfiles/{chrom}x10_{fraction}p.bed" #make depth value different?
		#all_chromosomes = "coverage_bedfiles/all_chromosomesx10_{fraction}p.bed"
	shell:"""
	gunzip -d {input.seq_summary}
	python topmed_rare.py coverage/{chrom}.BRAVO_TOPMed_coverage_hg38.txt {params.procent} > {output.bedfile}
	gzip coverage/{chrom}.BRAVO_TOPMed_coverage_hg38.txt
	"""
#echo {output.bedfile} > coverage_bedfiles/all_chromosomesx10_{wildcards.fraction}p.bed
###creatin vcf with alle indels under a give frequency

rule vcf_indel:
	input:
		raw_vcf = "raw_vcf/{chrom}.BRAVO_TOPMed_Freeze_8.vcf.gz"
	output:
		filtered = "raw_vcf/{chrom}_indel_{freq}.vcf.gz"
		#freqi = "raw_vcf/{chrom}_0.1.vcf.gz"
	shell:"""
	tabix -f -p vcf {input.raw_vcf}
	bcftools filter -O z -o {output.filtered} -i 'AF<{wildcards.freq} && VRT=2' {input.raw_vcf}
	"""
#bcftools filter -O z -o {output.indel} -i 'VRT=2' {output.freqi}
#cat {output.indel} raw_vcf/indel_0.1.vcf.gz > raw_vcf/indel_0.1.vcf.gz
rule catting_vcf:
	resources:
		threads=1,
		time=30,
		mem_mb=1000
	output:
		all_vcf = "raw_vcf/indel_{freq}.vcf.gz"
	shell:"""
	cat raw_vcf/*_indel_{wildcards.freq}.vcf.gz > {output.all_vcf}
	"""
rule catting_cov:
	resources:
		threads=1,
		time=30,
		mem_mb=1000
	output:
		all_cov = "coverage_bedfiles/all_chromosomesx10_{fraction}p.bed"
	shell:"""
	cat coverage_bedfiles/*x10_{wildcards.fraction}p.bed > {output.all_cov}
	"""
## a rule which makes the MegaBases bedfile 
## Creating regions which are to be investigated
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
		bedfiles = "{window_sizes}mb_windows/clean/{region}.bed"
	shell:"""
	printf '%s\t%s\t%s\n' {params.chrom} {params.start} {params.end} > {output.bedfiles}
	"""


rule checking_regions_with_coverage:
	input:
		bedfiles = "{window_sizes}mb_windows/clean/{region}.bed", #again do it by chromosome
		coverage_accepted = "coverage_bedfiles/all_chromosomesx10_{fraction}p.bed" #output?
	resources:
		threads=1,
		time=30,
		mem_mb=1000
	output:
		intersected_regions = "{window_sizes}mb_windows/intersected/intersected_{region}_{fraction}p.bed"
	shell:"""
	bedtools intersect -a {input.bedfiles} -b {input.coverage_accepted} > {output.intersected_regions}
	"""

rule checking_regions_with_blacklist:
	input:
		intersected_regions = "{window_sizes}mb_windows/intersected/intersected_{region}_{fraction}p.bed",
		blacklist = blacklist
	resources:
		threads=1,
		time=30,
		mem_mb=1000
	output:
		final_regions = "{window_sizes}mb_windows/final/final_{region}_{fraction}p.bed"
	shell:"""
	bedtools intersect -v -a {input.intersected_regions} -b {input.blacklist} > {output.final_regions}
	"""
rule accepting_regions:
	input:
		clean = "{window_sizes}mb_windows/clean/{region}.bed",
		final = "{window_sizes}mb_windows/final/final_{region}_{fraction}p.bed"
	resources:
		threads=1,
		time=30,
		mem_mb=1000
	output:
		number = "{window_sizes}mb_windows/accepted_regions/accepted_{region}_{fraction}p.bed"
	shell:"""
	tmp=`bedtools intersect -wo -a {input.clean} -b {input.final}| awk '{{s+=$7}} END {{print s}}'`
	num=$(expr {window_sizes} \* 1000000 / 2)
	if [[ $tmp -ge $num ]]
	then 
		cp {input.final} {output.number}
	else
		touch {output.number}
	fi
	"""

rule indel_background_counter: #im not sure this works
	input:
		accepted_regions = "{window_sizes}mb_windows/accepted_regions/accepted_{region}_{fraction}p.bed",
		genome = two_bit
	params:
		before_break  = lambda wildcards: int(creating_breakpoints(wildcards.kmer)[0]),
		after_break = lambda wildcards: int(creating_breakpoints(wildcards.kmer)[1])
	output:
		background = "{window_sizes}mb_windows/background_{kmer}mer/background_{region}_{kmer}mer_{fraction}p.bed",
		ss_background = "{window_sizes}mb_windows/background_{kmer}mer/final/background_{region}_{kmer}mer_{fraction}p.bed"
	shell:"""
	check=`cat {input.accepted_regions} | wc -l`
	if [[ $check -gt 0 ]]
	then 
		kmer_counter background --bed {input.accepted_regions} --before_after {params.before_break} {params.after_break} --reverse_complement_method both {input.genome} > {output.background}
		awk -v OFS='\t' '{{print "{wildcards.region}",$1,$2}}' {output.background} > {output.ss_background}
	else
		touch {output.background}
		touch {output.ss_background}
	fi
	"""

rule creating_indel_variants:
	input:
		accepted_regions = "{window_sizes}mb_windows/accepted_regions/accepted_{region}_{fraction}p.bed",
		vcf_file = "raw_vcf/indel_{freq}.vcf.gz" #change, make one or two rules that create this file(filter and ) output?
	output:
		variants =  "{window_sizes}mb_windows/variants/indels_{region}_{freq}_{fraction}p.bed", # can be removed # no it shouldnt 
		ins_variants = "{window_sizes}mb_windows/variants/ins_{region}_{freq}_{fraction}p.bed",
		del_variants = "{window_sizes}mb_windows/variants/del_{region}_{freq}_{fraction}p.bed"
	shell:"""
	check=`cat {input.accepted_regions} | wc -l`
	if [[ $check -gt 0 ]]
	then
		bedtools intersect -a {input.vcf_file} -b {input.accepted_regions} | awk -v OFS='\t' '{{print $1,$2,$4,$5}}' > {output.variants}
	else
		touch {output.variants}
	fi
	python creating_deletions.py {output.variants} > {output.del_variants}
	python creating_insertions.py {output.variants} > {output.ins_variants}
	"""

rule indel_variant_counter:
	input:
		ins_variants = "{window_sizes}mb_windows/variants/ins_{region}_{freq}_{fraction}p.bed",
		del_variants = "{window_sizes}mb_windows/variants/del_{region}_{freq}_{fraction}p.bed",
		genome = two_bit
	params:
		radius  = lambda wildcards: int(int(wildcards.kmer)/2)
	output:
		kmer_count_ins = "{window_sizes}mb_windows/indels_{kmer}mer/frequency_{freq}_at_{fraction}p/ins_counts_{region}_{kmer}mer.bed",
		kmer_count_del = "{window_sizes}mb_windows/indels_{kmer}mer/frequency_{freq}_at_{fraction}p/del_counts_{region}_{kmer}mer.bed",
		ss_ins = "{window_sizes}mb_windows/indels_{kmer}mer/frequency_{freq}_at_{fraction}p/final/ins_counts_{region}_{kmer}mer.bed",
		ss_del = "{window_sizes}mb_windows/indels_{kmer}mer/frequency_{freq}_at_{fraction}p/final/del_counts_{region}_{kmer}mer.bed"
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

### DOING METHYLATION DATA for plot

rule methylation_slicing: #im not sure this works
	input:
		accept_regions = "{window_sizes}mb_windows/accepted_regions/accepted_{region}_{fraction}p.bed",
		methylation = methylation_data
	resources:
		threads=1,
		time=60,
		mem_mb=1000
	output:
		#intersected_meth = "{window_sizes}mb_windows/methylation/intersected",
		ss_meth = "{window_sizes}mb_windows/methylation_{fraction}p/{region}_methylation.bed"
	shell:"""
	check=`cat {input.accept_regions}  | wc -l`
	if [[ $check -gt 0 ]]
	then 
		bedtools intersect -wa -a {input.methylation} -b {input.accept_regions} | awk -v OFS='\t' '{{print "{wildcards.region}",$2,$3,$5,$10,$11,$12,$13,$14}}' > {output.ss_meth}
	else
		printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' {wildcards.region} NA NA NA NA NA NA NA NA > {output.ss_meth}
	fi
	"""

##replication_time
rule replication_timing: #im not sure this works
	input:
		accept_regions = "{window_sizes}mb_windows/accepted_regions/accepted_{region}_{fraction}p.bed",
		rep_time = replication_time
	resources:
		threads=1,
		time=60,
		mem_mb=1000
	output:
		#intersected_meth = "{window_sizes}mb_windows/methylation/intersected",
		reptime = "{window_sizes}mb_windows/replication_timing_{fraction}p/{region}_replicationtime.bed"
	shell:"""
	check=`cat {input.accept_regions}  | wc -l`
	if [[ $check -gt 0 ]]
	then 
		bedtools intersect -wa -a {input.rep_time} -b {input.accept_regions} | awk -v OFS='\t' '{{print "{wildcards.region}",$0}}' > {output.reptime}
	else
		printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' {wildcards.region} NA NA NA NA NA > {output.reptime}
	fi
	"""

# DOING SIZE VARIANTS

rule size_indel_variants:
	input:
		ins_variants = "{window_sizes}mb_windows/variants/ins_{region}_{freq}_{fraction}p.bed",
		del_variants = "{window_sizes}mb_windows/variants/del_{region}_{freq}_{fraction}p.bed"
	params:
		#size = lambda wildcards: ",".join(wildcards.size_partition.split(","))
	output:
		insertion = "{window_sizes}mb_windows/size_difference_{freq}_{fraction}p/variant/ins_{region}_{size_partition}.bed",
		deletion = "{window_sizes}mb_windows/size_difference_{freq}_{fraction}p/variant/del_{region}_{size_partition}.bed"
	shell:"""
	check=`cat {input.ins_variants} | wc -l`
	if [[ $check -gt 0 ]]
	then
		python splitting_indel_size.py {input.ins_variants} {wildcards.size_partition} > {output.insertion} 
		python splitting_indel_size.py {input.del_variants} {wildcards.size_partition} > {output.deletion}
		
	else
		touch {output.insertion}
		touch {output.deletion}
	fi
	"""

rule indel_size_counter:
	input:
		insertion_size = "{window_sizes}mb_windows/size_difference_{freq}_{fraction}p/variant/ins_{region}_{size_partition}.bed",
		deletion_size = "{window_sizes}mb_windows/size_difference_{freq}_{fraction}p/variant/del_{region}_{size_partition}.bed",
		genome = two_bit
	params:
		radius  = lambda wildcards: int(int(wildcards.kmer)/2)
	output:
		insertion_sized = "{window_sizes}mb_windows/size_difference_{freq}_{fraction}p/tmp/ins_{region}_{size_partition}_{kmer}mer_tmp.bed",
		deletion_sized = "{window_sizes}mb_windows/size_difference_{freq}_{fraction}p/tmp/del_{region}_{size_partition}_{kmer}mer_tmp.bed",
		insertion_size_ss = "{window_sizes}mb_windows/size_difference_{freq}_{fraction}p/final/ins_{region}_{size_partition}_{kmer}mer_final.bed",
		deletion_size_ss = "{window_sizes}mb_windows/size_difference_{freq}_{fraction}p/final/del_{region}_{size_partition}_{kmer}mer_final.bed"
	shell:"""
	check=`cat {input.insertion_size} | wc -l`
	if [[ $check -gt 0 ]]
	then
		kmer_counter indel -r {params.radius} --sample {input.genome} {input.insertion_size} ins > {output.insertion_sized}
		kmer_counter indel -r {params.radius} --sample {input.genome} {input.deletion_size} del_start > {output.deletion_sized}
		awk -v OFS='\t' '{{print "{wildcards.region}","{wildcards.size_partition}",$1,$2}}' {output.insertion_sized} > {output.insertion_size_ss} 
		awk -v OFS='\t' '{{print "{wildcards.region}","{wildcards.size_partition}",$1,$2}}' {output.deletion_sized} > {output.deletion_size_ss}

	else
		touch {output.insertion_sized}
		touch {output.deletion_sized}
		touch {output.deletion_size_ss}
		touch {output.insertion_size_ss}
	fi
	"""

####COMPLEX STRUCTURES#####

rule creating_complex_structure:
	input:
		accepted_regions = "{window_sizes}mb_windows/accepted_regions/accepted_{region}_{fraction}p.bed",
		complex_structure = "nonBdna/{complex_structure}.bed" #change, make one or two rules that create this file(filter and ) output?
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

####Recombination rate######

rule recombination_rate:
	input:
		accepted_regions = "{window_sizes}mb_windows/accepted_regions/accepted_{region}_{fraction}p.bed",
		recombination_map = recombination #change, make one or two rules that create this file(filter and ) output?
	output:
		recomb =  "{window_sizes}mb_windows/recombination/{region}_recombination_{fraction}p.bed"
	shell:"""
	check=`cat {input.accepted_regions} | wc -l`
	if [[ $check -gt 0 ]]
	then
		bedtools intersect -a {input.recombination_map} -b {input.accepted_regions} | awk -v OFS='\t' '{{print "{wildcards.region}",$0}}' > {output.recomb}
	else
		touch {output.recomb}
	fi
	"""