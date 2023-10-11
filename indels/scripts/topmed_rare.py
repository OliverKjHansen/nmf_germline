import sys
import pandas as pd

#making file that transform coverage file from topmed into a bedfile with certain fraction over a certain threshold of coverage
# in bed file format. bedfiles are 0 indexed
#percent_of_individual with coverage over 

#header could be more elegant

def making_bed(file, percent_of_individual):
    with open(file, "r") as f:
        start = 0 
        end = 0
        pos_counter = 0
        for idx, line in enumerate(f):
            if idx == 2:
                chrom = int(line.split("\t")[0])
                break
        for line in f:
            if line.split("\t")[0] == "CHROM":
                continue
            pos = int(line.split("\t")[1])
            percent = float(line.split("\t")[6])
            if percent < float(percent_of_individual):
                if start != end:
                    end = pos    #stop region 
                    if (end-start) >= 1:
                        print("\t".join(["chr"+str(chrom),str((start-1)),str(end),str((end-start))]))
                        pos_counter += (end-start)
                    start = end
                if start == end: #continue searching for region 
                    continue
            if percent >= float(percent_of_individual):
                if (start != end) and (pos-old_pos == 1): # continue region
                    end = pos
                    old_pos = pos
                if start == end:
                    old_pos = pos # start region 
                    start = pos
                if (start != end) and (pos-old_pos > 1):
                    print("\t".join(["chr"+str(chrom),str((start-1)),str(end),str((end-start))])) #to make the output zerobased
                    pos_counter += (end-start)
                    start = end
        if (end-start != 0):
            print("\t".join(["chr"+str(chrom),str((start-1)),str(end),str((end-start))]))
            pos_counter += (end-start)
    print(pos_counter)
                
def counting_breakpoints(file, percent_of_individual):
    with open(file, "r") as f:
        counter = 0
        old_pos = 1
        jump_counter = 0
        for line in f:
            if line.split("\t")[0] == "CHROM":
                continue
            percent = float(line.split("\t")[6])
            new_pos = int(line.split("\t")[1])
            if (new_pos-old_pos) > 1:
                #print(new_pos,old_pos)
                jump_counter += (new_pos-old_pos)
            old_pos = int(line.split("\t")[1])
            if percent >= float(percent_of_individual):
                counter += 1
    print(counter, jump_counter)

if __name__ == '__main__':
    chr_cov = sys.argv[1]
    percent_of_individual = sys.argv[2]
    making_bed(chr_cov, percent_of_individual)
    counting_breakpoints(chr_cov, percent_of_individual)

