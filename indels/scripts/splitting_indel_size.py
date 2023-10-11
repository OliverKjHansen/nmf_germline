import sys

def splitting_into_sizes(file, interval):
    with open(file, "r") as f:
        for line in f:
            variant_indel = line.split("\n")[0].split("\t")
            if len(variant_indel[2]) > len(variant_indel[3]): # this is a deletion
                vs = len(variant_indel[2])-1
                if vs >= int(interval[0]):
                    if vs <= int(interval[1]):
                        print(line.split("\n")[0])
                    if int(interval[1]) == 0:
                        print(line.split("\n")[0])
            if len(variant_indel[2]) < len(variant_indel[3]): # this is a insertion
                vs = len(variant_indel[3])-1
                if vs >= int(interval[0]):
                    if vs <= int(interval[1]):
                        print(line.split("\n")[0])
                    if int(interval[1]) == 0:
                        print(line.split("\n")[0])
                
if __name__ == '__main__':
    str_list_input = sys.argv[2:]
    if len(str_list_input) > 1:
        int_list_input = list(str_list_input)
        size_interval = [int(x) for x in int_list_input]
    if len(str_list_input) == 1:
        str_list_input = ",".join(str_list_input[0].split())
        size_interval = [int(x) for x in str_list_input]
    variants = sys.argv[1]
    splitting_into_sizes(variants, size_interval)
