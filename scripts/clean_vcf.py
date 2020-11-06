import pysam
import sys

def read_map(map):
    markers = {}
    with open(map, "r") as mapfile:
        map_data = mapfile.readlines()
        i = 0
        for item in map_data: 
            if i > 7:
                try:
                    line = item.strip("\n").split(",")
                    markers[line[1]] = (line[9], line[10], line[3].split("/")[0].strip("["), line[3].split("/")[1].strip("]"))
                except IndexError:
                    pass
            i += 1
    return(markers)

def mod_vcf_in_place(vcf, markers):
    with open(vcf, "r") as vcffile:
        vcf_data = vcffile.readlines()
        for item in vcf_data:
            if item[0] == "#":
                print(item.strip("\n"))
                continue
            line = item.strip("\n").split("\t")
            try:
                new_data = markers[line[2]]
                if new_data[0] == "XY":
                    pass
                else:
                    line[0] = new_data[0]
                    line[1] = new_data[1]
            except KeyError:
                continue
            if str(line[0]) != "0" and str(line[1]) != "0" and str(line[3]) not in ["I", "D"]:
                print("\t".join(line))

if __name__ == "__main__":
    markers = read_map(sys.argv[1])
    mod_vcf_in_place(sys.argv[2], markers)
