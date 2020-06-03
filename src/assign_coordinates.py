#!/usr/bin python3
import pandas as pd
from tqdm import tqdm


# Declare chromosome numbers
chromosomes = {
    "I": 1,
    "II": 2,
    "III": 3
}

# File routes
data_file = "/media/danilexn/danilexn/Lab010 Alfonso/SGA_Alfonso/data/allPictures/screen02/1020/Day3/Copy2/S001.05d3.2_non-filtered.txt"
coor_file = "/home/danilexn/Escritorio/genes.gff3"

# Parse data and preprocess. Reduce dataframes
data = pd.read_csv(data_file, sep = "\t")
coor = pd.read_csv(coor_file, sep = "\t", header = None) # gff3 file without initial header
coor = coor[[0,3,4,8]].set_index(8)
genelist = [g for g in list(data['ORF'].values) if isinstance(g, str) and "wild-type" not in g]

# Main iterator
for i, gene in enumerate(tqdm(genelist)):
    try:
        pos = coor.index.get_loc(gene)
        chromo, start, end = coor.iloc[pos,[0,1,2]]
        data.at[i,'start'] = int(start)
        data.at[i,'end'] = int(end)
        data.at[i,'chromosome'] = int(chromosomes[chromo])
    except Exception as e:
        print("ORF {} could not be assigned".format(gene))

data.to_csv("/home/danilexn/Escritorio/crossed.tsv", sep = "\t")