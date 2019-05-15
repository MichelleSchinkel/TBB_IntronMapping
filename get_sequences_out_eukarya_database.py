# retrieving specific sequences from the eukarya4-database

# first open the file containing the OGs and their corresponding eukarya4-identifiers
tubulin_doc = open("/home/michelle/Documents/project/data/tubulin/tubulin_LECA-OGs.tsv", "r")

tubulin_dic ={}		# dictionary which will contain eukarya4-identifiers as keys and corresponding OGs as values
genes = []		# list in which all the OGs will be stored 
sequence_ids = []	# list in which all the eukarya4-identifiers will be stored

for line in tubulin_doc:
    line = line.strip()
    words = line.split()
    tubulin_dic[words[1]] =str(words[0])
    if words[0] not in genes:
        genes.append(words[0])
    if words[1] not in sequence_ids:
        sequence_ids.append(words[1])
        
tubulin_doc.close()

# if your want to do some checks if everything is in there
#print(genes)
#print(len(sequence_ids))

# open the eukarya4-database and create a file in which output will be written
eukarya = open("/home/julian/julian2/snel-clan-genomes/eukarya/eukarya_4.0.1/eukarya_4.0.1_proteomes_lt.fa", "r")
outfile = open("/home/michelle/Documents/project/data/tubulin/tubulin_sequences.fasta", "w")

# import the SeqIO module from Bio, make sure Bio is installed with pip install Bio
from Bio import SeqIO

# look into the eukarya4-database, 
# see if the header contains an OG present in our sequence_ids list
# if yes, print the header in the desired format >OG_eukarya4-identifier
# if yes, print under header the corresponding sequence
for record in SeqIO.parse(eukarya, 'fasta'):
    for item in sequence_ids:
        if item.strip() == record.id:
            outfile.write(">" + tubulin_dic[record.id] + "_" + record.id + "\n")
            outfile.write(str(record.seq) + "\n")

# close the output file
outfile.close()
