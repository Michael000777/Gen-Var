import bamnostic as bs
from collections import Counter

# Pileup Object 
class Pileup:
    """
    Contains a counter of nulceotides aligned at a genomic position and computes depth and consensus for a genomic position
    
    Attributes:
    counts (Counter):  counter of all nucleotides aligned at the genomic position
    depth (int): total number of reads aligned at the genomic position
    consensus (str): most common nucleotide
    
    Methods:
    update(str): update the counter with the given string which is a nucleotide or a sequence to add to the pileup              
    compute_frequency: compute the ratio between the frequency of the most common nucleotide and the total number of nucleotides
    """
    
    def __init__(self, counts = None):
        self.counts = counts
        if self.counts == None:
            self.counts = Counter()
        self.__depth = sum(self.counts.values())
        if self.__depth == 0:
            self.__consensus = None        
        else:
            self.__consensus = self.counts.most_common()[0][0]

    def __str__(self):
        return f"Pileup(counts = {self.counts})"
    
    def __repr__(self):
        return f"Pileup(counts = {self.counts})"
    
    @property # getter
    def consensus(self):
        """
        Get the consensus nucleotide for the pileup
        """
        return self.__consensus

    @property # getter
    def depth(self):
        """
        Get the dept for the pileup, that is the total number of nucleotides
        """
        return self.__depth      
        
    def update(self, seq):
        """
        Given new nucleotides to add to the pileup, update the counts depth and consensus nucleotide 
        """
        self.counts.update(seq)
        self.__depth = sum(self.counts.values())
        self.__consensus = self.counts.most_common()[0][0]
        
    def compute_frequency(self):
        """
        Compute the ratio between the frequency of the most common nucleotide and the total number of nucleotides
        Returns None if there are no nucleotides in the Pileup
        """
        if self.__depth != 0:
            return self.counts[self.consensus]/self.__depth
        


def initialize_positions(genome_filename):
    genome_positions = {'normal':[], 'tumor':[]}
    # FILL IN THE CODE 
    
    with open(genome_filename, "r") as ref_gen:
        for line in ref_gen:
            if not line.startswith(">"):
                for _ in line:
                    if _ != "\n":
                        genome_positions['normal'].append(Pileup())
                        genome_positions['tumor'].append(Pileup())
                    
                
    
        
    return genome_positions


def process_read(read_pos, read_seq, sample_name, genome_positions = None):
    # FILL IN THE CODE 
    pos = int(read_pos)

    for nuc in read_seq:
        if nuc != "\n":
            p =  genome_positions[sample_name][pos]
            p.update(nuc)

            pos += 1

    return genome_positions


def process_bam(filename, sample_name, genome_positions = None):
    with bs.AlignmentFile(filename) as bam:
        reads_no = 0
        for read in bam:
            process_read(read.pos, read.seq, sample_name, genome_positions)
            
            # we keep track of the number of read processed to display progress
            reads_no += 1
            if reads_no % 20000 == 0:
                    print(f"Processed {reads_no} reads!")
    return genome_positions


def process_genomic_data(genome_positions = None):
    variant_calls = []
    # FILL IN THE CODE 
        
    for pos, (normal, tumor) in enumerate(zip(genome_positions['normal'], genome_positions['tumor'])):
        
        x = tumor.consensus
        cons_tumor = x
        cons_normal = normal.consensus
        
        if tumor.depth != 0:
            if normal.depth != 0:
                if cons_normal != cons_tumor:

                    counts_variant_base_calls = tumor.counts[x]
                    total_base_counts = tumor.depth

                    alelle_freq = counts_variant_base_calls/total_base_counts
                    if alelle_freq > 0.5:
                        variant_calls.append((pos, cons_tumor, cons_normal, alelle_freq))
        
                

    return variant_calls

# Initializing the list
print("Initializing genome positions")
genome_positions = initialize_positions('b_subtilis_genome.fa')
print("done initialization")


# Process all the bam files
for filename in ('normal.bam', 'tumor.bam'):
    print("Processing bam", filename)
    genome_positions = process_bam(filename, filename.split('.')[0], genome_positions)
    print("done processing bam", filename)
    print()

print("Processing genome positions")
results = process_genomic_data(genome_positions)
print("done process genomic data")

#Visualizing first 20 results 
results[:20]
