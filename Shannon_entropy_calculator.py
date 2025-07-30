import math
import sys
from matplotlib import pyplot as plt
import gzip
import re

def shannon_entropy(column):
    """Calculate Shannon entropy for a column of nucleotides"""
    counts = {}
    total = 0
    for base in column:
        """I don't think N should be in this list as it already means any base so is uninformative"""
        if base not in 'ACGT-':
            continue  # skip invalid characters
        counts[base] = counts.get(base, 0) + 1
        total += 1
    if total == 0:
        return 0.0
    entropy = 0.0
    for count in counts.values():
        p = count / total
        entropy -= p * math.log2(p)

    return entropy

def read_alignment(filename):
    """Check if alignment file is gzipped or not"""
    if re.search(".gz", filename):
        with gzip.open(filename,"rt") as f:
            sequences = sequence_reader(f)
            return(sequences)
    else:
        with open(filename,"r") as f:
            sequences = sequence_reader(f)
            return(sequences)
    

def sequence_reader(f,sequences=[]):
    """Read alignment file in FASTA format and return list of sequences"""   
    seq = ''
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            if seq:
                sequences.append(seq)
                seq = ''
        else:
            seq += line.upper()
    if seq:
        sequences.append(seq)
    return sequences

def calculate_entropy_matrix(sequences):
    """Calculate entropy for each column in the alignment"""
    if not sequences:
        return []
    length = len(sequences[0])
    entropies = []
    for i in range(length):
        column = [seq[i] for seq in sequences if len(seq) > i]
        entropies.append(shannon_entropy(column))
    
    return entropies

def plotter(entropies,title,smooth):
    """Plotting the entropies in a line graph"""
    entropies_dict={}
    for n,i in enumerate(entropies):
        if not smooth:
            entropies_dict[n] = i
        else:
            """Smoothing works by averageing the entropy across a window of size smooth"""
            if n < smooth:
                entropies_dict[n] = (sum(entropies[0:(n+int(smooth/2))]))/(1+n+int(smooth/2))
            elif n > (len(entropies) - smooth):
                entropies_dict[n] = (sum(entropies[(n-int(smooth/2)):len(entropies)]))/((len(entropies)-n)+(int(smooth/2)))
            else:
                entropies_dict[n] = (sum(entropies[(n-int(smooth/2)):(n+int(smooth/2))]))/(smooth)
    
    plt.axes(xlabel="Amplified Region",ylabel="Entropy",ylim=(0,2.5))  
    plt.plot(entropies_dict.keys(),entropies_dict.values(),ls='-')
    plt.hlines(2.321,0,len(entropies_dict)) # max entropy for 5 items is 2.321
    plt.savefig(title)
    #plt.show()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python entropy_calculator.py <alignment.fasta> <title> <smooth_factor>")
        sys.exit(1)
    sequences = read_alignment(sys.argv[1])
    try: 
        title = sys.argv[2]
    except:
        title = sys.argv[1].split(".")[0]
    try: 
        smooth = int(sys.argv[3])
    except:
        smooth = 0
    
    entropies = calculate_entropy_matrix(sequences)
    plotter(entropies,title,smooth)
