import sys
from Bio import SeqIO, Seq
from collections import defaultdict

def kmers(str,k): # make kmers
    for i in range(len(str) - k + 1):
        yield str[i:i+k]


def BuildDeBrujinGraph(reads, k=25): #built a graph
    graph = defaultdict(list)
    all_kmers = set()
    for read in reads:
        for kmer in kmers(read, k):
            head = kmer[:-1]
            tail = kmer[1:]
            graph[head].append(tail)

    return graph

def Tour(start, graph, end=None): #generate a tour in graph
    tour = [start]

    if end is None:
        finish = start
    else:
        finish = end
    while True:
        if tour[-1] not in graph.keys():
            break

        nodes = graph[tour[-1]]
        if not nodes: # case not eulerian circle
            break

        tour.append(nodes.pop())
        if tour[-1] == finish: # case eulerian circle
            break

    offset = 0
    for i,step in enumerate(tour):
        try:
            nodes = graph[step]
            if nodes:
                tour_ = Tour(nodes.pop(), graph, step) # generate a subtour
                i += offset
                tour = tour[ : i + 1 ] + tour_ + tour[ i + 1 : ]
                offset += len(tour_)
                #print("-----------")

        except:
            continue

    return tour

def BuildGenomeFromDeBrujnGraph(graph):
    start = None
    for key in list(graph.keys()):
        if len(graph[key]) % 2 == 0 and len(graph[key]) != 0:
            start = key
            break
    print(f"START: {start}")
    tour = Tour(start, graph)
    return "".join([tour[0]] + [s[-1] for s in tour[1:]])


if __name__ == "__main__":
    reads_l = [seqrecord.seq._data.upper() for seqrecord in list(SeqIO.parse("Carsonella_ruddii_reads_paired_reads_left.fastq", "fastq"))]
    reads_r = [seqrecord.seq._data.upper() for seqrecord in list(SeqIO.parse("Carsonella_ruddii_reads_paired_reads_right.fastq", "fastq"))]
    all_reads = list(set(reads_l + reads_r))
    graph = dict(BuildDeBrujinGraph(all_reads))
    genome =  BuildGenomeFromDeBrujnGraph(graph)
    with open("build_genome.txt", "w") as f:
        f.write(genome)
