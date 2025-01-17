import DBG
import utils
import numpy as np
import math

def main(reads_path: str, output_path: str, k: int=19, expected_coverage: float=6.0):
    reads = utils.load_reads(reads_path)
    kmerhist = DBG.kmerHist(reads, k)  # generate k-mer histogram
    kmer_thresh = math.floor(np.mean(list(kmerhist.values())))  # calculate k-mer thershold for read correction
    reads_corrected = [DBG.correct1mm(read, k, kmerhist, kmer_thresh) for read in reads]  # change infrequent to frequent k-mers
    
    dbg = DBG.DeBruijnGraph(reads_corrected, k)

    weights = []
    for node in dbg.G:
        for neighbor in dbg.G[node]:
            weights.append(dbg.G[node][neighbor])

    weight_thresh = max(math.floor(np.mean(weights) - expected_coverage), 0)

    dbg.remove_tips(weight_thresh)
    contigs = DBG.extract_contigs_greedy(dbg)

    # graph = DBG.DeBruijnGraph2(reads_corrected, k).to_dot(weights=True)
    # graph.render(filename='graphs/reads1_graph_without_tips', format='png', cleanup=True)

    utils.save_contigs(contigs, output_path)
    
if __name__ == '__main__':
    main('training/reads/reads1.fasta', 'outs/reads1_contigs')
    main('training/reads/reads2.fasta', 'outs/reads2_contigs')
    main('training/reads/reads3.fasta', 'outs/reads3_contigs')