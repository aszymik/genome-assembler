import DBG
import utils
import numpy as np
import argparse

def main(reads_path: str, output_path: str, k: int=19, expected_coverage: float=6.0):

    # Load reads and correct them using k-mer frequencies
    reads = utils.load_reads(reads_path)
    kmerhist = DBG.kmerHist(reads, k)
    kmer_thresh = np.floor(np.mean(list(kmerhist.values())))  # calculate k-mer thershold for read correction
    reads_corrected = [DBG.correct1mm(read, k, kmerhist, kmer_thresh) for read in reads]  # change infrequent to frequent k-mers
    
    # Generate de Bruijn graph from corrected reads
    dbg = DBG.DeBruijnGraph(reads_corrected, k)

    # Remove tips from the graph
    mean_weight = dbg.get_mean_edge_weight()
    weight_thresh = max(np.floor(mean_weight - expected_coverage), 0)
    dbg.remove_tips(weight_thresh)

    # Extract contigs and save them to file
    contigs = DBG.extract_contigs_greedy(dbg)
    utils.save_contigs(contigs, output_path)
    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Assemble reads into contigs using a de Bruijn graph approach.")
    parser.add_argument("reads_path", type=str, help="Path to the input reads in FASTA format.")
    parser.add_argument("output_path", type=str, help="Path to the output contigs in FASTA format.")
    parser.add_argument("-k", type=int, default=19, help="k-mer size (default: 19).")
    parser.add_argument("-c", "--coverage", type=float, default=6.0, help="Expected coverage for tip removal (default: 6.0).")

    args = parser.parse_args()

    main(args.reads_path, args.output_path, args.k, args.coverage)