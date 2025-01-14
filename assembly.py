import DBG
import utils

def main(reads_path: str, output_path: str, k: int, kmer_thresh: int=3):
    reads = utils.load_reads(reads_path)
    kmerhist = DBG.kmerHist(reads, k)  # generate k-mer histogram
    reads_corrected = [DBG.correct1mm(read, k, kmerhist, kmer_thresh) for read in reads]  # change infrequent to frequent k-mers
    dbg = DBG.DeBruijnGraph(reads_corrected, k)

    # TODO:
    # remove tips
    # remove bubbles

    contigs = DBG.extract_contigs_greedy(dbg)
    utils.save_contigs(contigs, output_path)
    
if __name__ == '__main__':
    main('training/reads/reads1.fasta', 'outs/reads1_contigs', 21)