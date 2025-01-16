import DBG
import utils

def main(reads_path: str, output_path: str, k: int, kmer_thresh: int=2, weight_thresh: int=1, similarity_thresh: int=3):
    reads = utils.load_reads(reads_path)
    kmerhist = DBG.kmerHist(reads, k)  # generate k-mer histogram
    reads_corrected = [DBG.correct1mm(read, k, kmerhist, kmer_thresh) for read in reads]  # change infrequent to frequent k-mers
    
    dbg = DBG.DeBruijnGraph(reads_corrected, k)
    dbg.remove_tips(weight_thresh)
    # dbg.remove_bubbles()
    # dbg = DBG.remove_bubbles(reads, dbg, k, similarity_thresh)

    contigs = DBG.extract_contigs_greedy(dbg)

    # graph = DBG.DeBruijnGraph2(reads_corrected, k).to_dot(weights=True)
    # graph.render(filename='graphs/reads1_graph_without_tips', format='png', cleanup=True)

    utils.save_contigs(contigs, output_path)
    
if __name__ == '__main__':
    main('training/reads/reads3.fasta', 'outs/reads3_contigs', 15)