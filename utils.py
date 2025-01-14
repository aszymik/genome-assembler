from typing import List
from Bio import SeqIO

def load_reads(reads_path: str) -> List[str]:
    reads = []
    with open(reads_path, 'r') as reads_file:
        for record in SeqIO.parse(reads_file, 'fasta'):
            reads.append(str(record.seq))
    return reads

def save_contigs(contigs: List[str], output_path: str) -> None:
    with open(output_path, 'w') as out_file:
        for i, contig in enumerate(contigs):
            out_file.write(f'>contig_{i}\n')
            out_file.write(contig + '\n')
            