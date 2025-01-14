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
        for contig in contigs:
            out_file.write(contig + '\n')

print(load_reads('training/reads/reads1.fasta'))