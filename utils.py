from typing import List

def load_reads(reads_path: str) -> List[str]:
    reads = []
    with open(reads_path, 'r') as reads_file:
        for line in reads_file:
            if not line.startswith('>'):
                reads.append(line.strip())
    return reads