# Genome assembler

Genome assembly is a fundamental task in bioinformatics, enabling the reconstruction of complete DNA sequences from fragmented short reads produced by high-throughput sequencing technologies like Illumina sequencing.

This assembler employs a de Bruijn graph-based approach, which is particularly effective for short-read assembly. 
This workflow provides an efficient method for assembling sequences from short reads, incorporating error correction and graph simplification for improved accuracy.

## Pipeline

### Load reads

### K-mer correction

Replace rare k-mers with their more frequent version. To do that for each k-mer which don't
pass the expected frequency threshold we collect its neighbours (k-mers that have Hammming distance = 1)
and loop through it replacing it with first kmer that meet the criteria.

### Build DeBruijn graph

Loop through all reads and add they k-1 mers to the graph

### Remove divergent branches

Remove short divergent paths, which where caused by some more frequent sequencing error.

This approach removes "tip" nodes from a graph by iteratively simplifying the structure. A tip is defined as a node with either no outgoing edges (forward tip) or no incoming edges (reverse tip), where adjacent edge weights meet a specified threshold.

* Forward Tips: Nodes with zero outgoing edges (nout == 0) are identified. For these nodes, incoming edges with weights below or equal to a specified threshold are removed, and the nodes themselves are deleted from the graph.

* Reverse Tips: Similarly, nodes with zero incoming edges (nin == 0) are located. Outgoing edges from these nodes that meet the weight condition are removed, followed by the deletion of the nodes.

This process is repeated iteratively until no more tip nodes that satisfy the criteria remain. The method ensures the graph is simplified by removing low-weight, dangling nodes that do not contribute significantly to the overall structure. 

### Contig extraction

contigs are constructed from a de Bruijn graph by greedily traversing its nodes to form continuous sequences. The process involves the following steps:

* **Identify Starting Points**: Nodes in the graph are iterated over, and each unvisited node is treated as the potential start of a new contig.

* **Build Contigs**: From each starting node, the method extends the sequence by following outgoing edges to unvisited neighboring nodes. Each step adds the overlapping portion of the neighbor’s sequence to the contig.

* **Stop Conditions**: The traversal ends when a node has no unvisited outgoing edges or all its neighbors have already been included in another contig.

* **Repeat**: The process continues until all nodes in the graph are visited, ensuring every possible contig is extracted.

## Example usage

Before running make sure to install all requirements listed in `requirements.txt`

### Assemble
`./assambly <path_to_reads> <output_path>`

`python assambly.py <path_to_reads> <output_path>`

### Evaluate

1 `./evaluate.sh <path_to_file_containing_contigs>`