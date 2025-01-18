import numpy as np

def kmerHist(reads, k):
    """ Return k-mer histogram and average # k-mer occurrences """
    kmerhist = {}
    for read in reads:
        for kmer in [ read[i:i+k] for i in range(len(read)-(k-1)) ]:
            kmerhist[kmer] = kmerhist.get(kmer, 0) + 1
    return kmerhist


def neighbors1mm(kmer, alpha):
    """ Generate all neighbors at Hamming distance 1 from kmer """
    neighbors = []
    for j in range(len(kmer)-1, -1, -1):
        oldc = kmer[j]
        for c in alpha:
            if c == oldc: continue
            neighbors.append(kmer[:j] + c + kmer[j+1:])
    return neighbors

def correct1mm(read, k, kmerhist, thresh, alpha=['A','C', 'G','T']):
    """
    Return an error-corrected version of read.  

    Args:
        k: k-mer length.
        kmerhist: k-mer count map.
        thresh: count threshold above which k-mer is considered correct.
        alpha: alphabet.
    """
    # Iterate over k-mers in read
    for i in range(len(read)-(k-1)):
        kmer = read[i:i+k]
        # If k-mer is infrequent...
        if kmerhist.get(kmer, 0) <= thresh:
            # Look for a frequent neighbor
            for newkmer in neighbors1mm(kmer, alpha):
                if kmerhist.get(newkmer, 0) > thresh:
                    # Found a frequent neighbor; replace old kmer
                    # with neighbor
                    read = read[:i] + newkmer + read[i+k:]
                    break
    # Return possibly-corrected read
    return read

class DeBruijnGraph:
    """ De Bruijn directed multigraph built from a collection of
        strings. User supplies strings and k-mer length k.  Nodes
        are k-1-mers.  An Edge corresponds to the k-mer that joins
        a left k-1-mer to a right k-1-mer. """

    @staticmethod
    def chop(st, k):
        """ Chop string into k-mers of given length """
        for i in range(len(st)-(k-1)):
            yield (st[i:i+k], st[i:i+k-1], st[i+1:i+k])

    class Node:
        """ Node representing a k-1 mer.  Keep track of # of
            incoming/outgoing edges so it's easy to check for
            balanced, semi-balanced. """

        def __init__(self, km1mer):
            self.km1mer = km1mer
            self.nin = 0
            self.nout = 0
            self.parents = set()

        def isSemiBalanced(self):
            return abs(self.nin - self.nout) == 1

        def isBalanced(self):
            return self.nin == self.nout

        def __hash__(self):
            return hash(self.km1mer)

        def __str__(self):
            return self.km1mer

    def __init__(self, strIter, k, circularize=False):
        """ Build de Bruijn multigraph given string iterator and k-mer
            length k """
        self.G = {}     # Forward graph: {Node: {neighbor: weight}}
        self.reverse_G = {}  # Reverse graph: {Node: {parent: weight}}
        self.nodes = {} # Maps k-1-mers to Node objects
        
        for st in strIter:
            if circularize:
                st += st[:k - 1]
            for kmer, km1L, km1R in self.chop(st, k):
                nodeL = self.nodes.setdefault(km1L, self.Node(km1L))
                nodeR = self.nodes.setdefault(km1R, self.Node(km1R))
                
                # Update edge weights in forward graph
                if nodeL not in self.G:
                    self.G[nodeL] = {}
                self.G[nodeL][nodeR] = self.G[nodeL].get(nodeR, 0) + 1

                # Update edge weights in reverse graph
                if nodeR not in self.reverse_G:
                    self.reverse_G[nodeR] = {}
                self.reverse_G[nodeR][nodeL] = self.reverse_G[nodeR].get(nodeL, 0) + 1

                # Update node properties
                nodeL.nout += 1
                nodeR.nin += 1
                nodeR.parents.add(nodeL)

        # Iterate over nodes; tally # balanced, semi-balanced, neither
        self.nsemi, self.nbal, self.nneither = 0, 0, 0
        self.head, self.tail = None, None
        for node in self.nodes.values():
            if node.isBalanced():
                self.nbal += 1
            elif node.isSemiBalanced():
                if node.nin == node.nout + 1:
                    self.tail = node
                if node.nin == node.nout - 1:
                    self.head = node
                self.nsemi += 1
            else:
                self.nneither += 1

    def __str__(self):
        """ For debugging: print the forward and reverse graphs. """
        forward = "\n".join(
            f"{str(u)} -> {str(v)} [weight={weight}]"
            for u, neighbors in self.G.items()
            for v, weight in neighbors.items()
        )
        reverse = "\n".join(
            f"{str(v)} <- {str(u)} [weight={weight}]"
            for v, parents in self.reverse_G.items()
            for u, weight in parents.items()
        )
        return f"Forward Graph:\n{forward}\n\nReverse Graph:\n{reverse}"
    
    def get_mean_edge_weight(self):
        """ Return the mean of weights of all edges in the graph. """
        weights = []
        for node in self.G:
            for neighbor in self.G[node]:
                weights.append(self.G[node][neighbor])
        return np.mean(weights)

    def get_ancestors(self, node):
        """ Return the ancestors of a node. """
        return list(self.reverse_G.get(node, {}).keys()) 
    
    def remove_tips(self, weight_threshold):
        """ Remove tips from the graph.
        
        Args:
            weight_threshold (int): The maximum weight for an edge to be considered a tip.
        """
        removed = True  # to keep track of whether we removed a tip in the last iteration

        while removed:
            removed = False

            # Find and remove forward tips (nodes with 0 outgoing edges)
            forward_tips = [node for node in self.nodes.values() 
                            if node.nout == 0 and any(w <= weight_threshold for w in self.reverse_G.get(node, {}).values())]
            for tip in forward_tips:
                # Remove edges pointing to this tip
                for parent in self.get_ancestors(tip):
                    weight = self.reverse_G[tip].pop(parent, 0)
                    if weight:
                        self.G[parent][tip] -= weight
                        if self.G[parent][tip] == 0:
                            del self.G[parent][tip]
                        parent.nout -= 1
                # Remove the tip node
                try:
                    del self.reverse_G[tip]
                    del self.G[tip]
                    del self.nodes[tip.km1mer]
                except KeyError:
                    pass
                removed = True

            # Find and remove reverse tips (nodes with 0 incoming edges)
            reverse_tips = [node for node in self.nodes.values()
                            if node.nin == 0 and any(w <= weight_threshold for w in self.G.get(node, {}).values())]
            for tip in reverse_tips:
                # Remove edges originating from this tip
                for child, weight in list(self.G.get(tip, {}).items()):
                    if weight:
                        self.G[tip].pop(child, 0)
                        self.reverse_G[child][tip] -= weight
                        if self.reverse_G[child][tip] == 0:
                            del self.reverse_G[child][tip]
                        child.nin -= 1
                # Remove the tip node
                try:
                    del self.G[tip]
                    del self.reverse_G[tip]
                    del self.nodes[tip.km1mer]
                except KeyError:
                    pass
                removed = True

    def nnodes(self):
        """ Return # nodes """
        return len(self.nodes)

    def nedges(self):
        """ Return # edges """
        return len(self.G)

    def hasEulerianWalk(self):
        """ Return true iff graph has Eulerian walk. """
        return self.nneither == 0 and self.nsemi == 2

    def hasEulerianCycle(self):
        """ Return true iff graph has Eulerian cycle. """
        return self.nneither == 0 and self.nsemi == 0

    def isEulerian(self):
        """ Return true iff graph has Eulerian walk or cycle """
        # technically, if it has an Eulerian walk
        return self.hasEulerianWalk() or self.hasEulerianCycle()

    def eulerianWalkOrCycle(self):
        """ Find and return sequence of nodes (represented by
            their k-1-mer labels) corresponding to Eulerian walk
            or cycle """
        assert self.isEulerian()
        g = self.G
        if self.hasEulerianWalk():
            g = g.copy()
            g.setdefault(self.tail, []).append(self.head)
        # graph g has an Eulerian cycle
        tour = []
        src = next(iter(g.keys())) # pick arbitrary starting node

        def __visit(n):
            while len(g[n]) > 0:
                dst = g[n].pop()
                __visit(dst)
            tour.append(n)

        __visit(src)
        tour = tour[::-1][:-1] # reverse and then take all but last node

        if self.hasEulerianWalk():
            # Adjust node list so that it starts at head and ends at tail
            sti = tour.index(self.head)
            tour = tour[sti:] + tour[:sti]

        # Return node list
        return list(map(str, tour))


def extract_contigs_greedy(de_bruijn_graph):
    """
    Extracts contigs from a de Bruijn graph using a greedy approach.

    Args:
        de_bruijn_graph: A DeBruijnGraph object.

    Returns:
        A list of contig strings.
    """
    contigs = []
    visited_nodes = set()

    for node in de_bruijn_graph.nodes.values():
        if node not in visited_nodes:
            contig = node.km1mer
            visited_nodes.add(node)
            current_node = node

            while current_node in de_bruijn_graph.G and \
                  de_bruijn_graph.G[current_node]:  # Check if the node has outgoing edges
                # Get the first unvisited neighbor
                next_node = next(iter(de_bruijn_graph.G[current_node].keys()))
                if next_node in visited_nodes:
                    break
                contig += next_node.km1mer[-1]
                visited_nodes.add(next_node)
                current_node = next_node

            contigs.append(contig)
    return contigs
