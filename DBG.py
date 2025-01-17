# import graphviz

# TODO:
# improve extract_contigs
# correct remove_tips and remove_bubbles

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

    def get_edge_weight(self, nodeL, nodeR):
        """ Return the weight of the edge from nodeL to nodeR. """
        return self.G.get(nodeL, {}).get(nodeR, 0)

    def get_ancestors(self, node):
        """ Return the ancestors of a node. """
        return list(self.reverse_G.get(node, {}).keys())

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
    
    def remove_tips(self, weight_threshold):
        """ Remove tips from the graph.
        
        Args:
            weight_threshold (int): The maximum weight for an edge to be considered a tip.
        """
        removed = True  # To keep track of whether we removed a tip in the last iteration

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

    def remove_bubbles(self):
        """Remove bubbles from the graph by identifying diverging and converging paths."""
        removed = True  # Track if bubbles are removed

        while removed:
            removed = False

            for node in list(self.nodes.values()):
                # Check if the node is a divergence point (more than one outgoing edge)
                if len(self.G.get(node, {})) > 1:
                    # Find all paths starting from this node
                    paths = self._find_paths_from_node(node)
                    if not paths:
                        continue

                    # Identify the convergence point and remove redundant paths
                    convergence_point = self._find_convergence_point(paths)
                    if convergence_point:
                        # Select the best path and remove others
                        self._resolve_bubble(paths, convergence_point)
                        removed = True

    def _find_paths_from_node(self, start_node, max_depth=100):
        """Find all paths starting from a node up to a maximum depth."""
        paths = []
        stack = [(start_node, [start_node])]  # (current_node, path_so_far)

        while stack:
            current_node, path = stack.pop()

            # Stop exploring if we exceed max depth
            if len(path) > max_depth:
                continue

            # Check for convergence point
            if len(path) > 1 and len(self.G.get(current_node, {})) == 0:
                paths.append(path)
                continue

            # Add neighbors to the stack
            for neighbor in self.G.get(current_node, {}):
                if neighbor not in path:  # Avoid cycles
                    stack.append((neighbor, path + [neighbor]))

        return paths

    def _find_convergence_point(self, paths):
        """Find the first node that all paths re-converge to."""
        if not paths:
            return None

        # Use sets to track nodes in each path
        sets = [set(path) for path in paths]

        # Find common nodes among all paths
        common_nodes = set.intersection(*sets)
        if not common_nodes:
            return None

        # Return the first common node in any path as the convergence point
        for node in paths[0]:
            if node in common_nodes:
                return node
        return None

    def _resolve_bubble(self, paths, convergence_point):
        """Remove redundant paths in a bubble."""
        # Evaluate paths and select the best one
        path_scores = [(path, self._evaluate_path(path)) for path in paths]
        path_scores.sort(key=lambda x: x[1], reverse=True)  # Sort by score (higher is better)

        best_path = path_scores[0][0]
        redundant_paths = [path for path, _ in path_scores[1:]]

        # Remove edges of redundant paths
        for path in redundant_paths:
            for i in range(len(path) - 1):
                node, next_node = path[i], path[i + 1]
                weight = self.G[node].pop(next_node, 0)
                if weight:
                    self.reverse_G[next_node].pop(node, 0)
                    node.nout -= weight
                    next_node.nin -= weight
                if not self.G[node]:
                    del self.G[node]
                if not self.reverse_G[next_node]:
                    del self.reverse_G[next_node]

    def _evaluate_path(self, path):
        """Evaluate a path based on cumulative edge weight."""
        score = 0
        for i in range(len(path) - 1):
            node, next_node = path[i], path[i + 1]
            score += self.G.get(node, {}).get(next_node, 0)
        return score


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


# class DeBruijnGraph2(DeBruijnGraph):
#     def to_dot(self, weights=False):
#         """ Return string with graphviz representation.  If 'weights'
#             is true, label edges corresponding to distinct k-1-mers
#             with weights, instead of drawing separate edges for
#             k-1-mer copies. """
#         g = graphviz.Digraph(comment='DeBruijn graph')
#         for node in iter(self.G.keys()):
#             g.node(node.km1mer, node.km1mer)
#         for src, dsts in iter(self.G.items()):
#             if weights:
#                 weightmap = {}
#                 if weights:
#                     for dst in dsts:
#                         weightmap[dst] = weightmap.get(dst, 0) + 1
#                 for dst, v in weightmap.items():
#                     g.edge(src.km1mer, dst.km1mer, label=str(v))
#             else:
#                 for dst in dsts:
#                     g.edge(src.km1mer, dst.km1mer)
#         return g
    
# def extract_contigs_greedy(de_bruijn_graph):
#     """
#     Extracts contigs from a de Bruijn graph using a greedy approach.

#     Args:
#         de_bruijn_graph: A DeBruijnGraph object.

#     Returns:
#         A list of contig strings.
#     """
#     contigs = []
#     visited_nodes = set()

#     for node in de_bruijn_graph.nodes.values():
#         if node not in visited_nodes:
#             contig = node.km1mer
#             visited_nodes.add(node)
#             current_node = node

#             while current_node in de_bruijn_graph.G and \
#                   de_bruijn_graph.G[current_node] and \
#                   de_bruijn_graph.G[current_node][0] not in visited_nodes:
#                 next_node = de_bruijn_graph.G[current_node][0]
#                 contig += next_node.km1mer[-1]
#                 visited_nodes.add(next_node)
#                 current_node = next_node
#             contigs.append(contig)
#     return contigs

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

def remove_tips(de_bruijn_graph, tip_length_threshold):
    """
    Remove tips from a De Bruijn graph.

    Args:
        de_bruijn_graph: A DeBruijnGraph object.
        tip_length_threshold: Maximum length for a tip to be removed.

    Returns:
        Modified DeBruijnGraph object.
    """
    to_remove = set()
    for node in de_bruijn_graph.nodes.values():
        if node.nout == 0 or node.nin == 0:  # Dead-end nodes
            path_length = 0
            current_node = node
            while current_node and path_length <= tip_length_threshold:
                if current_node in to_remove:
                    break
                to_remove.add(current_node)
                if current_node.nout > 0:
                    current_node = de_bruijn_graph.G[current_node][0]
                else:
                    current_node = None
                path_length += 1
            if path_length > tip_length_threshold:
                to_remove.difference_update(to_remove)

    for node in to_remove:
        try:
            del de_bruijn_graph.G[node]
            del de_bruijn_graph.nodes[node.km1mer]
        except KeyError:
            pass
    return de_bruijn_graph


def remove_bubbles(reads, de_bruijn_graph, k, similarity_threshold):
    """
    Remove bubbles from a De Bruijn graph.

    Args:
        de_bruijn_graph: A DeBruijnGraph object.
        similarity_threshold: Maximum allowed distance between paths to be considered a bubble.

    Returns:
        Modified DeBruijnGraph object.
    """
    kmerhist = kmerHist(reads, k)
    for node in de_bruijn_graph.nodes.values():
        if len(de_bruijn_graph.G[node]) > 1:  # Multiple outgoing paths
            paths = []
            for dst in de_bruijn_graph.G[node]:
                current_path = dst.km1mer
                current_node = dst
                while current_node and current_node.nout == 1:
                    current_node = de_bruijn_graph.G[current_node][0]
                    current_path += current_node.km1mer[-1]
                paths.append(current_path)

            for i in range(len(paths)):
                for j in range(i + 1, len(paths)):
                    if hamming_distance(paths[i], paths[j]) <= similarity_threshold:
                        # Remove the less frequent path
                        freq_i = sum([kmerhist[kmer] for kmer in paths[i]])
                        freq_j = sum([kmerhist[kmer] for kmer in paths[j]])
                        to_remove = paths[i] if freq_i < freq_j else paths[j]
                        for kmer in to_remove:
                            try:
                                del de_bruijn_graph.nodes[kmer]
                            except KeyError:
                                pass
    return de_bruijn_graph


def hamming_distance(s1, s2):
    """
    Calculate the Hamming distance between two strings.
    """
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))
