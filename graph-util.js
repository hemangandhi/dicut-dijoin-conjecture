// Circular linked list as a queue.
var Queue = function() {
    var Node = function(data, next) {
	this.data = data;
	this.next = next;
    };

    var tail = null;
    this.enqueue = function(datum) {
	if (tail == null) {
	    tail = new Node(datum, null);
	    tail.next = tail;
	} else {
	    tail.next = new Node(datum, tail.next);
	}
    };
    this.dequeue = function() {
	var ret = tail.next.data;
	if (tail.next != tail) {
	    tail.next = tail.next.next;
	} else {
	    tail = null;
	}
	return ret;
    };
    this.empty = function() { return tail == null; };
};

// I'm such a nerd.
var Heap = function (cmp) {
    var list = [];
    this.push = function(datum) {
	list.push(datum);
	for (var i = list.length - 1; i > 0; i = Math.floor(i/2)) {
	    if (cmp(list[i], list[Math.floor(i/2)])) {
		var tmp = list[i];
		list[i] = list[Math.floor(i/2)];
		list[Math.floor(i/2)] = tmp;
	    }
	}
    };
    this.pop = function () {
	if (list.length == 0) return null;
	if (list.length == 1) return list.pop();
	var max = list[0];
	list[0] = list.pop();
	var i, l, r;
	for (i = 0, l = 1, r = 2; r < list.length; l = i * 2 + 1, r = i * 2 + 2) {
	    var new_i = (cmp(list[l], list[r])) ? l : r;
	    if (!cmp(list[i], list[new_i])) break;
	    var tmp = list[i];
	    list[i] = list[new_i];
	    list[new_i] = tmp;
	    i = new_i;
	}
	// Consider popping off of a size 3 heap (like [3, 2, 1] as a min heap).
	if (l < list.length && cmp(list[i], list[l])) {
	    var tmp = list[i];
	    list[i] = list[l];
	    list[l] = tmp;
	}
	return max;
    };
    this.empty = function() {
	return list.length == 0;
    };
};

// Graphs.
//
// They are either represented as a Map<K, Set<K>> if there's no weight or
// Map<K, Map<K, W>> if it is weighted (K is the key type, usually strings, and
// W is the weight type, usually numeric).

var edgeSet = function(graph) {
    var sources = new Set();
    for (const [source, dests] of graph) {
	sources.add(source);
	for(const dest of dests.keys()) {
	    sources.add(dest);
	}
    }
    return sources;
};

var getARootVertex = function(dag) {
    var maybe_roots = new Set(dag.keys());
    for (const [r, neighbours] of dag) {
	for (const neighbour of neighbours.keys()) {
	    maybe_roots.delete(neighbour);
	}
    }
    return maybe_roots.keys().next().value;
};

var isWeaklyConnected = function(graph) {
    var reversed_edges = new Map();
    var unvisiteds = new Set();
    for (const [vertex, neighbors] of graph) {
	unvisiteds.add(vertex);
	for (const neighbor of neighbors) {
	    unvisiteds.add(neighbor);
	    if (reversed_edges.has(neighbor)) {
		reversed_edges.get(neighbor).add(vertex);
	    } else {
		var r_neighs = new Set();
		r_neighs.add(vertex);
		reversed_edges.set(neighbor, r_neighs);
	    }
	}
    }

    var recurse = function(v) {
	if (graph.has(v)) {
	    for (const neighbor of graph.get(v)) {
		if (unvisiteds.delete(neighbor)) {
		    recurse(neighbor);
		}
	    }
	}
	if (reversed_edges.has(v)) {
	    for (const neighbor of reversed_edges.get(v)) {
		if (unvisiteds.delete(neighbor)) {
		    recurse(neighbor);
		}
	    }
	}
    }

    var some_vertex = unvisiteds.values().next().value;
    unvisiteds.delete(some_vertex);
    recurse(some_vertex);
    return unvisiteds.size == 0;
};

// https://en.wikipedia.org/wiki/Tarjan%27s_strongly_connected_components_algorithm
// Returns a Map<K, int> where each int identifies a component.
var getSCCs = function(graph) {
    var index = 0;
    var stack = [];
    var sccs = new Map();
    var scc_id = 0;
    var extra_info = new Map();

    var recurse = function(vertex) {
	extra_info.set(vertex, {
	    index: index,
	    lowlink: index,
	    on_stack: true
	});
	index++;
	stack.push(vertex);
	if (graph.has(vertex)) {
	    for(const neighbor of graph.get(vertex)) {
		if (!extra_info.has(neighbor)) {
		    recurse(neighbor);
		    extra_info.get(vertex).lowlink = Math.min(
			extra_info.get(vertex).lowlink,
			extra_info.get(neighbor).lowlink);
		}else if (extra_info.get(neighbor).on_stack) {
		    extra_info.get(vertex).lowlink = Math.min(
			extra_info.get(vertex).lowlink,
			extra_info.get(neighbor).lowlink);
		}
	    } // for each neighbor
	}

	var curr_info = extra_info.get(vertex);
	if (curr_info.index != curr_info.lowlink) return;
	var memb;
	do {
	    memb = stack.pop();
	    sccs.set(memb, scc_id);
	    extra_info.get(memb).on_stack = false;
	} while(memb != vertex);
	scc_id++;
    }; // function recurse(vertex)
    
    for (const vertex of graph.keys()) {
	if (!extra_info.has(vertex)) {
	    recurse(vertex);
	    scc_id++;
	}
    }
    return sccs;
}; // function getSCCs(graph)

// Returns two representations of the SCCs (as a pair):
//  1) a DAG with just edge weights corresponding to the number of
//     edges between the components.
//  2) a DAG where the edge weights are the set of edges in the original
//     graph (where the edge is a [source, target] pair).
var getDAGsOfSCCs = function(graph, sccs) {
    var scc_graph = new Map();
    var scc_edges = new Map();
    for (const [vertex, neighbours] of graph) {
	var vertex_scc = sccs.get(vertex);
	if (!scc_graph.has(vertex_scc)) {
	    scc_graph.set(vertex_scc, new Map());
	    scc_edges.set(vertex_scc, new Map());
	}
	var scc_neighbors = scc_graph.get(vertex_scc);
	var neighbor_edges = scc_edges.get(vertex_scc);
	for (const neighbor of neighbours) {
	    var neighbor_scc = sccs.get(neighbor);
	    if (neighbor_scc == vertex_scc) continue;
	    if (!scc_neighbors.has(neighbor_scc)) {
		neighbor_edges.set(neighbor_scc, [[vertex, neighbor]]);
		scc_neighbors.set(neighbor_scc, 1);
	    } else {
		scc_neighbors.set(neighbor_scc, scc_neighbors.get(neighbor_scc) + 1);
		neighbor_edges.get(neighbor_scc).push([vertex, neighbor]);
	    }
	}
    }
    return [scc_graph, scc_edges];
};

// SCC graph is the weighted DAG of SCCs (Map<K, Map<K, int>>).
// SCC edges is the same DAG, but each edge knows what edges it corresponds to in the
// original graph (Map<K1, Map<K1, List<Pair<K2, K2>>) where K1 is the SCC id and
// K2 is the original vertex type).
//
// Returns both the size of the dicut and (a guess of) the edges in the original graph
// that would be cut.
var dicut = function(scc_graph, scc_edges) {
    // Converts a dicut expressed as cuts in the SCC DAG
    var augmentEdgeSet = function(scc_edges, edge_set) {
	var aug = [];
	for (const [source_scc, dest_scc] of edge_set) {
	    for (const og_edge of scc_edges.get(source_scc).get(dest_scc)) {
		aug.push(og_edge);
	    }
	}
	return aug;
    };
    
    // If a path from source to target exists, this will return a tree representing the path.
    // That is, a graph with only reversals of edges that were traversed by the BFS.
    // If no path is found, this function returns null.
    var bfs = function (graph, source, target) {
	var parents = new Map();
	var visited = new Set();
	var q = new Queue();
	q.enqueue(source);
	visited.add(source);
	while(!q.empty()) {
	    var nd = q.dequeue();
	    for (const [neighbor, w] of graph.get(nd)) {
		if (visited.has(neighbor) || (w <= 0 && w != -1)) {
		    continue;
		}
		q.enqueue(neighbor);
		visited.add(neighbor);
		parents.set(neighbor, nd);
	    }
	} // while (!q.empty())
	if (visited.has(target) && parents.size > 0) return parents;
	return null;
    };

    // NOTE: this is Ford Fulkerson++ (persay) -- it adds infnite capacity reversed
    // edges and then computes the max flow.
    // https://en.wikipedia.org/wiki/Ford%E2%80%93Fulkerson_algorithm
    //
    // Returns a pair of the weight and (a guess of) the cut edges (as a List<Pair<K, K>>).
    var fordFulkerson = function(scc_graph, source, target) {
	// Some jank (because I hate floating point "inf"): -1 is infinity.
	var path_flow = new Map();
	for (const [vertex, neighbors] of scc_graph) {
	    if (!path_flow.has(vertex)) path_flow.set(vertex, new Map());
	    path_flow_map = path_flow.get(vertex);
	    for (const [neighbor, weight] of neighbors) {
		path_flow_map.set(neighbor, weight);
		if (!path_flow.has(neighbor)) {
		    var ns = new Map();
		    ns.set(vertex, -1);
		    path_flow.set(neighbor, ns);
		} else {
		    path_flow.get(neighbor).set(vertex, -1);
		}
	    }
	}

	var max_flow = 0;
	var cut_edges = [];
	var i = 0;
	for (var path = bfs(path_flow, source, target); path != null && i < 3;
	     path = bfs(path_flow, source, target), i++) {
	    var min_wt = -1;
	    var min_edge = null;
	    for (var s = target; s != source && path.has(s); s = path.get(s)) {
		var weight = path_flow.get(path.get(s)).get(s);
		if (min_wt == -1 || (weight < min_wt && weight > 0)) {
		    min_wt = weight;
		    min_edge = [path.get(s), s];
		}
	    }
	    // The only flow is on reversed edges -- let's not cut those.
	    // And, yeah, fuck it, we'll just abandon our search at that point.
	    if (min_wt <= 0) {
		// OK, this is an error, but it's not one I
		// want a user to understand.
		console.log("Path only passes reversed edges");
		console.log(path);
		break;
	    }
	    max_flow += min_wt;
	    cut_edges.push(min_edge);
	    for (var s = target; s != source && path.has(s); s = path.get(s)) {
		var p = path.get(s);
		var weight = path_flow.get(p).get(s);
		if (weight > 0) {
		    path_flow.get(p).set(s, weight - min_wt);
		}
		// rev = reversed along the s-t path (not necessarily reversed
		// relative to the original DAG).
		var rev_weight = path_flow.get(s).get(p);
		if (rev_weight > 0) {
		    path_flow.get(s).set(p, weight + min_wt);
		}
	    }
	} // for each BFS path
	return [max_flow, cut_edges];
    };

    var min_dicut = -1;
    var cut_edges = [];
    var some_source = getARootVertex(scc_graph);
    for (some_dest of edgeSet(scc_graph).keys()) {
	if (some_dest == some_source) continue;
	var flow = fordFulkerson(scc_graph, some_source, some_dest);
	if (min_dicut == -1 || flow[0] < min_dicut) {
	    min_dicut = flow[0];
	    cut_edges = flow[1];
	}
    }
    return [min_dicut, augmentEdgeSet(scc_edges, cut_edges)];
};

// Conjecture: the Woodhall conjecture holds and a set of disjoint dijoins of
// the required size can be computed greedily by expanding the set of edges
// the furthest away from being a dijoin by adding the set that connects it
// the most (of edges not already in another dijoin) until all the edges in
// the dicut have expanded to form the disjoint dijoins.
var getDisjointDijoins = function(graph, scc_dag, sccs, dicut_edge_set) {
    // Set of edges used by some dijoin in this pool
    var pool_allocation = new Set();
    var addEdgeToPool = function (edge) {
	var og_edge_id = edge[0] + edge[1];
	if (pool_allocation.has(og_edge_id)) {
	    return false;
	}
	pool_allocation.add(og_edge_id);
	return true;
    };

    // EXPECTATION: all edges added to a Dijoin are of type Pair<string, string>.
    var GrowingDijoin = function(init_edge) {
	var edge_array = [];
	// TODO: is there a faster shallow copy?
	var my_sccs = new Map();
	for (const [vertex, scc_id] of sccs) {
	    my_sccs.set(vertex, scc_id);
	}
	var neediness = -1;

	// Returns what would happen if src and dest were connected
	// in my_scc.
	var supposeConnected = function(og_edge, curr_dag) {
	    var src_scc_class = my_sccs.get(og_edge[0]);
	    var dest_scc_class = my_sccs.get(og_edge[1]);
	    var weightless_dag = new Map();
	    for (const [vertex, neighbors] of curr_dag) {
		var nbd = new Set();
		for (const neighbor of neighbors.keys()) {
		    nbd.add(neighbor);
		}
		weightless_dag.set(vertex, nbd);
	    }
	    if (weightless_dag.has(dest_scc_class)) {
		weightless_dag.get(dest_scc_class).add(src_scc_class);
	    } else {
		var n = new Set();
		n.add(src_scc_class);
		weightless_dag.set(dest_scc_class, n);
	    }
	    var new_sccs = getSCCs(weightless_dag);
	    // If you pretend that my_sccs: original vertex set -> old SCC IDs
	    // and new_sccs : old SCC IDs -> new SCC IDs, this is new_sccs \circ my_sccs
	    var composition = new Map();
	    var max_id = -1;
	    for (const [og_vertex, old_scc_id] of my_sccs) {
		var new_scc_id = new_sccs.get(old_scc_id);
		composition.set(og_vertex, new_scc_id);
		if (new_scc_id > max_id) {
		    max_id = new_scc_id;
		}
	    }
	    return {
		neediness: max_id,
		sccs: composition,
		edge: og_edge
	    };
	};
	
	var addEdge = function(edge) {
	    if (!addEdgeToPool(edge)) return false;
	    // NOTE: if they're equal, I don't know what we're doing, but
	    // oh well.
	    var src_scc_class = my_sccs.get(edge[0]);
	    var dest_scc_class = my_sccs.get(edge[1]);
	    if (src_scc_class != dest_scc_class) {
	        var curr_dag = getDAGsOfSCCs(graph, my_sccs)[0];
		var new_state = supposeConnected(edge, curr_dag);
		neediness = new_state.neediness;
		my_sccs = new_state.sccs;
	    }
	    edge_array.push(edge);
	    return true;
	};
	this.error = "";
	if (!addEdge(init_edge)) {
	    this.error = "Could not add dicut edge " + init_edge;
	}

	this.is_dijoin = function() {
	    return neediness == 0;
	};
	this.neediness = function() {
	    return neediness;
	};
	this.edges = function () {
	    return edge_array;
	};
	this.debugLog = function () {
	    console.log({
		edges: this.edges(),
		neediness: this.neediness(),
		error: this.error
	    });
	    console.log(this.edges().join(", "));
	}; 

	this.grow = function() {
	    if (neediness == 0) return;
	    var [scc_dag, grouped_edges] = getDAGsOfSCCs(graph, my_sccs);
	    var curr_min = null;
	    for (const [src, src_n] of grouped_edges) {
		for (const [dest, og_edges] of src_n) {
		    for (const og_edge of og_edges) {
			var og_edge_id = og_edge[0] + og_edge[1];
			if (pool_allocation.has(og_edge_id)) continue;
			var supposition = supposeConnected(og_edge, scc_dag);
			if (curr_min == null || supposition.neediness < curr_min.neediness) {
			    curr_min = supposition;
			}
			// By grouping into the edges in the SCC DAG, only checking
			// one feasible edge is (should be) enough.
			// Really, this is probably a crux of my conjecture:
			// if my earlier conjecture is correct, it's good enough to break here.
			break;
		    }
		}
	    }
	    if (curr_min == null) {
		this.error = "No feasible edges to grow along";
		return;
	    } else if (!addEdgeToPool(curr_min.edge)) {
		this.error = "Edge " + curr_min.edge + " pooled at a weird time";
		return;
	    }
	    edge_array.push(curr_min.edge);
	    neediness = curr_min.neediness;
	    my_sccs = curr_min.sccs;
	};
    };
    
    var growing_dijoins = new Heap(function (l, r){
	return l.neediness() > r.neediness();
    });
    var done_dijoins = [];
    for (const cut_edge of dicut_edge_set) {
	var dj = new GrowingDijoin(cut_edge);
	if (dj.error.length > 0) {
	    return [false, dj.error];
	}
	if (dj.is_dijoin()) {
	    done_dijoins.push(dj.edges());
	} else {
	    growing_dijoins.push(dj);
	}
    }
    while (!growing_dijoins.empty()) {
	var dj = growing_dijoins.pop();
	dj.grow();
	if (dj.error.length > 0) {
	    return [false, dj.error];
	}
	if (dj.is_dijoin()) {
	    var edges = dj.edges();
	    done_dijoins.push(edges);
	} else {
	    growing_dijoins.push(dj);
	}
    }
    return [true, done_dijoins];
};

var getGraphConjectureInfo = function(graph) {
    if (graph.size == 0) return {
	is_error: true,
	error_msg: "You seem to have emptied the graph out. Please add nodes and edges to verify the conjecture."
    };
    if (!isWeaklyConnected(graph)) return {
	is_error: true,
	error_msg: "Graph isn't connected yet. Please add more edges to weakly connect the graph."
    };
    var sccs = getSCCs(graph);
    if ((new Set(sccs.values())).size == 1) return {
	is_error: true,
	error_msg: "The graph seems to be an SCC, so the conjecture is trivially true."
    };
    var [dag, edges] = getDAGsOfSCCs(graph, sccs);
    var [dicut_size, dicut_edges] = dicut(dag, edges);
    var [is_ok, val] = getDisjointDijoins(graph, dag, sccs, dicut_edges);
    if (!is_ok) {
	return {
	    is_error: true,
	    error_msg: val
	};
    }
    return {
	is_error: false,
	dicut_edges: dicut_edges,
	disjoint_dijoins: val
    };
};
