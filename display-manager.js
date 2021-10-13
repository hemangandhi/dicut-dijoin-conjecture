var sanitizeNodeToNodeMap = function(sanitizer, map) {
    var clean = new Map();
    for (const [s, t] of map) {
	clean.set(sanitizer(s), sanitizer(t));
    }
    return clean;
}

var sanitizeEdgePair = function(sanitizer, pair) {
    return [sanitizer(pair[0]), sanitizer(pair[1])];
};

var sanitizeSetOfEdges = function (sanitizer, edge_set) {
    return edge_set.map(function(e) { return sanitizeEdgePair(sanitizer, e); });
};

var DisplayManager = function(easel, notifier, dicut_display, dijoins_display, debug_checkbox) {
    var nodes = new vis.DataSet([]);
    var edges = new vis.DataSet([]); 
    var data = {
	nodes: nodes,
	edges: edges
    };
    var graph = new Map();
    var tmp_expected_node_id = null;
    var node_id_map = new Map();
    // Needed to map from edge IDs to the actual
    // data for the edge.
    var edge_map = new Map();
    var tmp_expected_edge_info = null;
    // setup for node naming
    nodes.on('add', function(e, props, s) {
	if (tmp_expected_node_id === null) return;
	node_id_map.set(props.items[0], tmp_expected_node_id);
	tmp_expected_node_id = null;
    });
    edges.on('add', function(e, props, s) {
	if (tmp_expected_edge_info === null) return;
	edge_map.set(props.items[0], tmp_expected_edge_info);
	tmp_expected_edge_info = null;
    });
    var addEdgeToGraph = function(edge) {
	var src = edge.from;
	var dest = edge.to;
	if (!graph.has(src)) {
	    graph.set(src, new Set());
	}
	graph.get(src).add(dest);
	tmp_expected_edge_info = edge;
    };
    var removeEdgeByVisId = function(vis_id) {
	var edge = edge_map.get(vis_id);
	if (!graph.has(edge.from)) return;
	var out_edges = graph.get(edge.from);
	out_edges.delete(edge.to);
	if (out_edges.size == 0) {
	    graph.delete(edge.from);
	}
    };

    this.makeNodeSanitizer = function(sccs) {
	var reversed_map = new Map();
	for (const [node_id, scc_id] of sccs) {
	    if (!reversed_map.has(scc_id)) {
		reversed_map.set(scc_id, [node_id]);
	    } else {
		reversed_map.get(scc_id).push(node_id);
	    }
	}
	return function(node_id) {
	    if (!debug_checkbox.checked) return "(debugging disabled)";
	    if (!reversed_map.has(node_id)) return node_id + " (not found)";
	    return "the SCC containing nodes {" + reversed_map.get(node_id)
		.map(function (i) { return node_id_map.get(i); }).join(", ") + "}";
	};
    };

    var held_dicut = [];
    this.setDicutState = function(dicut) { held_dicut = dicut; };
    this.getDicutState = function() { return held_dicut; };
    var held_dijoins = [];
    this.setDijoinState = function(dijoins) { held_dijoins = dijoins };
    var state_stack = [];
    this.pushState = function(state) { if (debug_checkbox.checked) state_stack.push(state); };
    this.popState = function() { if (debug_checkbox.checked) return state_stack.pop(); };
    var held_error = false;
    this.setError = function(error) {
	if (debug_checkbox.checked) {
	    console.log(state_stack.join("\n"));
	    console.log("Error found in the above: " + error);
	}
	held_error = error;
    };
    this.getError = function() { return held_error; };

    var _this = this;
    // Conjecture state.
    var updateConjectureState = function(graph) {
	getGraphConjectureInfo(graph, _this);
	notifier.innerHTML = "";
	dicut_display.innerHTML = "";
	dijoins_display.innerHTML = "";
	var dicut_list = document.createElement('ol');
	if (held_error) {
	    notifier.innerText = held_error;
	} else if (held_dicut.length == held_dijoins.length) {
	    notifier.innerText = "Barring bugs, this does not look like a counter-example. Please verify the examples below.";
	} else {
	    notifier.innerText = "Barring bugs, this looks like a counter-example. Please verify the examples below.";
	}
	for (const [src, dest] of held_dicut) {
	    var li = document.createElement('li');
	    li.innerText = "The edge from " + node_id_map.get(src) + " to " + node_id_map.get(dest);
	    dicut_list.appendChild(li);
	}
	dicut_display.appendChild(dicut_list);
	var dijoin_list = document.createElement('ol');
	for (const dijoin of held_dijoins) {
	    var li = document.createElement('li');
	    var span = document.createElement('span');
	    span.innerText = "Dijoin that adds edges: ";
	    li.appendChild(span);
	    var per_dijoin = document.createElement('ol');
	    for (const [src, target] of dijoin) {
		var li2 = document.createElement('li');
		li2.innerText = "From " + node_id_map.get(target) + " to " + node_id_map.get(src);
		per_dijoin.appendChild(li2);
	    }
	    li.appendChild(per_dijoin);
	    dijoin_list.appendChild(li);
	}
	dijoins_display.appendChild(dijoin_list);
    };

    var node_ctr = 0;
    var options = {
	manipulation: {
	    enabled: true,
	    addNode: function(node, done) {
		node.label = '' + node_ctr;
		node_ctr++;
		tmp_expected_node_id = node.label;
		done(node);
	    },
	    addEdge: function(edge, done) {
		addEdgeToGraph(edge);
		updateConjectureState(graph);
		done(edge);
	    },
	    editEdge: function(edge, done) {
		removeEdgeByVisId(edge.id);
		addEdgeToGraph(edge);
		updateConjectureState(graph);
		done(edge);
	    },
	    deleteEdge: function(edge, done) {
		removeEdgeByVisId(edge.edges[0]);
		updateConjectureState(graph);
		done(edge);
	    },
	    deleteNode: function(selection, done) {
		for(const node_id of selection.nodes) {
		    node_id_map.delete(node_id);
		}
		for(const edge_id of selection.edges) {
		    removeEdgeByVisId(edge_id);
		}
		updateConjectureState(graph);
		done(selection);
	    }
	},
	edges: {
	    arrows: 'to'
	},
    };

    var network = new vis.Network(easel, data, options);
};
