# What?

A place where you can draw graphs and the site will tell you
whether it's a counter-example to the
[Woodhall conjecture](http://www.openproblemgarden.org/op/woodalls_conjecture).

## The code.

Unreadable.

### index.html

There's a small bit of code here that maps from the vis.js graph structure and events to the internal structure.
Only that and nothing more. (Well, and the trivial HTML rendering of the data for dicuts and disjoint dijoins.)

### graph-util.js

- A circular linked-list queue
- An array-based priority queue (heap)
- Tarjan's SCC finding algorithm (this partitions the graph into equivalence classes)
- A conversion from an assignment of vertex to SCC to two DAGs -- one with edges weighted by how many edges
  were in the original graph and one weighted by the edges themselves.
- A function to get dicuts (really sort of a module) that contains:
  - A BFS.
  - An implementation of Ford Fulkerson
- A function to get disjoint dijoins (assuming the conjecture holds and a greedy sort of algorithm is fine).
  Also sort of a module containing:
  - A dijoin object that tracks how many SCCs it has yet to merge (and where they are), the edges it's reversed,
    and figures out an edge (not in the pool mentioned below) that it can use to merge some SCCs.
  - A universal pool of edges for all the dijoin objects in question so that they stay disjoint.

## Stuff I wanna fix.

- Can I present this better (ok, but ugh CSS)?
  - Maybe showing the DAG of SCCs can be handy?
  - Coloring the main graph and the cut edges?
- Look into incremental implementations of this to really save time.
  - (I think incrementally knowing SCCs is a pain and so have yet to try.)
- Clean up the code. (Because who doesn't want that?)