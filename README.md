# C++ practice repository

## Intro

This is some coursework (happens to be mine: Janos Brezniczky), please do not 
plagiarise, but feel free to check it out otherwise.
While reviewing others' solutions, I had some thoughts in retrospect, maybe it
can give some ideas for you, too.

# About the solution

This is a variant of Dijkstra's shortest path algorithm, which leads to an 
empiric (Monte-Carlo style) estimation of the mean shortest path between two 
nodes in a graph. A more proper solution could benefit from a binomial heap, 
although with small matrices, the classes in the module should work fine.

Varying things:
* sides, between 1 and 10
* connectivity, namely the density of the graph is set at 0.2 in the first round,
  then at 0.4

