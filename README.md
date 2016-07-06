# C++ practice repository

## Intro

This is some coursework (happens to be mine: Janos Brezniczky), please do not 
plagiarise, but otherwise if you are a classmate, feel free to check it out.

Before, but also after reviewing others' solutions, I had some thoughts -- so
the hope is maybe, if you are a fellow coursetaker, it can give you some ideas 
as well, or reassure you about concepts you didn't have time to try out.

## About the solution

Homework2

This is a Monte-Carlo style estimation of the average shortest path length in
graphs of certain properties.

The method:

* a number of random graphs with certain properties are generated
* for each graph, shortest path lengths from a certain node to all of the other 
  nodes are calculated
* a mean shortest path is calculated for each graph, and printed
* a mean of the above means is calculated and printed

This is performed with graphs of 0.2 density, then of 0.4 (density: probability 
of two nodes being connected).

Homework3

This is the previous thing slightly extended to also implement a minimum
spanning tree algorithm, namely Prim's, and runs that with one certain input
(at hard-coded path at the time writing, so if it doesn't run, worth a check).


