# C++ practice repository

## Intro

This is some coursework (happens to be mine: Janos Brezniczky), please do not 
plagiarise, but feel free to check it out otherwise.

Before, but also after reviewing others' solutions, I had some thoughts -- so
the hope is maybe, if you are a fellow coursetaker, it can give you some ideas 
as well, or reassure you about concepts you didn't have time to try out.

## About the solution

This is a Monte-Carlo style estimation of the average shortest path length in
graphs of certain properties.

The method:

* a number of random graphs with certain properties are generated
* for each graph, shortest path lengths from a certain node to all of the other 
  nodes are calculated
* a mean shortest path is calculated for each graph, and printed
* a mean of the above means is calculated and printed

This is repeated with graphs of 0.2 and 0.4 density (density: probability of two
nodes being connected).
