/*
  Implementation of Dijsktra's algorithm for shortest path in a graph, with a
  focus on estimating average path length between nodes.

  Main classes:

  Matrix
  Graph
  ShortestPathFinder

  The classes have the "<<" operator defined for verification purposes.

  Future improvements:
  - define (x, y) operator on the matrix class

  - the matrix could be triangular, since it is symmetric - memory consumption
    could be reduced

  - create generic mean or sum calculation solution
*/

#include <iostream>
#include <vector>
#include <random>
#include <map>

using namespace std;


/* Infinite value constant. */
const double INF = numeric_limits<double>::infinity();


class Matrix {
  /* Class representing a matrix of floating point values. */
private:
  double* values;
  int rows;
  int cols;

  /* Quickly retrieves the index of the memory slot for a given cell. */
  inline int get_index(int row, int col) {
    return(row * cols + col);
}

public:
  /* Creates and a matrix of the given size.
     The values are allocated but undefined. */
  Matrix(int rows, int cols) {
    this->rows = rows;
    this->cols = cols;
    values = new double[rows * cols];
  }

  ~Matrix() {
    delete(values);
  }

  /* Gets the value of a given cell. */
  double inline get(int row, int col) {
    return(values[get_index(row, col)]);
  }

  /* Sets the value of a given cell. */
  void inline set(int row, int col, const double value) {
    values[get_index(row, col)] = value;
  }

  /* Tells how many rows the matrix has. */
  inline int get_rows() {
    return(rows);
  }

  /* Tells how many columns the matrix has. */
  inline int get_cols() {
    return(cols);
  }
};


class Graph {
  /* A simple graph class allowing directed and undirected weighted graphs to be
     specified, and for random graph generation.

     Nodes are represented by integer values ranging from 0..(n - 1), where n is
     the number of nodes in the graph.
  */
private:
  /* Adjacency (aka. connectivity) matrix of the graph. Weights represent
     distances. "inf" weights are to be used between the non-adjacent nodes. */
  Matrix* adjacency;

  /* Initializes a symmetric adjacency matrix with random weights. */
  void random_init(double density, double min_dist, double max_dist, int seed) {

    default_random_engine generator(seed);
    uniform_real_distribution<double> dist_distribution(min_dist, max_dist);
    uniform_real_distribution<double> is_connected_dist(0.0, 1.0);

    for(int i = 0; i < adjacency->get_rows(); ++i) {
      for(int j = i + 1; j < adjacency->get_cols(); ++j) {
        // the adjacency matrix of an undirected graph is symmetric, set the
        // values equally for both directions
        if (is_connected_dist(generator) < density) {
          double dist = dist_distribution(generator);
          adjacency->set(i, j, dist);
          adjacency->set(j, i, dist);
        } else {
          adjacency->set(i, j, INF);
          adjacency->set(j, i, INF);
        }
      }
    }
  }

public:
  /* Creates the graph with random density and distances as described by the
     given parameters. */
  Graph(int nodes = 50, double density = 0.1,
        double min_dist = 1.0, double max_dist = 10.0, int seed = 0) {

    adjacency = new Matrix(nodes, nodes);;
    random_init(density, min_dist, max_dist, seed);
  }

  /* Destrutor - frees up allocated memory. */
  ~Graph() {
    delete adjacency;
  }

  /* Tells the number of nodes in the graph. */
  int get_node_count() {
    return(adjacency->get_rows());
  }

  /* Returns a reference to the adjacency matrix allowing external
     modifications. */
  Matrix& get_adjacency_matrix() {
    return(*adjacency);
  }

  /* Returns the length of a given edge. */
  double get_edge_length(int node1, int node2) {
    return(adjacency->get(node1, node2));
  }

  /* Returns all the edges originating from a given node. */
  vector<pair<int, double>> get_edges(int node) {

    vector<pair<int, double>> neighbours = vector<pair<int, double>>();

    for(int i = 0; i < adjacency->get_cols(); i++) {
      double dist = adjacency->get(node, i);
      if (dist < INF) {
        neighbours.push_back(pair<int, double>(i, dist));
      }
    }

    return(neighbours);
  }
};


template<class Key, class Priority> class SimplePriorityQueue {
  /* A map-based simple priority queue implementation.

     Scales poorly for larger input, but is a useful placeholder while it needs
     no replacement. A more efficient (but also involving) version can be
     created then. Suffices for the current input set.
  */
private:
  map<Key, Priority> priority_by_key;
public:
  /* Adds an item with the given priority. */
  void push(const Key &key, const Priority &priority) {
    priority_by_key.emplace(key, priority);
  }

  /* Removes the item with the given key. */
  bool erase(const Key key) {
    return(priority_by_key.erase(key) > 0);
  }

  /* Removes and retrieves the highest priority item
     (i.e. lowest priority value). */
  bool pop_min(Key &key, Priority &priority) {
    if (priority_by_key.empty())
      return(false);

    typename map<Key, Priority>::iterator it = priority_by_key.begin();

    key = it->first;
    priority = it->second;

    while(it != priority_by_key.end()) {
      if (it->second < priority) {
        key = it->first;
        priority = it->second;
      }
      ++it;
    }

    if (!erase(key)) {
      cout << "could not delete key " << key << endl;
      erase(key);
    }

    return(true);
  }

  /* Changes the priority associated with the item as specified. */
  void change_priority(const Key& key, const Priority& new_priority) {
    erase(key);
    push(key, new_priority);
  }

  /* Tells if the queue is empty. */
  bool empty() {
    return(priority_by_key.empty());
  }

  /* Allows priority queues to be printed to output streams in a human readable
     way. */
  ostream& operator<<(ostream& out) {
    typename map<Key, Priority>::iterator it = priority_by_key.begin();

    while(it != priority_by_key.end()) {
      cout << it->first << ":" << it -> second << endl;
      ++it;
    }
  }
};


class ShortestPathFinder {
  /* Class implementing Dijsktra's shortest path algorithm, with focusing only
     on the average path lengths from a given node. */
private:
  Graph* graph;
public:
  /* Creates the instance, preparing for operations on the given graph. */
  ShortestPathFinder(Graph &g) {
    graph = &g;
  }

  /* Returns all the lengths of shortest paths between a given node and all
     the nodes, including itself.

     A simplified version of the full Dijsktra algorithm, based on the
     pseudo-code at https://en.wikipedia.org/wiki/Dijkstra%27s_algorithm */
  vector<double> get_shortest_dists(int from_node) {
    int node_count = graph->get_node_count();
    vector<double> dist = vector<double>();

    SimplePriorityQueue<int, double> q = SimplePriorityQueue<int, double>();

    for(int i = 0; i < node_count; ++i) {
      if (i != from_node)
        dist.push_back(INF);
      else
        dist.push_back(0);
      q.push(i, dist[i]);
    }

    dist[from_node] = 0;

    while (!q.empty()) {
      int act_node;
      double act_dist;
      q.pop_min(act_node, act_dist);

      // iterate over all neighbours of the node
      vector<pair<int, double>> edges = graph->get_edges(act_node);

      for(int j = 0; static_cast<size_t>(j) < edges.size(); ++j) {
        double alt_dist = act_dist + edges.at(j).second;
        int v = edges.at(j).first;
        if (alt_dist < dist[v]) {
          dist[v] = alt_dist;
          q.change_priority(v, alt_dist);
        }
      }
    }

    return(dist);
  }

  /* Returns the average length of shortest paths originating from the given
     node to each of the other nodes, omitting the ones that cannot be reached.

     Remark:
     -1 is returned if the node is disconnected from other nodes.
  */
  double get_mean_shortest_path_length(int from_node = 0) {
    vector<double> dists = get_shortest_dists(from_node);

    double sum = 0;
    double count = 0;
    for(int i = 1; i < graph->get_node_count(); ++i) {
      if (dists[i] < INF) {
        sum += dists[i];
        count += 1;
      }
    }

    if (count > 0) {
      return(sum / count);
    }
    else {
      return(-1);
    }
  }
};


/* Allows the matrix instances to be printed to output streams in a human
   readable way. */
ostream& operator<<(ostream& out, Matrix& mat) {
  for(int i = 0; i < mat.get_rows(); ++i) {
    for(int j = 0; j < mat.get_cols(); ++j) {
      out << mat.get(i, j);
      if (j < (mat.get_cols() - 1)) {
        out << ", ";
      }
      else {
        out << endl;
      }
    }
  }

  return(out);
}


/* Allows graphs to be printed to output streams in a human readable way. */
ostream& operator <<(ostream &out, Graph &grr) {
  out << grr.get_adjacency_matrix();

  return(out);
}


/* Main function.
   Calculates typical values for 0.2 and 0.4 density random 50x50 matrices with
   1-10 long edges. */
int main() {
  const double densities[] = {0.2, 0.4};

  for(int graph_index = 0; graph_index < 2; ++graph_index) {
    cout << "Graphs with density " << densities[graph_index] << endl;
    cout << "---------------------------------------" << endl;

    double sum = 0.0;
    int count = 0;

    /* By changing the seed different, but reproducible random graphs are
       generated, over which the mean shortest path is calculated. */
    for(int seed = 0; seed < 20; ++seed) {
      Graph graph(50, densities[graph_index], 1.0, 10.0, seed);

      ShortestPathFinder path_finder = ShortestPathFinder(graph);
      double act_mean_length = path_finder.get_mean_shortest_path_length();
      cout << "average dist:" << act_mean_length << endl;
      sum += act_mean_length;
      ++count;
    }
    cout << "---------------------------------------" << endl;
    cout << "Summary: mean was " << sum / count << endl << endl;
  }

  return 0;
}
