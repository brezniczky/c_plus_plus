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

  - handling of exceptions/unsuccessful execution

  - use priority_queue standard template via inheritance as mentioned at
    http://stackoverflow.com/questions/19467485/c-priority-queue-removing-element-not-at-top
    hm... looks like this approach will rebuild the heap on removal, making that
    an o(n) operation

    here is another approach, which is indeed fast
    http://stackoverflow.com/questions/3076163/stl-priority-queue-deleting-an-item
    but may pollute memory (fits well though in my case)

    yet another: use multiset (can't remember where it was suggested)
    http://en.cppreference.com/w/cpp/container/multiset
    "a sorted set of objects of type Key"
    "Search, insertion, and removal operations have logarithmic complexity."

*/

#include <iostream>
#include <vector>
#include <random>
#include <map>


using namespace std;


/* Infinite value constant. */
const double kInf = numeric_limits<double>::infinity();


/* Class representing a matrix of floating point values. */
class Matrix {
 private:
  double* values_;
  int rows_;
  int cols_;

  /* Quickly retrieves the index of the memory slot for a given cell. */
  inline int GetIndex(int row, int col) {
    return(row * cols_ + col);
  }

 public:
  /* Creates and a matrix of the given size.
  The values are allocated but undefined. */
  Matrix(int rows_, int cols_) {
    this->rows_ = rows_;
    this->cols_ = cols_;
    values_ = new double[rows_ * cols_];
  }

  ~Matrix() {
    delete(values_);
  }

  /* Gets the value of a given cell. */
  double inline Get(int row, int col) {
    return(values_[GetIndex(row, col)]);
  }

  /* Sets the value of a given cell. */
  void inline Set(int row, int col, const double value) {
    values_[GetIndex(row, col)] = value;
  }

  /* Tells how many rows the matrix has. */
  inline int GetRows() {
    return(rows_);
  }

  /* Tells how many columns the matrix has. */
  inline int GetCols() {
    return(cols_);
  }
};


/* A simple graph class allowing directed and undirected weighted graphs to be
   specified, and for random graph generation.

   Nodes are represented by integer values ranging from 0..(n - 1), where n is
   the number of nodes in the graph.
*/
class Graph {
 private:
  /* Adjacency (aka. connectivity) matrix of the graph. Weights represent
     distances. "inf" weights are to be used between the non-adjacent nodes. */
  Matrix* adjacency_;

  /* Initializes a symmetric adjacency matrix with random weights. */
  void RandomInit(double density, double min_dist, double max_dist, int seed) {

    default_random_engine generator(seed);
    uniform_real_distribution<double> edge_distribution(min_dist, max_dist);
    uniform_real_distribution<double> is_connected_dist(0.0, 1.0);

    for(int i = 0; i < adjacency_->GetRows(); ++i) {
      for(int j = i + 1; j < adjacency_->GetCols(); ++j) {
        // the adjacency matrix of an undirected graph is symmetric, set the
        // values equally for both directions
        if (is_connected_dist(generator) < density) {
          double dist = edge_distribution(generator);
          adjacency_->Set(i, j, dist);
          adjacency_->Set(j, i, dist);
        } else {
          adjacency_->Set(i, j, kInf);
          adjacency_->Set(j, i, kInf);
        }
      }
    }
  }

 public:
  /* Creates the graph with random density and distances as described by the
     given parameters. */
  Graph(int nodes = 50, double density = 0.1,
        double min_dist = 1.0, double max_dist = 10.0, int seed = 0) {

    adjacency_ = new Matrix(nodes, nodes);;
    RandomInit(density, min_dist, max_dist, seed);
  }

  /* Destructor - frees up allocated memory. */
  ~Graph() {
    delete adjacency_;
  }

  /* Tells the number of nodes in the graph. */
  int GetNodeCount() {
    return(adjacency_->GetRows());
  }

  /* Returns a reference to the adjacency matrix allowing external
     modifications. */
  Matrix& GetAdjacencyMatrix() {
    return(*adjacency_);
  }

  /* Returns the length of a given edge. */
  double GetEdgeLength(int node1, int node2) {
    return(adjacency_->Get(node1, node2));
  }

  /* Returns all the edges originating from a given node. */
  vector<pair<int, double>> GetEdgesFrom(int node) {

    vector<pair<int, double>> neighbours = vector<pair<int, double>>();

    for(int i = 0; i < adjacency_->GetCols(); i++) {
      double dist = adjacency_->Get(node, i);
      if (dist < kInf) {
        neighbours.push_back(pair<int, double>(i, dist));
      }
    }

    return(neighbours);
  }
};


/* A map-based simple priority queue implementation.

   Scales poorly for larger input, but is a useful placeholder while it needs
   no replacement. A more efficient (but also involving) version can be
   created then. Suffices for the current input set. */
template<class Key, class Priority> class SimplePriorityQueue {
 private:
  map<Key, Priority> priority_by_key_;
 public:
  /* Adds an item with the given priority. */
  void Push(const Key& key, const Priority& priority) {
    priority_by_key_.emplace(key, priority);
  }

  /* Removes the item with the given key. */
  bool Erase(const Key key) {
    return(priority_by_key_.erase(key) > 0);
  }

  /* Removes and retrieves the highest priority item
     (i.e. lowest priority value). */
  bool PopMin(Key& key, Priority& priority) {
    if (priority_by_key_.empty())
      return(false);

    typename map<Key, Priority>::iterator it = priority_by_key_.begin();

    key = it->first;
    priority = it->second;

    while(it != priority_by_key_.end()) {
      if (it->second < priority) {
        key = it->first;
        priority = it->second;
      }
      ++it;
    }

    Erase(key);

    return(true);
  }

  /* Changes the priority associated with the item as specified. */
  void ChangePriority(const Key& key, const Priority& new_priority) {
    Erase(key);
    Push(key, new_priority);
  }

  /* Tells if the queue is empty. */
  bool Empty() {
    return(priority_by_key_.empty());
  }

  /* Allows priority queues to be printed to output streams in a human readable
     way. */
  ostream& operator<<(ostream& out) {
    typename map<Key, Priority>::iterator it = priority_by_key_.begin();

    while(it != priority_by_key_.end()) {
      cout << it->first << ":" << it -> second << endl;
      ++it;
    }
  }
};


/* Class implementing Dijsktra's shortest path algorithm, with focusing only on
   the average path lengths from a given node. */
class ShortestPathFinder {
 private:
  Graph* graph_;
 public:
  /* Creates the instance, preparing for operations on the given graph. */
  ShortestPathFinder(Graph& g) {
    graph_ = &g;
  }

  /* Returns all the lengths of shortest paths between a given node and all
     the nodes, including itself.

     A simplified version of the full Dijsktra's algorithm, based on the
     pseudo-code at https://en.wikipedia.org/wiki/Dijkstra%27s_algorithm */
  vector<double> GetShortestPathLengths(int from_node) {
    int node_count = graph_->GetNodeCount();
    vector<double> dist = vector<double>();

    SimplePriorityQueue<int, double> q = SimplePriorityQueue<int, double>();

    for(int i = 0; i < node_count; ++i) {
      if (i != from_node)
        dist.push_back(kInf);
      else
        dist.push_back(0);
      q.Push(i, dist[i]);
    }

    dist[from_node] = 0;

    while (!q.Empty()) {
      int act_node;
      double act_dist;
      q.PopMin(act_node, act_dist);

      // iterate over all neighbours of the node
      vector<pair<int, double>> edges = graph_->GetEdgesFrom(act_node);

      for(int j = 0; static_cast<size_t>(j) < edges.size(); ++j) {
        double alt_dist = act_dist + edges.at(j).second;
        int v = edges.at(j).first;
        if (alt_dist < dist[v]) {
          dist[v] = alt_dist;
          q.ChangePriority(v, alt_dist);

          // keep going: this allows to check the paths to all other nodes
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
  double GetMeanShortestPathLength(int from_node = 0) {
    vector<double> dists = GetShortestPathLengths(from_node);

    double sum = 0;
    double count = 0;
    for(int i = 1; i < graph_->GetNodeCount(); ++i) {
      if (dists[i] < kInf) {
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
  for(int i = 0; i < mat.GetRows(); ++i) {
    for(int j = 0; j < mat.GetCols(); ++j) {
      out << mat.Get(i, j);
      if (j < (mat.GetCols() - 1)) {
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
ostream& operator <<(ostream& out, Graph& g) {
  out << g.GetAdjacencyMatrix();

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
      double act_mean_length = path_finder.GetMeanShortestPathLength();
      cout << "average dist: " << act_mean_length << endl;
      sum += act_mean_length;
      ++count;
    }
    cout << "---------------------------------------" << endl;
    cout << "Summary: mean was " << sum / count << endl << endl;
  }

  return 0;
}
