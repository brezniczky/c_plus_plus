/*
  Implementation of Prim's (and Dijsktra's) algorithm.

  The demonstrated method is Prim's, over a graph read in from a file.

  Main classes:

  Matrix
  Graph
  (SimplePriorityQueue)
  (ShortestPathFinder)
  MinSpanningTreeFinder


  Most classes have the "<<" operator defined mainly for verification purposes.

  Future improvements:
  - define (x, y) operator on the matrix class

  - the matrix could be triangular, since it is symmetric - memory consumption
    could be reduced

  - further handling of exceptions/unsuccessful execution

  - better encapsulate the boolean array and the queue or replace the queue
    class with a proper minheap solution
*/

#include <iostream>
#include <vector>
#include <random>
#include <map>
#include <fstream>


using namespace std;

/* Name of the file to read in the edges of the searched graph from. */
const string kInputFilename =
  "/home/janca/c++/C++ for Programmers/Homework3/sample_test_data.txt";

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
  Matrix(int rows, int cols) {
    rows_ = rows;
    cols_ = cols;
    values_ = new double[rows * cols];
  }

  /* Creates and a matrix of the given size.
     The values are allocated and initialized as specified. */
  Matrix(int rows, int cols, const double initial_weight) : Matrix(rows, cols) {
    int ncells = rows * cols;

    for(int i = 0; i < ncells; ++i) {
      values_[i] = initial_weight;
    }
  }

  /* Destructor: frees up dynamically allocated memory. */
  ~Matrix() {
    delete[](values_);
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


/* An edge selection class for symmetric edges.

   Stores edges in (node1, node2) form, where node1 < node2.

   The declaration also allows to overload the << operator for this specific
   class.
*/
class EdgeSelection {
 private:
   multimap<int, int> items_;
 public:
  /* Adds a (node1, node2) to the map, while ensuring node1 < node2 to enhance
     readability. */
  void emplace(int node1, int node2) {
    if (node1 > node2) {
      items_.emplace(node2, node1);
    } else {
      items_.emplace(node1, node2);
    }
  }

  /* Allows access to the contained items but prevents inconsistent
     modifications. */
  const multimap<int, int>& GetItems() {
    return(items_);
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

  /* Sets an edge in the graph enforcing that the graph stays undirected. */
  inline void Set(int start_node, int target_node, double weight) {
    adjacency_->Set(start_node, target_node, weight);
    adjacency_->Set(target_node, start_node, weight);
  }

  /* Initializes a symmetric adjacency matrix with random weights. */
  void RandomInit(double density, double min_dist, double max_dist, int seed) {

    default_random_engine generator(seed);
    uniform_real_distribution<double> edge_distribution(min_dist, max_dist);
    uniform_real_distribution<double> is_connected_dist(0.0, 1.0);

    for(int i = 0; i < adjacency_->GetRows(); ++i) {
      Set(i, i, 0);
      for(int j = i + 1; j < adjacency_->GetCols(); ++j) {
        // the adjacency matrix of an undirected graph is symmetric, set the
        // values equally for both directions
        if (is_connected_dist(generator) < density) {
          double dist = edge_distribution(generator);
          Set(i, j, dist);
        } else {
          Set(i, j, kInf);
        }
      }
    }
  }

 public:
  /* Creates the graph with random density and distances as described by the
     given parameters. */
  Graph(int nodes = 50, double density = 0.1,
        double min_dist = 1.0, double max_dist = 10.0, int seed = 0) {

    adjacency_ = new Matrix(nodes, nodes, 0.0);
    RandomInit(density, min_dist, max_dist, seed);
  }

  /* Creates the graph according to the contents found in the file.

     A simple file structure is mandated:

     - first line: number of nodes
     - any further line: defines an edge in the form "node1 node2 weight"

     The edges not defined are considered disonnected (i.e. get an infinite
     weight assigned). */
  Graph(string filename) {
    ifstream filestream(filename);
    int nodes;
    filestream >> nodes;
    adjacency_ = new Matrix(nodes, nodes, kInf);

    while(!filestream.eof()) {
      int start_node, target_node;
      double weight;
      filestream >> start_node >> target_node >> weight;
      adjacency_->Set(start_node, target_node, weight);
    }
  };

  /* Destructor: frees up dynamically allocated memory. */
  ~Graph() {
    delete adjacency_;
  }

  const int kInvalidNode = -1;

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

  /* Returns all edges in the graph. */
  void GetEdges(vector<pair<int, int>>& edges) {
    for(int i = 0; i < adjacency_->GetRows(); i++) {
      for(int j = 0; j < adjacency_->GetCols(); j++) {
        if (adjacency_->Get(i, j) < kInf) {
          edges.push_back(pair<int, int>(i, j));
        }
      }
    }
  }

  /* Returns all the edges originating from a given node and their weights. */
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
  multimap<Priority, Key> item_by_priority_;
 public:
  /* Adds an item with the given priority. */
  void Push(const Key item, const Priority priority) {
    item_by_priority_.emplace(priority, item);
  }

  /* Removes the item identified by the given key. */
  bool Erase(const Key item) {
    typename multimap<Priority, Key>::iterator it;
    for(it = item_by_priority_.begin(); it != item_by_priority_.end(); ++it) {
      if (it->second == item) {
        item_by_priority_.erase(it);
        return(true);
      }
    }
    return(false);
  }

  /* Faster overload of removal, knowing the priority value allows to leverage
     the priority-based ordering. Useful if priorities tend to be diverse. */
  bool Erase(const Key item, const Priority priority) {
    for(auto it = item_by_priority_.find(priority);
      (it != item_by_priority_.end()) && (it->first==priority);
      ++it) {

      if (it->second == item) {
        item_by_priority_.erase(it);
        return(true);
      }
    }
    return(false);
  }

  /* Removes and retrieves the highest priority item
     (i.e. lowest priority value). */
  bool PopMin(Key& key, Priority& priority) {
    auto it = item_by_priority_.begin();
    if (it != item_by_priority_.end()) {
      priority = it->first;
      key = it->second;
      item_by_priority_.erase(it);
      return(true);
    } else {
      return(false);
    }
  }

  /* Changes the priority associated with the item as specified. */
  void ChangePriority(const Key& key, const Priority& new_priority) {
    Erase(key);
    Push(key, new_priority);
  }

  /* Faster overload: changes the priority associated with the item as
     specified but requires knowledge of the old priority. */
  void ChangePriority(const Key& key,
    const Priority& old_priority, const Priority& new_priority) {

    Erase(key, old_priority);
    Push(key, new_priority);
  }

  /* Tells if the queue is empty. */
  bool Empty() {
    return(item_by_priority_.empty());
  }

  /* Allows priority queues to be printed to output streams in a human readable
     way. */
  ostream& operator<<(ostream& out) {
    for(auto it : item_by_priority_) {
      cout << it.first << ":" << it.second << endl;
    }
    return(out);
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
      if (i != from_node) {
        dist.push_back(kInf);
      } else {
        dist.push_back(0);
      }
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
          q.ChangePriority(v, dist[v], alt_dist);
          dist[v] = alt_dist;

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
    } else {
      return(-1);
    }
  }
};

/* Class to find minimum spanning trees (MSTs) in a given graph. */
class MinSpanningTreeFinder {
 private:
  Graph* graph_;

  double* min_costs;   // ~ C
  int* min_prev_nodes; // ~ E, only store the previous node index
  SimplePriorityQueue<int, double> free_vertices; // ~ Q
  bool* is_free_vertex; // speeds up checking if a node is enqueued

  /* Initializes member variables for the search. */
  void InitializeSearch(double& total_cost) {
    int nodes = graph_->GetNodeCount();
    for(int i = 0; i < nodes; ++i) {
      min_costs[i] = kInf;
      min_prev_nodes[i] = graph_->kInvalidNode; // just in case - never checked
      free_vertices.Push(i, kInf);
      is_free_vertex[i] = true;
    }
    total_cost = 0;
  }

  /* Performs necessary updates to members when a new node joins the tree. */
  void UpdateOnNewNode(int new_node) {
    is_free_vertex[new_node] = false;

    // from the new node, further nodes may become available
    vector<pair<int, double>> new_edges = graph_->GetEdgesFrom(new_node);
    for(size_t i = 0; i < new_edges.size(); ++i) {
      int target = new_edges[i].first;
      int target_cost = new_edges[i].second;

      if (is_free_vertex[target] && (target_cost < min_costs[target])) {
        free_vertices.ChangePriority(target, min_costs[target], target_cost);
        min_prev_nodes[target] = new_node;
        min_costs[target] = target_cost;
      }
    }
  }

 public:
  /* Creates the instance, preparing for operations on the given graph. */
  MinSpanningTreeFinder(Graph& g) {
    graph_ = &g;

    // RAII
    int nodes = graph_->GetNodeCount();

    min_costs = new double[nodes];
    min_prev_nodes = new int[nodes];
    is_free_vertex = new bool[nodes];
  }

  /* Destructor: extends default behaviour to free up dynamically allocated
     memory. */
  ~MinSpanningTreeFinder() {
    delete[](min_costs);
    delete[](min_prev_nodes);
    delete[](is_free_vertex);
  }

  /* Searches for an MST in the graph. Implements Prim's algorithm
     (see https://en.wikipedia.org/wiki/Prim%27s_algorithm). */
  bool PerformSearch(EdgeSelection& edges, double& total_cost) {
    // edges ~ F

    // initialize
    InitializeSearch(total_cost);

    int new_node;
    double min_length;
    bool is_first = true;

    // look for the closest nodes to be added: try to find one
    while(free_vertices.PopMin(new_node, min_length)) {

      if (min_costs[new_node] < kInf) {
        // take note of the edge pointing to easiest to reach node
        int new_prev_node = min_prev_nodes[new_node];
        edges.emplace(new_prev_node, new_node);
        total_cost += min_costs[new_node];
      } else if (is_first) {
        // for the first node it is okay - every node is disconnected from an
        // empty tree
        is_first = false;
      } else {
        // an infinitely long hop (unless for the first node) means a
        // disconnected graph
        return(false);
      }

      UpdateOnNewNode(new_node);
    }

    return(true);
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
      } else {
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

/* Allows edge selections to be printed to output streams in a human readable
   way. */
ostream& operator <<(ostream& out, EdgeSelection& edges) {
  const multimap<int, int>* items = &edges.GetItems();
  multimap<int, int>::const_iterator iter = items->begin();
  while (iter != items->end()) {
    out << iter->first << "\t" << iter-> second << endl;
    ++iter;
  }
  return(out);
}


/* Main function.
   Finds the minimum spanning tree for a graph read in from the input file. */
int main() {

  // please modify kInputFilename according to your location of the file as
  // needed
  Graph gr(kInputFilename);

  MinSpanningTreeFinder tree_finder(gr);
  double total_cost;
  EdgeSelection min_spanning_tree;

  if (tree_finder.PerformSearch(min_spanning_tree, total_cost)) {
    cout << "Prim's algorithm finished successfully" << endl;
    cout << "Total cost: " << total_cost << endl;
    cout << "Edges:" << endl;
    cout << min_spanning_tree;
  } else {
    cout << "Prim's algorithm failed, the graph is disconnected." << endl;
  }

  return 0;
}
