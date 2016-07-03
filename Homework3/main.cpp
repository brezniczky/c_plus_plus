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

  - examine multiset as a PriorityQueue implementation option (can't remember
    where it was suggested)
    http://en.cppreference.com/w/cpp/container/multiset
    "a sorted set of objects of type Key"
    "Search, insertion, and removal operations have logarithmic complexity."
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

  Matrix(int rows, int cols, const double initial_weight) : Matrix(rows, cols) {
    int ncells = rows * cols;

    for(int i = 0; i < ncells; ++i) {
      values_[i] = initial_weight;
    }
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


/* An edge selection class.

   Stores edges in (node1, node2) form, where node1 < node2.

   The declaration also allows to overload the << operator for this specific
   class.
*/
class EdgeSelection {
 private:
   multimap<int, int>* items_;
 public:
  EdgeSelection() {
    items_ = new multimap<int, int>();
  }

  ~EdgeSelection() {
    delete(items_);
  }

  /* Adds a (key, value) to the map, while ensuring key < value to enhance
     readability. */
  void emplace(int key, int value) {
    if (key > value) {
      items_->emplace(value, key);
    } else {
      items_->emplace(key, value);
    }
  }

  /* Allows access to the contained items but prevents inconsistent 
     modifications */
  const multimap<int, int>& GetItems() {
    return(*items_);
  }

  /* TODO: to be removed */
  const void Exclude(EdgeSelection& edges) {
    cout << "count before:" << items_->size() << endl;
    auto items = edges.GetItems();
    for(auto edge : items) {
      items_->erase(edge.first, edge.second);

      items_->erase()
    }
    cout << "count after:" << items_->size() << endl;
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

    adjacency_ = new Matrix(nodes, nodes);;
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

  /* Returns all edges in the graph. */
  void GetEdges(EdgeSelection& edges) {
    for(int i = 0; i < adjacency_->GetRows(); i++) {
      for(int j = 0; j < adjacency_->GetCols(); j++) {
        if (adjacency_->Get(i, j) < kInf) {
          edges.emplace(i, j);
        }
      }
    }
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

  /* TODO: to be removed */
  /* Finds the minimum cost edge originating from a given node, considering only
     the allowed nodes. The i'th node is allowed if allowed_nodes[i] is true,
     allowed_nodes must contain an item per each node.

     Returns false iff there was no connected and allowed node. */
  bool GetMinCostEdgeFrom(
    int from_node, const vector<bool> allowed_nodes,
    int& min_target, double& min_cost) {

    min_cost = kInf;
    min_target = -1;

    for(int i = 0; static_cast<size_t>(i) < allowed_nodes.size(); ++i) {

      if (allowed_nodes[i] &&
          (adjacency_->Get(from_node, i) < min_cost)) {
        min_cost = adjacency_->Get(from_node, i);
        min_target = i;
      }
    }

    return(min_cost < kInf);
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


class MinSpanningTreeFinder {
 private:
  Graph* graph_;

  /* Removes unused edges from the graph */
  void RemoveUnusedEdges(EdgeSelection& edges,
                         const vector<bool> allowed_nodes) {

    EdgeSelection unused;

    for(auto iter in edges.GetItems()) {
      if (!allowed_nodes[iter->first] && !allowed_nodes[iter->second]) {
        unused.emplace(iter->first, iter->second);
      }
    }

    edges.Exclude(unused);
  }

  /* Finds the edge with a minimum weight pointing out from the known vertices
     (i.e. to those marked with true in allowed_nodes[]. */
  bool GetMinCostEdgeFrom(
    EdgeSelection& edges, const vector<bool> allowed_nodes,
    int& min_target, double& min_cost) {

    min_cost = kInf;
    min_target = -1;

    EdgeSelection unused;

    multimap<int, int> items;
    items = edges.GetItems();

    for(auto item : items) {
      // the edge should connect an already included node and one that isn't,
      // should also be so far minimal
      if (allowed_nodes[item.first] == allowed_nodes[item.second]) {
        if (!allowed_nodes[item.first]) {
          /* connects two nodes which have both already been included - mark it
             for removal */
          unused.emplace(item.first, item.second);
        }
      } else if (graph_->GetEdgeLength(item.first, item.second) < min_cost) {
        min_cost = graph_->GetEdgeLength(item.first, item.second);
        min_target = allowed_nodes[item.first] ? item.first : item.second;
      }
    }

    edges.Exclude(unused);

    return(min_cost < kInf);
  }

 public:
  /* Creates the instance, preparing for operations on the given graph. */
  MinSpanningTreeFinder(Graph& g) {
    graph_ = &g;
  }

  /* Searches the graph for a minimum spanning tree according to Prim's
     algorithm (https://en.wikipedia.org/wiki/Prim%27s_algorithm).

     Returns true iff the MST was found. */
  bool RunPrim(EdgeSelection& edges, double& total_cost) {
    int nodes = graph_->GetNodeCount();

    // boolean matrix: which nodes are yet to reach
    vector<bool> node_allowed(nodes);
    for(int i = 1; i < nodes; ++i) {
      node_allowed[i] = true;
    }

    // start with an arbitrary node: e.g. the first
    node_allowed[0] = false;
    int act_node = 0;

    EdgeSelection all_edges;
    graph_->GetEdges(all_edges);

    /* build a spanning tree:
       in each iteration add the closest reachable vertex until all of them are
       included in the spanning tree */
    for(int nodes_left = nodes - 1; nodes_left > 0; --nodes_left) {
      int min_target;
      double min_cost;

      if (GetMinCostEdgeFrom(all_edges, node_allowed,
                             min_target, min_cost)) {

        edges.emplace(act_node, min_target);
        total_cost += min_cost;

        node_allowed[min_target] = false;
        act_node = min_target;
      } else {
        // the graph is disconnected
        return(false);
      };
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

/* Allows edge selections to be printed to output streams in a human readable
   way. */
ostream& operator <<(ostream& out, EdgeSelection& edges) {
  const multimap<int, int>* items = &edges.GetItems();
  multimap<int, int>::const_iterator iter = items->begin();
  while (iter != items->end()) {
    cout << iter->first << "\t" << iter-> second << endl;
    ++iter;
  }
  return(out);
}


/* Main function.
   Finds the minimum spanning tree for a graph read in from the input file. */
int main() {
  Graph gr(kInputFilename);

  MinSpanningTreeFinder tree_finder(gr);
  double total_cost;
  EdgeSelection min_spanning_tree;

  if (tree_finder.RunPrim(min_spanning_tree, total_cost)) {
    cout << "Prim's algorithm finished successfully" << endl;
    cout << "Total cost:" << total_cost << endl;
    cout << "Edges:" << endl;
    cout << min_spanning_tree;
  }
  else {
    cout << "Prim's algorithm failed, the graph is disconnected." << endl;
  }

  return 0;
}
