#include <iostream>
#include <vector>
#include <random>
// #include <queue>
#include <map>

/*
ideas:
- define >> operator again, on the matrix class
  useful for manual verification

- define (x, y) operator on the matrix class

- the matrix could be triangular, since it is symmetric

- use size_t where int is an indexer

- new TODO:
  priority queue: always proceed with the closest node from "u"

- priority queue:
  Compare argument

mention:
- node indexes are zero-based
*/

/*

potential partial interface definition for a Graph could be:

Class Graph:

V (G): returns the number of vertices in the graph
E (G): returns the number of edges in the graph
adjacent (G, x, y): tests whether there is an edge from node x to node y.
neighbors (G, x): lists all nodes y such that there is an edge from x to y.
add (G, x, y): adds to G the edge from x to y, if it is not there.
delete (G, x, y): removes the edge from x to y, if it is there.
get_node_value (G, x): returns the value associated with the node x.
set_node_value( G, x, a): sets the value associated with the node x to a.
get_edge_value( G, x, y): returns the value associated to the edge (x,y).
set_edge_value (G, x, y, v): sets the value associated to the edge (x,y) to v.

Class PriorityQueue

chgPrioirity(PQ, priority): changes the priority (node value) of queue element.
minPrioirty(PQ): removes the top element of the queue.
contains(PQ, queue_element): does the queue contain queue_element.
Insert(PQ, queue_element): insert queue_element into queue
top(PQ):returns the top element of the queue.
size(PQ): return the number of queue_elements.

Class ShortestPath

vertices(List): list of vertices in G(V,E).
path(u, w): find shortest path between u-w and returns the sequence of vertices representing shorest path u-v1-v2-…-vn-w.
path_size(u, w): return the path cost associated with the shortest path.

*/

using namespace std;

const double INF = numeric_limits<double>::infinity();

typedef int node_id_t;

class Matrix {
private:
  double* values;
  int rows;
  int cols;

  inline int get_index(int row, int col) {
    return(row * cols + col);
  }

public:
  Matrix(int rows, int cols) {
    this->rows = rows;
    this->cols = cols;
    values = new double[rows * cols];
    // TODO: initialize to false?
  }

  ~Matrix() {
    delete values;
  }

  double inline get(int row, int col) {
    return(values[get_index(row, col)]);
  }

  // TODO: define the = operator instead of set()
  void inline set(int row, int col, const double value) {
    values[get_index(row, col)] = value;
  }

  inline int get_rows() {
    return(rows);
  }

  inline int get_cols() {
    return(cols);
  }
};

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


class Graph {
private:
  // adjacency (aka. connectivity) matrix of the graph
  //
  Matrix* adjacency;

  void random_init(double density, double mindist, double maxdist, int seed) {

    default_random_engine generator(seed);
    uniform_real_distribution<double> dist_distribution(mindist, maxdist);
    uniform_real_distribution<double> is_connected_dist(0.0, 1.0);

    for(int i = 0; i < adjacency->get_rows(); ++i) {
      for(int j = i + 1; j < adjacency->get_cols(); ++j) {
        if (is_connected_dist(generator) < density) {
          double dist = dist_distribution(generator);

          // the adjacency matrix of an undirected graph is symmetric
          adjacency->set(i, j, dist);
          adjacency->set(j, i, dist);
        } else {
          adjacency->set(i, j, INF);
          adjacency->set(j, i, INF);
        }
      }
    }
  }

  void init_adjacency(int nodes) {
    adjacency = new Matrix(nodes, nodes);
  }
public:
  /*  The random graph procedure should have edge density
  as a parameter and distance range as a parameter. */
  // TODO: What is the distribution supposed to be?
  //       I'll go with uniform first.
  Graph(int nodes = 50, double density = 0.1,
        double mindist = 1.0, double maxdist = 10.0, int seed = 0) {
    // TODO: do I have perform this init here or is the parameterless autocalled?
    init_adjacency(nodes);
    random_init(density, mindist, maxdist, seed);
  }

  ~Graph() {
    delete adjacency;
  }

  int get_node_count() {
    return(adjacency->get_rows());
  }

  Matrix& get_adjacency_matrix() {
    return(*adjacency);
  }

  double get_edge_length(int node1, int node2) {
    return(adjacency->get(node1, node2));
  }

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

//class Path:vector<int> {
//public:
//  double get_length(Graph g) {
//    double sum = 0.0;
//    for(int i = (this->size() - 1); i > 1; --i) {
//      sum += g.get_edge_length((*this)[i], ((*this)[i - 1]));
//    }
//
//    return(sum);
//  };
//};


template<class Key, class Priority> class SimplePriorityQueue {
  /* first implementation: an array-based version for verifying the rest of the
     algorithm and to reach a prototype implementation

     slow and bad, but it is only to give correct results while the real
     solution is developed
  */
private:
  map<Key, Priority> priority_by_key;
public:
  Key get_max_item() {
    if (priority_by_key.empty()) {
      return(NULL);
    };

    typename map<Key, Priority>::iterator it = priority_by_key.begin();
    typename map<Key, Priority>::iterator max_item;

    for(; it != priority_by_key.end(); ++it)
      if (it->second > max_item->second)
        max_item = it;

    return(max_item->first);
  }

  void push(const Key &key, const Priority &priority) {
    priority_by_key.emplace(key, priority);
  }

  void dump() {
    typename map<Key, Priority>::iterator it = priority_by_key.begin();

    while(it != priority_by_key.end()) {
      cout << it->first << ":" << it -> second << endl;
      ++it;
    }
  }

  bool erase(const Key key) {
    return(priority_by_key.erase(key) > 0);
  }

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

//    dump();
    if (!erase(key)) {
      cout << "could not delete key " << key << endl;
//      dump();
      erase(key);
    }

    return(true);
  }

  void change_priority(const Key& key, const Priority& priority) {
    erase(key);
    push(key, priority);
  }

  bool empty() {
    return(priority_by_key.empty());
  }
};

class ShortestPathFinder {
private:
  Graph* graph;
public:
  ShortestPathFinder(Graph &g) {
    graph = &g;
  }

  vector<double> get_shortest_dists(int from_node) {
    int node_count = graph->get_node_count();
    vector<double> dist = vector<double>();
//    double prev[node_count];

    SimplePriorityQueue<int, double> q = SimplePriorityQueue<int, double>();

    for(int i = 0; i < node_count; ++i) {
      if (i != from_node)
        dist.push_back(INF);
      else
        dist.push_back(0);
//      prev[i] = -1;
      q.push(i, dist[i]);
    }

    dist[from_node] = 0;

    while (!q.empty()) {
      int act_node;
      double act_dist;
      q.pop_min(act_node, act_dist);

//      cout << "checking next node:" << act_node << endl;
//      cout << "current dist:" << act_dist << endl;

      // iterate over all neighbours of the node
      vector<pair<int, double>> edges = graph->get_edges(act_node);

      for(int j = 0; static_cast<size_t>(j) < edges.size(); ++j) {
        double alt_dist = act_dist + edges.at(j).second;
        int v = edges.at(j).first;
        if (alt_dist < dist[v]) {
          dist[v] = alt_dist;
          q.change_priority(v, alt_dist);
//            prev[v] = u;
        }
      }
    }

    return(dist);
  }

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

ostream& operator <<(ostream &out, Graph &grr) {
  out << grr.get_adjacency_matrix();

  return(out);
}

void test() {
  // 1->2: 2 units
  // 1->3: 4 units
  Graph g = Graph(4, 0);

  cout << g.get_adjacency_matrix() << "---" << endl;

  // 1->3, 3->2, 2->4 transitions cost 1, 2, 3 each
  //
  // thus,
  // path to #2: 1
  // path to #3: 1 + 2
  // path to #4: 1 + 2 + 3
  //
  // mean: 10 / 3 = 3.333

  g.get_adjacency_matrix().set(0, 2, 1.0);
  g.get_adjacency_matrix().set(2, 0, 1.0);

  g.get_adjacency_matrix().set(2, 1, 2.0);
  g.get_adjacency_matrix().set(1, 2, 2.0);

  g.get_adjacency_matrix().set(1, 3, 3.0);
  g.get_adjacency_matrix().set(3, 1, 3.0);

  cout << g.get_adjacency_matrix();

  // mean shortest path length: 3 units
  ShortestPathFinder sp = ShortestPathFinder(g);

  if (sp.get_mean_shortest_path_length() != (10.0 / 3)) {
    cout << "test 1 failed : " << sp.get_mean_shortest_path_length() << endl;
  }

  Graph g2 = Graph(4, 0);
  g2.get_adjacency_matrix().set(0, 1, 1.0);
  g2.get_adjacency_matrix().set(1, 0, 1.0);

  g2.get_adjacency_matrix().set(2, 1, 3.0);
  g2.get_adjacency_matrix().set(1, 2, 3.0);

  g2.get_adjacency_matrix().set(2, 3, 3.0);
  g2.get_adjacency_matrix().set(3, 2, 3.0);

  ShortestPathFinder sp2 = ShortestPathFinder(g2);
  if (sp2.get_mean_shortest_path_length() != 4) {
    cout << "test 2 failed : " << sp.get_mean_shortest_path_length() << endl;
  }
}

int main() {
  vector<Graph*> graphs = vector<Graph*>();
  graphs.push_back(new Graph(50, 0.2, 1.0, 10.0, 0));
  graphs.push_back(new Graph(50, 0.4, 1.0, 10.0, 0));

  cout << "graphs have been created" << endl;

  for(size_t i = 0; i < graphs.size(); ++i) {
    ShortestPathFinder path_finder = ShortestPathFinder(*(graphs[i]));
    cout << "average dist:" <<
            path_finder.get_mean_shortest_path_length() << endl;
  }

  return 0;
}
