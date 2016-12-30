

#ifndef MLG_H
#define MLG_H

#include <boost/graph/adjacency_list.hpp>

#include <string>

using namespace std;
using namespace boost;

/*---------------- Graph --------------------------------*/
struct VertexProperties  {
  VertexProperties() : name("0"){}
  VertexProperties(string const & name) : name(name){}
  string name;  // node label, to be distinguished from the node id
};

//--------------------------------------------------------------------
typedef adjacency_list<vecS, vecS,/*directedS,undirectedS bidirectionalS*/bidirectionalS,VertexProperties
		       /*, EdgeProperties*/> GraphBase;   // The first vecS is for the edge list. It accoutns form multiedges (parallel edges btw two distinct nodes) and self-loops.

typedef adjacency_list<vecS, vecS,/*directedS,undirectedS bidirectionalS*/undirectedS,VertexProperties
		       /*, EdgeProperties*/> GraphBase_undirected;   // The first vecS is for the edge list. It accoutns form multiedges (parallel edges btw two distinct nodes) and self-loops.


typedef graph_traits<GraphBase>::vertex_iterator vertex_iterator;
typedef graph_traits<GraphBase>::out_edge_iterator edge_iterator;
typedef graph_traits<GraphBase>::in_edge_iterator in_edge_iterator;
typedef graph_traits<GraphBase>::edge_iterator graph_edge_iterator;
typedef graph_traits<GraphBase>::edge_descriptor Edge;
typedef graph_traits<GraphBase>::vertex_descriptor Vertex;

typedef graph_traits<GraphBase_undirected>::vertex_iterator vertex_iterator_und;
typedef graph_traits<GraphBase_undirected>::out_edge_iterator g_edge_iterator_und;
typedef graph_traits<GraphBase_undirected>::edge_iterator edge_iterator_und;
typedef graph_traits<GraphBase_undirected>::edge_descriptor Edge_und;
typedef graph_traits<GraphBase_undirected>::vertex_descriptor Vertex_und;

struct Graph : public GraphBase {};
struct Graph_undirected : public GraphBase_undirected {};


#endif