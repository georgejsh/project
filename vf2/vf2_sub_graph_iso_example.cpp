//=======================================================================
// Copyright (C) 2012 Flavio De Lorenzi (fdlorenzi@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//=======================================================================

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/vf2_sub_graph_iso.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/mcgregor_common_subgraphs.hpp> 
using namespace boost;
int cou=0;
template <typename Graph1,
          typename Graph2>
  struct print_callback {

    print_callback(const Graph1& graph1, const Graph2& graph2) 
      : graph1_(graph1), graph2_(graph2) {}

    template <typename CorrespondenceMap1To2,
              typename CorrespondenceMap2To1>
    bool operator()(CorrespondenceMap1To2 f, CorrespondenceMap2To1) const {

      // Print (sub)graph isomorphism map
      //BGL_FORALL_VERTICES_T(v, graph1_, Graph1) 
       // std::cout << '(' << get(map1_, v) << ", " 
         //         << get(map2_, get(f, v)) << ") ";

      //std::cout << std::endl;
      cou++;
      return true;
    }

  private:
    const Graph1& graph1_;
    const Graph2& graph2_;
    //const VertexIndexMap1& map1_;
    //const VertexIndexMap2& map2_;
  };

int main() {

  typedef adjacency_list<setS, vecS, bidirectionalS> graph_type;
  
  // Build graph1
  int num_vertices1 = 8; 
  std::cin>>num_vertices1;
  graph_type graph1(num_vertices1);
  int num_edges;
  std::cin>>num_edges;
  for(int i=0;i<num_edges;i++){
    int a,b;
    std::cin>>a>>b;
    add_edge(a,b,graph1);
  }
  /*add_edge(0, 6, graph1); add_edge(0, 7, graph1);
  add_edge(1, 5, graph1); add_edge(1, 7, graph1);
  add_edge(2, 4, graph1); add_edge(2, 5, graph1); add_edge(2, 6, graph1);
  add_edge(3, 4, graph1);
  */
  // Build graph2
  int num_vertices2 = 9; 
  std::cin>>num_vertices2;
  graph_type graph2(num_vertices2);
  std::cin>>num_edges;
  for(int i=0;i<num_edges;i++){
    int a,b;
    std::cin>>a>>b;
    add_edge(a,b,graph2);
  }/*add_edge(0, 6, graph2); add_edge(0, 8, graph2);
  add_edge(1, 5, graph2); add_edge(1, 7, graph2);
  add_edge(2, 4, graph2); add_edge(2, 7, graph2); add_edge(2, 8, graph2);
  add_edge(3, 4, graph2); add_edge(3, 5, graph2); add_edge(3, 6, graph2);
  */
  // Create callback to print mappings
 //vf2_print_callback<graph_type, graph_type> callback(graph1, graph2);
 print_callback<graph_type, graph_type> callback(graph1, graph2);
  // Print out all subgraph isomorphism mappings between graph1 and graph2.
  // Vertices and edges are assumed to be always equivalent.
  vf2_subgraph_iso(graph1, graph2, callback);
  std::cout<<cou<<"\n";
  return 0;
}