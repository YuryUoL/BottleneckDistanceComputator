#ifndef HEADER_HPP_INCLUDED
#define HEADER_HPP_INCLUDED


#include <vector>

// just for output
#include <iostream>
#include <boost/foreach.hpp>

#include <boost/geometry/geometry.hpp>
#include <boost/geometry/geometries/adapted/std_array.hpp>
#include <boost/graph/undirected_graph.hpp>
#include <boost/graph/exterior_property.hpp>
#include <boost/graph/floyd_warshall_shortest.hpp>
//#include <fmt/ranges.h>

namespace constants
{
const int Dim = 4;
const int null_point_index = -1;
// ... other related constants

}

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

BOOST_GEOMETRY_REGISTER_STD_ARRAY_CS(bg::cs::cartesian)

typedef std::array<double, constants::Dim> point;
typedef bg::model::box<point> box;
typedef std::pair<point, unsigned> value;

typedef double t_weight;

// define the graph type
typedef boost::property<boost::edge_weight_t, t_weight> EdgeWeightProperty;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
                              boost::no_property, EdgeWeightProperty> Graph;
typedef boost::graph_traits<Graph>::adjacency_iterator AdjacencyIterator;

typedef boost::property_map<Graph, boost::edge_weight_t>::type WeightMap;

// Declare a matrix type and its corresponding property map that
// will contain the distances between each pair of vertices.
typedef boost::exterior_vertex_property<Graph, t_weight> DistanceProperty;
typedef DistanceProperty::matrix_type DistanceMatrix;
typedef DistanceProperty::matrix_map_type DistanceMatrixMap;





#endif // HEADER_HPP_INCLUDED
