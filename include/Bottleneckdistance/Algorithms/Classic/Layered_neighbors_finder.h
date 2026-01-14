#ifndef LAYERED_NEIGHBORS_FINDER_H
#define LAYERED_NEIGHBORS_FINDER_H

#include <unordered_set>
#include <vector>
#include <algorithm>  // for std::max
#include <cmath>  // for std::abs
#include <cassert>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include "header.hpp"


//#include "KDTreeController.h"
//#include "Dataset_Graph.h"

template<typename graphType,typename treeController>
class Layered_neighbors_finder
{
public:

    Layered_neighbors_finder(graphType & gg, double rr):g(gg),r(rr){}

    inline void add(int v_point_index, int vlayer)
    {
        for (int l = neighbors_finder.size(); l <= vlayer; l++)
        {
            neighbors_finder.push_back(treeController(this->g, r));
        }
        neighbors_finder.at(vlayer).add(v_point_index);
    }

    inline int pull_near(int u_point_index, int vlayer)
    {
        if (static_cast<int> (neighbors_finder.size()) <= vlayer)
            return constants::null_point_index;
        return neighbors_finder.at(vlayer).pull_near(u_point_index);
    }

    inline int vlayers_number() const
    {
        return static_cast<int> (neighbors_finder.size());
    }
//    std::vector<int> pull_all_near(int u_point_index)
//    {
//        std::vector<int> all_pull;
//        int last_pull = pull_near(u_point_index);
//        while (last_pull != null_point_index())
//        {
//            all_pull.emplace_back(last_pull);
//            last_pull = pull_near(u_point_index);
//        }
//        return all_pull;
//    }

    //! Delete all the shit
//    virtual ~Layered_neighbors_finder() {
//        for (auto it = neighbors_finder.begin(); it != neighbors_finder.end(); it++)
//        {
//            delete *it;
//        }
//
//    }

protected:

private:
    double r;
    graphType & g;
    std::vector<treeController> neighbors_finder;

};

#endif // LAYERED_NEIGHBORS_FINDER_H
