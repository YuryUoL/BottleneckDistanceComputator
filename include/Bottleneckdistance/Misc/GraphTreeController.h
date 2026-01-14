#ifndef GRAPHTREECONTROLLER_H
#define GRAPHTREECONTROLLER_H

#include <list>
#include "Bipartie_Graph.h"
#include <unordered_set>


class GraphTreeController
{
public:
    GraphTreeController(Bipartie_Graph & g, double r):graph(g),r(r) {}
    virtual ~GraphTreeController() {}

    void add(int v_point_index)
    {
        NodesInLayer.push_back(v_point_index);
        // NodesInLayer.insert(v_point_index);
    }


//       for (auto it = NodesInLayer.begin(); it != NodesInLayer.end(); it++)
//        {
//
//
//
//            if (graph.distance(u_point_index, *it) <= this->r)
//            {
//                int ii = *it;
//             //!   NodesInLayer.erase(it);
//                NodesInLayer.erase(it);
//                return ii;
//
//
//            }
//        }
//
    int pull_near(int u_point_index)
    {

        //! Find smallest
        std::list<int>::iterator iterasable;
        double minndistt = std::numeric_limits<double>::max();

        for (auto it = NodesInLayer.begin(); it != NodesInLayer.end(); it++)
        {


        //! graph.distance(u_point_index, *it)
        //! *it,u_point_index
            if (graph.distance(u_point_index, *it) <= minndistt)
            {
                minndistt = graph.distance(u_point_index, *it);
                iterasable = it;

            }
        }
        if (minndistt > this->r)
        {
            return constants::null_point_index;

        }
        else
        {

            int ii = *iterasable;
            NodesInLayer.erase(iterasable);
            return ii;


        }



        //return constants::null_point_index;
    }


    std::vector<int> pull_all_near(int u_point_index)
    {
//       int summ = 0;
//       for (auto it = NodesInLayer.begin(); it != NodesInLayer.end(); it++)
//       {
//           if (graph.distance(*it, u_point_index) <= this->r)
//           {
//               summ = summ + 1;
//            }
//
//       }
//         std::cout << "Pull all near method : " << summ << std::endl;
        std::vector<int> all_pull;
        int last_pull = pull_near(u_point_index);
        while (last_pull != constants::null_point_index)
        {
            all_pull.emplace_back(last_pull);
            last_pull = pull_near(u_point_index);
        }
//        std::cout << "However found: " << all_pull.size() << std::endl;
//        for (int i = 0; i < all_pull.size();i++)
//        {
//            std::cout << all_pull[i] << ", ";
//
//        }
//        std::cout << std::endl;
        return all_pull;
    }


protected:


private:
    std::list<int> NodesInLayer;
//! std::unordered_set<int> NodesInLayer;
    Bipartie_Graph & graph;
    double r;
};

#endif // GRAPHTREECONTROLLER_H
