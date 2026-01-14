#ifndef GRAPH_MATCHING_H
#define GRAPH_MATCHING_H

#include "Layered_neighbors_finder.h"

template<typename graphType,typename treeController>
class Graph_Matching
{
public:
    Graph_Matching(graphType& g):gp(&g),r(0.),v_to_u(g.size(),constants::null_point_index),unmatched_in_u(g.size())
    {
        for (int u_point_index = 0; u_point_index < g.size(); ++u_point_index)
        {
            unmatched_in_u.insert(u_point_index);
        }
    }

    bool perfect()
    {
        return unmatched_in_u.empty();
    }
    void set_r(double r)
    {
        this->r = r;
    }

    bool augment(Layered_neighbors_finder<graphType,treeController> & layered_nf, int u_start_index, int max_depth)
    {
        // V vertices have at most one successor, thus when we backtrack from U we can directly pop_back 2 vertices.
        std::vector<int> path;
        path.emplace_back(u_start_index);
        //   std::cout << "Starting path with " << u_start_index << std::endl;
        //  std::cout << "max depth is now " << max_depth << std::endl;
        do
        {
            if (static_cast<int> (path.size()) > max_depth)
            {
                path.pop_back();
                path.pop_back();
            }
            if (path.empty())
                return false;
            int targetlayer = static_cast<int> (path.size()) / 2;
            int backofpath = path.back();
            int result = layered_nf.pull_near(backofpath, targetlayer);
            //   std::cout << "node at :" << backofpath << "of layer " << targetlayer << " resulted into " << result << std::endl;
            path.emplace_back(result);
            while (path.back() == constants::null_point_index)
            {
                path.pop_back();
                path.pop_back();
                if (path.empty())
                    return false;
                path.pop_back();
                path.emplace_back(layered_nf.pull_near(path.back(), path.size() / 2));
            }
            path.emplace_back(v_to_u.at(path.back()));
        }
        while (path.back() != constants::null_point_index);
        // if v_to_u.at(path.back()) has no successor, path.back() is an exposed vertex
        path.pop_back();
        //   std::cout << "Size of path is now " << path.size() << std::endl;
        update(path);
        return true;
    }

    void update(std::vector<int>& path)
    {
        // Must return 1.
        this->unmatched_in_u.erase(path.front());
        for (auto it = path.cbegin(); it != path.cend(); ++it)
        {
            // Be careful, the iterator is incremented twice each time
            int tmp = *it;
            v_to_u[*(++it)] = tmp;
        }
    }

    bool multi_augment()
    {
        //std::cout << "unmatched currently: " << this->unmatched_in_u.size() << std::endl;
        if (perfect())
        {
         //   std::cout << "Shutting down because perfect" << std::endl;
            return false;
        }
        Layered_neighbors_finder<graphType,treeController> layered_nf(layering());
        int max_depth = layered_nf.vlayers_number()*2 - 1;
        double rn = sqrt(gp->size());
        //std::cout << "max depth is : " << max_depth << std::endl;
        //std::cout << "unmatched at this point: " << this->unmatched_in_u.size() << std::endl;
        //std::cout << "rn is : " << rn << std::endl;
        // verification of a necessary criterion in order to shortcut if possible
        if (max_depth < 0 || (unmatched_in_u.size() > rn && max_depth >= rn))
        {
          //  std::cout << "Shutting down because verification" << std::endl;
            return false;
        }
        bool successful = false;
        std::vector<int> tries(unmatched_in_u.cbegin(), unmatched_in_u.cend());
        int itar = 0;
        for (auto it = tries.cbegin(); it != tries.cend(); it++)
            // 'augment' has side-effects which have to be always executed, don't change order
        {
          // std::cout << "Starting augementation with index " << itar << " And point: " << *it << std::endl;
           bool kkk = augment(layered_nf, *it, max_depth);
          // std::cout << "Local success is: " << kkk << std::endl;
            successful = kkk || successful;
        //std::cout << "Success is: " << successful << std::endl;

             itar = itar + 1;
        }

      //  std::cout << "Final success: " << successful << std::endl;
        return successful;
    }

    Layered_neighbors_finder<graphType,treeController> layering() const
    {
        std::vector<int> u_vertices(unmatched_in_u.cbegin(), unmatched_in_u.cend());
      //  std::cout << "U vertices size: " << u_vertices.size() << std::endl;
      //  std::cout << "Radius: " << r << std::endl;
        std::vector<int> v_vertices;
        treeController nf(*gp, r);
        for (int v_point_index = 0; v_point_index < gp->size(); ++v_point_index)
            nf.add(v_point_index);
        Layered_neighbors_finder<graphType,treeController> layered_nf(*gp, r);
        for (int layer = 0; !u_vertices.empty(); layer++)
        {
        //        std::cout << "Launching layer " << layer << std::endl;
            // one layer is one step in the BFS
            for (auto it1 = u_vertices.cbegin(); it1 != u_vertices.cend(); ++it1)
            {
               // std::cout << "Iterating over u vertex: " << *it1 << std::endl;
                std::vector<int> u_succ(nf.pull_all_near(*it1));
              //  std::cout << "Size of u_succ: " << u_succ.size() << std::endl;
                for (auto it2 = u_succ.begin(); it2 != u_succ.end(); ++it2)
                {
                  //        std::cout << "Adding " << *it2 << std::endl;
                    layered_nf.add(*it2, layer);
                    v_vertices.emplace_back(*it2);
                }
            }
            // When the above for finishes, we have progress of one half-step (from U to V) in the BFS
            u_vertices.clear();
            bool end = false;
            for (auto it = v_vertices.cbegin(); it != v_vertices.cend(); it++)
                if (v_to_u.at(*it) == constants::null_point_index)
                {
                    // we stop when a nearest exposed V vertex (from U exposed vertices) has been found
                    //     std::cout << *it << " No neighbor, exiting " << std::endl;
                    end = true;
                }
                else
                    u_vertices.emplace_back(v_to_u.at(*it));
            // When the above for finishes, we have progress of one half-step (from V to U) in the BFS
            if (end)
                return layered_nf;
            v_vertices.clear();
        }
        return layered_nf;
    }
    std::vector<int> returnMatching()
    {
        return v_to_u;
    }




    virtual ~Graph_Matching() {}

protected:

private:
    graphType* gp;
    double r;
    /** \internal \brief Given a point from V, provides its matched point in U, null_point_index() if there isn't. */
    std::vector<int> v_to_u;
    /** \internal \brief All the unmatched points in U. */
    std::unordered_set<int> unmatched_in_u;
};

#endif // GRAPH_MATCHING_H
