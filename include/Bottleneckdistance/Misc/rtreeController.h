#ifndef RTREECONTROLLER_H
#define RTREECONTROLLER_H

//#include <boost/geometry.hpp>
//#include <boost/geometry/geometries/point.hpp>
//#include <boost/geometry/geometries/box.hpp>

//#include <boost/geometry/index/rtree.hpp>

// to store queries results
#include <vector>

// just for output
#include <iostream>
#include <boost/foreach.hpp>
#include <header.hpp>

#include "rtree_graph.h"
//#include "rtree/rtree_graph.h"


class rtreeController
{

public:
    rtreeController(rtree_graph& g, double r):graph(g),r(r),rtree()
    {}
    virtual ~rtreeController() {}

    void add(int v_point_index)
    {

        auto p = this->graph.getVPoint(v_point_index);
        // std::cout << "Adding point " << v_point_index << std::endl;
        //    std::cout << "Inserted point into KD-tree : " << p.vec[0] << " , " << p.vec[1] << std::endl;
      //  std::cout << "Inserting: " << v_point_index << std::endl;
        rtree.insert(std::make_pair(p,v_point_index));
    }

    int pull_near(int u_point_index)
    {
        point u_point = graph.getUPoint(u_point_index);
        //  std::cout << "Searching for: " << u_point.vec[0] << " , " << u_point.vec[1] << std::endl;
        //   std::cout << "Radius is set to: " << r << std::endl;
        //  kd_t.build();

        point min_corner;
        point max_corner;
        for (int i = 0; i < constants::Dim; i++)
        {
            double dam = u_point[i] - r;
            double dum = u_point[i] + r;
            min_corner[i] = dam;
            max_corner[i] = dum;

        }
        //std::cout << "mincorner: first: " << min_corner[0] << std::endl;
        // std::cout << "max first: " << max_corner[0] << std::endl;

        box query_box(min_corner, max_corner);
        std::vector<value> result_s;
        //! rtree.query(bgi::covered_by(query_box), std::back_inserter(result_s));
        rtree.query(bgi::covered_by(query_box), std::back_inserter(result_s));
       //! rtree.query(bgi::intersects(query_box), std::back_inserter(result_s));
        //std::cout << "Intersection points: " << result_s.size() << std::endl;

        //   Cube_query query{u_point, r};
        // std::cout << kd_t.capacity() << std::endl;
        //!   auto neighbor = kd_t.search_any_point(query);
        if (result_s.size() == 0)
        {
            // std::cout << "No neighbors found" << std::endl;
            //    std::cout << "No neighbor found" << std::endl;
            return constants::null_point_index;
        }
        else
        {
            //    std::cout << "Searching u point: " << u_point_index << std::endl;
            //   std::cout << "Radius is set to " << r << std::endl;
            //    std::cout << "Neighbor found " << std::endl;



            //   int tmp = result_s[0].;
            //  auto point = graph.getVPoint(tmp);
            //  int idx = point.point_index;
               // std::cout << "Before removal " << rtree.size() <<std::endl;
            rtree.remove(result_s[0]);
          //   std::cout << "After removal " << rtree.size() <<std::endl;

            return result_s[0].second;
        }


    }
    std::vector<int> pull_all_near(int u_point_index)
    {
        //   std::cout << "Pull all near : " << u_point_index << std::endl;

        point u_point = graph.getUPoint(u_point_index);
        //  std::cout << "Searching for: " << u_point.vec[0] << " , " << u_point.vec[1] << std::endl;
        //   std::cout << "Radius is set to: " << r << std::endl;
        //  kd_t.build();

        point min_corner;
        point max_corner;
        for (int i = 0; i < constants::Dim; i++)
        {
            double dam = u_point[i] - r;
            double dum = u_point[i] + r;
            min_corner[i] = dam;
            max_corner[i] = dum;

        }

        box query_box(min_corner, max_corner);
        std::vector<value> result_s;
     //!   rtree.query(bgi::intersects(query_box), std::back_inserter(result_s));
        rtree.query(bgi::covered_by(query_box), std::back_inserter(result_s));





        std::vector<int> all_pull;

        for (int i = 0; i < result_s.size(); i++)
        {
          //  std::cout << "In company removal" << result_s[i].second << std::endl;
         //   std::cout << "Before removal " << rtree.size() <<std::endl;
            rtree.remove(result_s[i]);
         //   std::cout << "After removal " << rtree.size() <<std::endl;
            all_pull.push_back(result_s[i].second);
            //return result_s[0].second;


        }

        //   int last_pull = pull_near(u_point_index);
        //   while (last_pull != constants::null_point_index)
        //  {
        //      all_pull.emplace_back(last_pull);
        //      last_pull = pull_near(u_point_index);
        //   }
        return all_pull;
    }


protected:

private:
    bgi::rtree< value, bgi::quadratic<4>> rtree;
    rtree_graph& graph;
    double r;
};

#endif // RTREECONTROLLER_H
