#ifndef BOTTLENECKCALCULATOR_H
#define BOTTLENECKCALCULATOR_H

#include "Graph_Matching.h"

template<typename graphType,typename treeController>
class BottleneckCalculator
{
public:
    BottleneckCalculator() {}




    double bottleneck_between_clouds(arma::mat & A, arma::mat & B, std::vector<int> & matching)
    {
        graphType g(A,B);
        return bottleneck_distance_exact_with_matching(g, matching);
    }

    double bottleneck_on_bipartie_graph(std::vector<std::vector<double>> & distMat, std::vector<int> & matching)
    {
        graphType g(distMat);
        return bottleneck_distance_exact_with_matching(g, matching);
    }

    void truncateMatching(std::vector<int> & matchingIn, std::vector<int> & matchingOut)
    {
      //  std::vector<int> tmpd = matching;
      matchingOut = std::vector<int>(matchingIn.size());
        for (int i = 0; i < matchingIn.size(); i++)
        {
            if (matchingIn[i] < matchingIn.size())
            {
            matchingOut[matchingIn[i]] = i;
            }

        }

    }


    double bottleneck_distance_exact_with_matching(graphType & gg, std::vector<int> & matching)
    {
        graphType g = gg;
        std::vector<int> matchingTmp;
        //     std::cout << "Starting method of bottleneck computation" << std::endl;
        std::vector<double> sd = g.sorted_distances();
        long lower_bound_i = 0;
        long upper_bound_i = sd.size() - 1;
        const double alpha = std::pow(g.size(), 1. / 5.);
        Graph_Matching<graphType, treeController> m(g);
        Graph_Matching<graphType, treeController> biggest_unperfect(g);
//        std::cout << "Sd size: " << sd.size() << std::endl;
//            for (int i = 0; i < sd.size(); i++)
//           {
//           std::cout << sd[i] << ", ";
//        }
//            std::cout << std::endl;
        //  std::cout << "Starting bottleneck computator " << std::endl;

        while (lower_bound_i != upper_bound_i)
        {

            long step = lower_bound_i + static_cast<long> ((upper_bound_i - lower_bound_i - 1) / alpha);
            //     std::cout << "Current step: " << step << std::endl;
            m.set_r(sd.at(step));
            //  std::cout << "Radius set to : " << sd.at(step) << std::endl;
            while (m.multi_augment()) {}  // compute a maximum matching (in the graph corresponding to the current r)
            if (m.perfect())
            {
                m = biggest_unperfect;
                upper_bound_i = step;
                //   std::cout << "Perfect found!" << std::endl;
                //    std::cout << "Upper bound is now " << upper_bound_i << std::endl;
            }
            else
            {
                biggest_unperfect = m;
                lower_bound_i = step + 1;
                //  std::cout << "UNPERFECT RESULT" << std::endl;
                // std::cout << "Low bound is now: " << lower_bound_i << std::endl;
            }
        }
        Graph_Matching<graphType, treeController> mm(g);
       //! std::cout << " Starting most important run " << std::endl;
        mm.set_r(sd.at(lower_bound_i));
        while (mm.multi_augment()) {};
        matchingTmp = mm.returnMatching();


        //   double maxxDist = -1;

        //  std::cout << "Final: " << sd.at(lower_bound_i) << std::endl;

        // for (int i = 0; i < matching.size(); i++)
        //   {
        //    maxxDist = std::max(maxxDist,g.distance(i, matching[i]));

        // }
       // truncateMatching(matchingTmp);
 //!       matching = matchingTmp;

        truncateMatching( matchingTmp, matching);

        return sd.at(lower_bound_i);
        // return maxxDist;
    }



protected:

private:
};

#endif // BOTTLENECKCALCULATOR_H
