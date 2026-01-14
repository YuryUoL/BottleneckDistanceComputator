#ifndef BOOSTBOTTLENECK_H
#define BOOSTBOTTLENECK_H

#include <iostream>
#include <string>
#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <cassert>
#include <boost/graph/max_cardinality_matching.hpp>
#include <armadillo>

using namespace boost;

typedef adjacency_list<hash_setS, vecS, undirectedS > my_graph;


class BoostBottleneck
{
public:
    BoostBottleneck() {}

    void modifyGraph(my_graph & G, int thresh, int oldthresh, int sizeCloud, std::vector<std::pair<double, std::vector<int>>> & matBoy )
    {
        if (thresh > oldthresh)
        {
            for (int i = oldthresh+1; i <= thresh; i++)
            {
                int first = matBoy[i].second[0];
                int second = matBoy[i].second[1] + sizeCloud;
                boost::add_edge(first,second, G);
            }
        }
        else
        {
            for (int i = oldthresh; i > thresh; i--)
            {
                int first = matBoy[i].second[0];
                int second = matBoy[i].second[1] + sizeCloud;
                boost::remove_edge(first,second, G);
            }



        }


    }

    double bottleneck_between_clouds(arma::mat & A, arma::mat & B, std::vector<int> & matching)
    {


        std::vector<std::vector<double>> distMatrix = std::vector<std::vector<double>>(A.n_cols, std::vector<double>(B.n_cols));
        //!distMatrix.resize(A.n_cols, B.n_cols);
        for (int i = 0 ; i < A.n_cols ; i++)
        {
            for (int j = 0; j < B.n_cols ; j++)
            {
                //arma::norm(u.col(i)-v.col(j), "inf");
                distMatrix[i][j] = arma::norm(A.col(i)-B.col(j), "inf");
            }


        }

//        for (int i = 0; i < distMatrix.size(); i++)
//        {
//            for (int j = 0; j < distMatrix.size(); j++)
//            {
//
//                std::cout << distMatrix[i][j] << ", ";
//
//            }
//
//            std::cout << std::endl;
//
//
//        }

        return  bottleneck_on_bipartie_graph(distMatrix, matching);

    }

    double bottleneck_on_bipartie_graph(std::vector<std::vector<double>> & matrix, std::vector<int> & matching)
    {

        std::vector<std::pair<double, std::vector<int>>> matBoy;
        // std::vector<std::pair<double, std::vector<int>> matBoy;
        for (int i = 0; i < matrix.size(); i++)
        {
            for (int j = 0; j < matrix[0].size(); j++)
            {
                matBoy.push_back(std::make_pair(matrix[i][j], std::vector<int> {i,j}));
            }

        }



        std::ranges::sort(matBoy, [](const std::pair<double, std::vector<int>> &a, const  std::pair<double, std::vector<int>>  &b)
        {
            return a.first < b.first;
        });

//
   //     std::cout << "Vector is: " << std::endl;
//
//        for (int i = 0; i < matBoy.size(); i++)
//        {
//            std::cout << matBoy[i].first << ", ";
//
//        }

      //  std::cout << std::endl;


        int sizeCloud = matrix.size();
        my_graph G(2*sizeCloud);
        int oldThresh = -1;
        int lower_bound_i = 0;
        int upper_bound_i = matBoy.size()-1;

        double alpha = std::pow(2*sizeCloud, 1. / 5.);

        int biggest_unperfect = -1;

        std::vector< graph_traits< my_graph >::vertex_descriptor > mateSafe(2*sizeCloud);
        my_graph Gsafe;

      //  std::cout << "Size cloud is " << sizeCloud << std::endl;

        while (lower_bound_i != upper_bound_i)
        {

            std::vector< graph_traits< my_graph >::vertex_descriptor> mate(2*sizeCloud);

            //!m =  std::vector< graph_traits< my_graph >::vertex_descriptor>(2*sizeCloud);

            long step = lower_bound_i + static_cast<long> ((upper_bound_i - lower_bound_i - 1) / alpha);
            modifyGraph(G, step, oldThresh, sizeCloud, matBoy );

            double dudu =  matBoy[step].first;
            //std::cout << "Currently " << lower_bound_i << " , " << upper_bound_i << " step is  " << step << " value: " << dudu <<  std::endl;
          //  printGraph(G);

            oldThresh= step;

            //checked_edmonds_maximum_cardinality_matching(G, &mate[0]);
            edmonds_maximum_cardinality_matching(G, &mate[0]);
            if (matching_size(G, &mate[0]) == sizeCloud)
            {

                mateSafe = mate;
                Gsafe = G;
                upper_bound_i = step;
                //std::cout << "Perfect matching! " << std::endl;

            }
            else
            {
                //std::cout << "not perfect matching! " << std::endl;
                lower_bound_i = step + 1;

            }

        }

        //!return 0;

      //  std::cout << "lower bound is : " << lower_bound_i << std::endl;
        matching = std::vector<int>(sizeCloud);
       // std::cout << "The matching is:" << std::endl;

        graph_traits< my_graph >::vertex_iterator vi, vi_end;
        for (boost::tie(vi, vi_end) = vertices(Gsafe); vi != vi_end; ++vi)
        {
            if (mateSafe[*vi] != graph_traits< my_graph >::null_vertex()
                    && *vi < mateSafe[*vi])
                    {
                        matching[*vi] = mateSafe[*vi];

                    }
             //   std::cout << "{" << *vi << ", " << mateSafe[*vi] << "}" << std::endl;

                }

     //   std::cout << std::endl;


        for (int i = 0; i < matching.size() ;i++)
        {
            matching[i] = matching[i] - matrix.size();

        }

        return matBoy[lower_bound_i].first;

    }

//
//    double computeBottleNeckMatching(std::vector<std::vector<double>> & matrix)
//    {
//
//        std::vector<std::pair<double, std::vector<int>>> matBoy;
//        // std::vector<std::pair<double, std::vector<int>> matBoy;
//        for (int i = 0; i < matrix.size(); i++)
//        {
//            for (int j = 0; j < matrix[0].size(); j++)
//            {
//                matBoy.push_back(std::make_pair(matrix[i][j], std::vector<int> {i,j}));
//            }
//
//        }
//        std::ranges::sort(matBoy, [](const std::pair<double, std::vector<int>> &a, const  std::pair<double, std::vector<int>>  &b)
//        {
//            return a.first < b.first;
//        });
//
//        int sizeCloud = matrix.size();
//        my_graph G(2*sizeCloud);
//        int oldThresh = -1;
//        int lower_bound_i = 0;
//        int upper_bound_i = matBoy.size()-1;
//
//        double alpha = std::pow(2*sizeCloud, 1. / 5.);
//
//        int biggest_unperfect = -1;
//
//        std::vector< graph_traits< my_graph >::vertex_descriptor > mate(2*sizeCloud);
//        while (lower_bound_i != upper_bound_i)
//        {
//            mate =  std::vector< graph_traits< my_graph >::vertex_descriptor>(2*sizeCloud);
//
//            long step = lower_bound_i + static_cast<long> ((upper_bound_i - lower_bound_i - 1) / alpha);
//            BoostBottleneck::modifyGraph(G, step, oldThresh, sizeCloud, matBoy );
//            oldThresh= step;
//
//            checked_edmonds_maximum_cardinality_matching(G, &mate[0]);
//            if (matching_size(G, &mate[0]) == sizeCloud)
//            {
//
//
//                upper_bound_i = step;
//
//            }
//            else
//            {
//
//                lower_bound_i = step + 1;
//
//            }
//
//        }
//
//        //!return 0;
//
//        return matBoy[lower_bound_i].first;
//
//    }



protected:

private:
};

#endif // BOOSTBOTTLENECK_H
