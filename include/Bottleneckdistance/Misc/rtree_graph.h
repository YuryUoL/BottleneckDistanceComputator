#ifndef RTREE_GRAPH_H
#define RTREE_GRAPH_H


#include <cstdlib>
#include <iostream>
//#include "Internal_point.h"
#include <armadillo>
//#include "experiments/vectorTools.h"
//#include "rtreeController.h"

class rtree_graph
{
    public:


    rtree_graph(arma::mat & diag1, arma::mat & diag2):u(diag1),v(diag2)
    {
       // std::cout << "Diagram 1 : " << diag1.n_cols << " , " << diag1.n_rows << std::endl;
        //std::cout << "Diagram 2 : " << diag2.n_cols << " , " << diag2.n_rows << std::endl;

        //:out(outputpath)
        //     std::vector<std::vector<double>> dg1;
        //  std::vector<std::vector<double>> dg2;
        //   vectorTools::mat_to_std_vec(diag1,dg1);
        //   vectorTools::mat_to_std_vec(diag2,dg2);
        //  this->PrintMatrixToFile("diagram1",dg1);
        // this->PrintMatrixToFile("diagram2",dg2);
        //  this->Transformator(dg1, this->u, 0);
        //  this->Transformator(dg2, this->v, 1);
        //Dataset_Graph(dg1, dg2);

    }

//    Dataset_Graph(std::vector<std::vector<double>> & diag1, std::vector<std::vector<double>> & diag2 )
//    {
//        this->Transformator(diag1, this->u, 0);
//        this->Transformator(diag2, this->v, 1);
//
//    }
    double distance(int i, int j) const
    {
      //  std::cout << i << " , " << j << std::endl;
       // std::cout << "Compute distance" << std::endl;
        double dd = arma::norm(u.col(i)-v.col(j), "inf");
      //  std::cout << "Compute distance finnsih" << std::endl;
        return dd;
        //   double maxx = 0;

        //  for (int k = 0 ; k < constants::Dim; k++)
        //  {
        //      maxx = std::max(maxx,std::abs(u[i].get<k>()-v[j].get<k>()));

        // }
        //     std::cout << "Result: " << maxx << std::endl;
        //  return maxx;

    }
    point getVPoint(int indxV) const
    {
       // std::cout << "Get V point" << std::endl;
        point a = point();
        for (int j = 0; j < constants::Dim; j++)
        {
            a[j] = v.col(indxV)(j);
        }
        //  std::cout << "Get V point finnish" << std::endl;
        return a;


    }
    point getUPoint(int indxU) const
    {
      //   std::cout << "Get U point" << std::endl;
        point a = point();
        for (int j = 0; j < constants::Dim; j++)
        {
            a[j] = u.col(indxU)(j);
        }
      //   std::cout << "Get U point finnish " << std::endl;
        return a;
    }
    inline int size() const
    {
        return u.n_cols;
        //  return static_cast<int> (u.size() + v.size());
    }

    std::vector<double> sorted_distances()
    {
        //  out << "Distances: " << std::endl;
        //  out << "U size: " << u.size() << std::endl;
        //  out << "V size: " << v.size() << std::endl;
        std::vector<double> distances;
        distances.push_back(0.);  // for empty diagrams
        for (int u_point_index = 0; u_point_index < u.n_cols; ++u_point_index)
        {
            //distances.push_back(distance(u_point_index, corresponding_point_in_v(u_point_index)));
            for (int v_point_index = 0; v_point_index < v.n_cols; ++v_point_index)
            {
                double dtmp = distance(u_point_index, v_point_index);
                //  out << u_point_index << " , " << v_point_index << " : " << dtmp << std::endl;
                distances.push_back(dtmp);
            }
        }
//#ifdef GUDHI_USE_TBB
//       tbb::parallel_sort(distances.begin(), distances.end());
//#else
        std::sort(distances.begin(), distances.end());
//#endif
        return distances;
    }




    virtual ~rtree_graph()
    {
        //out.close();
    }

protected:
//    void Transformator(std::vector<std::vector<double>> & diag1, std::vector<Internal_point> & output, int color)
//    {
//        // out << "Transformator activated with: " << diag1.size() << std::endl;
//        for (int i = 0; i < diag1.size() ; i++ )
//        {
//            //double* a = &diag1[i][0];
//            //   std::cout << "*** Adding point " << diag1[i][0] << ", "  << diag1[i][1] << std::endl;
//            output.push_back(Internal_point(diag1[i], i, color));
//        }





private:
    arma::mat u;
    arma::mat v;
};

#endif // RTREE_GRAPH_H
