#ifndef BIPARTIE_GRAPH_H
#define BIPARTIE_GRAPH_H

//#import "experiments/CrystData.h"


class Bipartie_Graph
{
public:

    virtual ~Bipartie_Graph() {}

    Bipartie_Graph(std::vector<std::vector<double>> & distMatrixx):distMatrix(distMatrixx)
    {


    }

    double computeLinf(arma::vec p, arma::vec q)
    {
        double maxxVal = 0;
        for (int i = 0; i < p.size(); i++)
        {
            maxxVal = std::max(maxxVal, abs(p(i)-q(i)));


        }
        return maxxVal;


    }

    Bipartie_Graph(arma::mat & A, arma::mat & B)
    {
        distMatrix = std::vector<std::vector<double>>(A.n_cols, std::vector<double>(B.n_cols));
        //!distMatrix.resize(A.n_cols, B.n_cols);
        for (int i = 0 ; i < A.n_cols ; i++)
        {
            for (int j = 0; j < B.n_cols ; j++)
            {
                //arma::norm(u.col(i)-v.col(j), "inf");
               //! distMatrix[i][j] = arma::norm(A.col(i)-B.col(j), "inf");
                distMatrix[i][j] = computeLinf(A.col(i), B.col(j));
            }


        }



    }
    double distance(int i, int j) const
    {
        return distMatrix[i][j];
//        double maxx = 0;
////        for (int k = 0 ; k < constants::Dim; k++)
//        {
//            maxx = std::max(maxx,std::abs(u[i].coord(k)-v[j].coord(k)));
//
//        }
//        return maxx;

    }
//    Internal_point getVPoint(int indxV) const
//    {
//        return v[indxV];
//    }
//    Internal_point getUPoint(int indxU) const
//    {
//        return u[indxU];

    inline int size() const
    {
        return distMatrix[0].size();
        //  return static_cast<int> (u.size() + v.size());
    }

    template<typename T>
    static void PrintVector(std::vector<T> & vec)
    {
        for (int i = 0; i < vec.size(); i++)
        {
            std::cout << vec[i] << ", ";


        }
        std::cout << std::endl;


    }

    inline std::vector<double> sorted_distances() const
    {
        double errorEpsilon = 0.001;
        std::vector<double> distances;
        distances.push_back(0.);  // for empty diagrams
        for (int i = 0; i < distMatrix.size(); i++)
        {
            //distances.push_back(distance(u_point_index, corresponding_point_in_v(u_point_index)));
            for (int j = 0; j < distMatrix[0].size(); j++)
            {
                distances.push_back(distance(i,j));
                //    distances.push_back(distance(i,j) + errorEpsilon);
                //    distances.push_back(distance(i,j) - errorEpsilon);
            }
        }

//#ifdef GUDHI_USE_TBB
//       tbb::parallel_sort(distances.begin(), distances.end());
//#else
        std::sort(distances.begin(), distances.end());
        //! Bipartie_Graph::PrintVector<double>(distances);
//#endif
        return distances;
    }

protected:

private:
    std::vector<std::vector<double>> distMatrix;

};

#endif // BIPARTIE_GRAPH_H
