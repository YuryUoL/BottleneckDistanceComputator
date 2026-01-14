#ifndef BD_GENERAL_INTERFACE_H
#define BD_GENERAL_INTERFACE_H

template<typename cloud_compute_engine_type, typename graph_compute_engine_type>
class bd_general_interface
{
public:
    bd_general_interface():ct(cloud_compute_engine_type()),gt(graph_compute_engine_type()) {}
    virtual ~bd_general_interface() {}


    double computeLinfty(arma::mat & dataA, arma::mat & dataB)
    {

        double Maxx = 0;
        //     return arma::norm(dataA - dataB,"inf");
        for (int i = 0; i < dataA.n_rows; i++)
        {
            for (int j = 0; j < dataA.n_cols; j++)
            {

                double distt = std::abs(dataA(i,j) - dataB(i,j));
                Maxx = std::max(distt, Maxx);
            }
        }
        return Maxx;

    }

    double computeBetweenInvariantDistancePair(arma::mat & distA, arma::mat & distB, arma::mat & ocdA, arma::mat & ocdB, double & d, std::vector<int> & matching)
    {

        //!  double bda = ct.bottleneck_between_clouds(distA, distB, matching);
        //! double bda = arma::norm(distA-distB,"inf");
        double bda = computeLinfty(distA, distB);
        //!    matching.clear();

//         std::vector<int> matching;
         cloud_compute_engine_type ct_tmp;

        double bdb = ct_tmp.bottleneck_between_clouds(ocdA, ocdB, matching);
        d = bdb;
        return std::max(bda, bdb);

    }

    double computeBetweenInvariantDistancePair(arma::mat & distA, arma::mat & distB, arma::mat & ocdA, arma::mat & ocdB, double & d, double & d2, std::vector<int> & matching)
    {

        //!  double bda = ct.bottleneck_between_clouds(distA, distB, matching);
        //! double bda = arma::norm(distA-distB,"inf");
        double bda = computeLinfty(distA, distB);
        //!    matching.clear();

//         std::vector<int> matching;
         cloud_compute_engine_type ct_tmp;

        double bdb = ct_tmp.bottleneck_between_clouds(ocdA, ocdB, matching);
        d = bdb;
        d2 = bda;
        return std::max(bda, bdb);

    }


    double computeBetweenCloudPair(arma::mat & ocdA, arma::mat & ocdB, std::vector<int> & matching)
    {
        double bdb = ct.bottleneck_between_clouds(ocdA, ocdB, matching);
        return bdb;


    }

      double computeBetweenInvariantDistancePair(arma::mat & distA, arma::mat & distB, arma::mat & ocdA, arma::mat & ocdB)
    {
        std::vector<int> matching;
         cloud_compute_engine_type ct_tmp;

        //!  double bda = ct.bottleneck_between_clouds(distA, distB, matching);
        //! double bda = arma::norm(distA-distB,"inf");
        double bda = computeLinfty(distA, distB);
        //!    matching.clear();
        double bdb = ct_tmp.bottleneck_between_clouds(ocdA, ocdB, matching);
        return std::max(bda, bdb);

    }

    double computeForGraph(std::vector<std::vector<double>> & distMat,std::vector<int> & matching)
    {

//!           std::vector<int> matching;

//        std::cout << "Distance matrix was : " << std::endl;
//
//        for (int i = 0; i < distMat.size(); i++)
//        {
//            for (int j = 0; j < distMat.size(); j++)
//            {
//                std::cout << distMat[i][j] << ", ";
//            }
//            std::cout << std::endl;
//
//
//        }


        double bda = gt.bottleneck_on_bipartie_graph(distMat, matching);

    //   std::cout << "Returnable bda is : " << bda << std::endl;

        return bda;

    }

protected:

private:
    cloud_compute_engine_type ct;
    graph_compute_engine_type gt;

};

#endif // BD_GENERAL_INTERFACE_H
