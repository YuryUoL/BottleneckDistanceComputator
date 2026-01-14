#ifndef PBI_ORIGINAL_H
#define PBI_ORIGINAL_H

#include "vectorTools.h"
#include <set>

class PBI_Original
{
public:
    PBI_Original() {}

    int NumElem(std::vector<std::vector<double>> & cloudAunc)
    {
        return cloudAunc.size();

    }
    int Dimension(std::vector<std::vector<double>> & cloudAunc)
    {

        return cloudAunc[0].size();

    }

    //! double initilize(arma::mat & cloudAunc, arma::mat & cloudA)


    double initilize(std::vector<std::vector<double>> & cloudIN, std::vector<std::vector<double>> & cloudOUT)
    {
        double sizeofcloudin = cloudIN.size();
        //std::cout << "Size of cloud in: " << sizeofcloudin << std::endl;
        std::vector<double> centroid(cloudIN[0].size());
        double radius = 0;
        for (int i = 0; i < cloudIN.size(); i++)
        {
            vectorTools::vectorAddition(centroid, cloudIN[i]);
        }
        double division = -1/sizeofcloudin;
        //  std::cout << "Division size: " << division << std::endl;
        // std::cout << "Centroid: "  << centroid[0] << ", " << centroid[1] << ", " << centroid[2] << std::endl;
        vectorTools::vectorScalarMultiplication(centroid,division);
        //   std::cout << "Centroid: "  << centroid[0] << ", " << centroid[1] << ", " << centroid[2] << std::endl;
        cloudOUT = cloudIN;
        for (int i = 0; i < cloudOUT.size(); i++)
        {
            //     std::cout << "cloudOut: "  << cloudOUT[i][0] << ", " << cloudOUT[i][1] << ", " << cloudOUT[i][2] << std::endl;
            vectorTools::vectorAddition(cloudOUT[i], centroid);
            //   std::cout << "cloudOut: (After addition) "  << cloudOUT[i][0] << ", " << cloudOUT[i][1] << ", " << cloudOUT[i][2] << std::endl;
            radius = std::max(radius, sqrt(vectorTools::SquaredNorm(cloudOUT[i])));
        }
        return radius;
    }


    void computeBasisMatrix(std::vector<int> & indices,std::vector<std::vector<double>> & cloud, std::vector<std::vector<double>> & basisCollection)
    {
        std::vector<double> c_P(cloud[indices[0]].size());
        // std::vector<std::vector<double>> basisCollection;
        for (int i = 0; i < indices.size(); i++)
        {
            basisCollection.push_back(cloud[indices[i]]);
        }
        vectorTools::multiDcrossProduct(basisCollection, c_P);
        //vectorTools::crossProduct(cloud[indices[0]], cloud[indices[1]], c_P);

        basisCollection.push_back(c_P);



    }

    void ConstructDistanceMatrix(std::vector<std::vector<double>> & basisCollection, std::vector<std::vector<double>> & matrix)
    {
        //!! Assume that dimension = 3
        matrix = std::vector<std::vector<double>>(basisCollection.size(), std::vector<double>(basisCollection.size()));

        for (int i = 0; i < basisCollection.size(); i++)
        {
            for (int j = 0; j < basisCollection.size(); j++)
            {
                std::vector<double> vecU = basisCollection[j];
                vectorTools::vectorScalarMultiplication(vecU,-1);
                std::vector<double> vecV = basisCollection[i];
                vectorTools::vectorAddition(vecV, vecU);
                double norm;
                matrix[i][j] = vectorTools::SquaredNorm(vecV);
            }
        }
    }

    void ConstructSRCloud( std::vector<std::vector<double>> & cloud,std::vector<std::vector<double>> & basis, std::vector<int> & indices,std::vector<std::vector<double>> & SrCloud)
    {
        std::vector<double> c_P(cloud[0].size());
        std::set<int> s(indices.begin(), indices.end());

        // std::cout << "First: " << cloud[i][0] << ", "<< cloud[i][1] << ", " << cloud[i][2] << std::endl;
        // std::cout << "Second: " << cloud[j][0] << ", "<< cloud[j][1] << ", " << cloud[j][2] << std::endl;
//        vectorTools::crossProduct(cloud[i], cloud[j], c_P);
        // std::cout << "Final: " << c_P[0] << ", "<< c_P[1] << ", " << c_P[2] << std::endl;
        for (int k = 0; k < cloud.size() ; k++)
        {
            if (s.find(k) == s.end())
            {
                std::vector<double> newVector;
                for (int i = 0; i < basis.size(); i++)
                {
                    newVector.push_back(vectorTools::dotProduct(cloud[k],basis[i]));
                }

                SrCloud.push_back(newVector);
            }
        }

        //CrystData::PrintMatrix(SrCloud);

        // std::cout << "SrCloudMat dimensions: " << SrCloud.size() << " , " << SrCloud[i].size() << std::endl;
    }

    void ConstructCloudInvariantCombo(std::vector<std::vector<double>> & cloud, std::vector<std::vector<int>> & combos, std::vector<arma::mat> & SrCloudCollection,
                                      std::vector<arma::mat> & distMatrix)
    {

        for (int i = 0; i < combos.size(); i++)
        {
            std::vector<std::vector<double>> basisCollection;
            computeBasisMatrix(combos[i],cloud,basisCollection);
            std::vector<std::vector<double>> SrCloud;
            std::vector<std::vector<double>> DistMat;
            ConstructSRCloud(cloud,basisCollection,combos[i],SrCloud);
            ConstructDistanceMatrix(basisCollection,DistMat);

            //CrystData::ConstructSRCloud(i,j,cloud,SrCloud);
            //    CrystData::ConstructDistanceMatrix(indices,cloud,DistMat);
            arma::mat A;
            arma::mat B;
            vectorTools::vec_to_mat(SrCloud, A);
            vectorTools::vec_to_mat(DistMat, B);

            SrCloudCollection.push_back(A);
            distMatrix.push_back(B);
//            if (i == 1)
//            {
//                std::cout << "Special prrint" << std::endl;
//
//
//                std::cout << "Anotehrrr" << std::endl;
//
//                PrintMatrix(SrCloudCollection[0]);
//            }




        }
    }

        virtual ~PBI_Original() {}

protected:

private:
    };

#endif // PBI_ORIGINAL_H
