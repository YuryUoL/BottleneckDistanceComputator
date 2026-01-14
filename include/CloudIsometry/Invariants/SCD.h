#ifndef SCD_H
#define SCD_H

#include <armadillo>
#include "CloudIsometry/Nested_Base.h"
class SCD
{
public:
    SCD():constC(0.43) {}
    SCD(double dt):constC(dt){}
    virtual ~SCD() {}


    int NumElem(std::pair<arma::mat, std::vector<std::vector<bool>>> & cloudAunc)
    {
        return cloudAunc.first.n_cols;

    }
    int Dimension(std::pair<arma::mat, std::vector<std::vector<bool>>> & cloudAunc)
    {

        return cloudAunc.first.n_rows;

    }

    double initilize(std::pair<arma::mat, std::vector<std::vector<bool>>> & cloudAunc, std::pair<arma::mat, std::vector<std::vector<bool>>> & cloudA)
    {

        arma::mat cloudOUT;
        double sizeofcloudin = cloudAunc.first.n_cols;

        arma::vec centroid(cloudAunc.first.n_rows);
        double radius = 0;
        for (int i = 0; i < cloudAunc.first.n_cols; i++)
        {
            arma::vec tmpZ =  cloudAunc.first.col(i);
            centroid = centroid + tmpZ;
        }
        double division = -1/sizeofcloudin;

        centroid = division * centroid;

        cloudOUT = cloudAunc.first;
        for (int i = 0; i < cloudOUT.n_cols; i++)
        {

            arma::vec blablaz = cloudAunc.first.col(i);
            arma::vec tmpd = blablaz + centroid;
            cloudOUT.col(i) = tmpd;

        }

        cloudA.second = cloudAunc.second;
        cloudA.first = cloudOUT;

     //   std::cout << "Initial cloud : " << std::endl;

     //   std::cout << cloudAunc.first << std::endl;

      //  std::cout << "Transformed cloud: " << std::endl;

      //  std::cout << cloudA.first << std::endl;

        return radius;
    }

    void selectMatrix(std::vector<int> & indices, arma::mat & input, arma::mat & output)
    {
        output.resize(input.n_rows, indices.size());
        for (int i = 0; i < indices.size(); i++)
        {
            output.col(i) = input.col(indices[i]);
        }

    }
    void theOthetMatrix(std::vector<int> & indices, std::vector<int> & counterCombos, arma::mat & input, arma::mat & output)
    {
        std::set<int> s(indices.begin(), indices.end());
        output.resize(input.n_rows, input.n_cols - indices.size());
        int j = 0;
        for (int i = 0 ; i < input.n_cols ; i++)
        {
            if (s.find(i) == s.end())
            {
                counterCombos.push_back(i);
                output.col(j) = input.col(i);
                j++;
            }
        }


    }



    void ConstructCloudInvariantCombo(std::pair<arma::mat, std::vector<std::vector<bool>>> & cloudA,
                                      std::vector<std::vector<int>> & combos, std::vector<arma::mat> &  OCDCloudCollection, std::vector<arma::mat> & distMatrixCollection)
    {
        for (int i = 0; i < combos.size(); i++)
        {
            //  std::cout << "Staritng combo " << i << std::endl;

            arma::mat selectedPoints;
            arma::mat notSelectedpoints;

            //CrystData::ConstructSRCloud(cloud,basisCollection,combos[i],SrCloud);
            // CrystData::ConstructDistanceMatrix(basisCollection,DistMat);
            std::vector<int> counterCombos;
          //   std::cout << "Selecting matrices" << std::endl;
            this->selectMatrix(combos[i],cloudA.first,selectedPoints);
            this->theOthetMatrix(combos[i],counterCombos,cloudA.first,notSelectedpoints);
            arma::mat ocdCloud;
            arma::mat distMat;
           //  std::cout << "Matrix selection succesful" << std::endl;
            //  off << "Combo: " << i << std::endl;
            this->ProduceSCDMatrixM(selectedPoints,notSelectedpoints,ocdCloud, cloudA.second, combos[i], counterCombos, cloudA.first.n_cols);
           //  std::cout << "Produce OCD succeful" << std::endl;
            this->ProduceDistanceMat(selectedPoints,distMat, cloudA.second, combos[i], cloudA.first.n_cols);
          //  std::cout << "Produce distMat succesful" << std::endl;


            OCDCloudCollection.push_back(ocdCloud);
            distMatrixCollection.push_back(distMat);
        }

    }


    void ProduceSCDMatrixM(arma::mat & selectedPoints, arma::mat & notselectedpoints, arma::mat & ocdCloud,
                           std::vector<std::vector<bool>> & graph, std::vector<int> & combos, std::vector<int> & counterCombos, int NumPoints)
    {
//        std::cout << "Selectedpoints n cols: " << selectedPoints.n_cols << std::endl;
//        std::cout << "Not selected points n cols: " << notselectedpoints.n_cols << std::endl;
//        std::cout << "Combos: " << std::endl;
//        for (int i = 0; i < combos.size(); i++)
//        {
//            std::cout << combos[i] << ", ";
//        }
//        std::cout << std::endl;
//
//        std::cout << "Countercombos: " << std::endl;
//        for (int i = 0; i < counterCombos.size(); i++)
//        {
//            std::cout << counterCombos[i] << ", ";
//        }
//        std::cout << std::endl;

        ocdCloud.resize(selectedPoints.n_cols+2,notselectedpoints.n_cols);
        arma::vec origin(selectedPoints.n_rows, arma::fill::zeros);

        for (int j = 0; j < notselectedpoints.n_cols; j++)
        {
            arma::mat finalDeterminant(notselectedpoints.n_rows, selectedPoints.n_rows);
            for (int i = 0 ; i < selectedPoints.n_cols; i++)
            {

                double dist = arma::norm(notselectedpoints.col(j)-selectedPoints.col(i));
                finalDeterminant.col(i) = notselectedpoints.col(j)-selectedPoints.col(i);
                if (graph[combos[i]][counterCombos[j]] == false)
                {
                    dist = -1 * dist;
                }
                ocdCloud(i,j) = dist;
            }

            double tmpdist = arma::norm(origin-notselectedpoints.col(j));

            if (graph[NumPoints][counterCombos[j]] == false)
            {
                tmpdist = -1 * tmpdist;
            }

            ocdCloud(selectedPoints.n_cols,j) = tmpdist;


            finalDeterminant.col(selectedPoints.n_cols) = notselectedpoints.col(j);


            double dett = arma::det(finalDeterminant);
            double answer = 0;

            if (dett < 0)
            {
                answer = -1;
            }
            if (dett > 0)
            {
                answer = +1;
            }

            arma::mat volumeObj(selectedPoints.n_rows, selectedPoints.n_cols + 2, arma::fill::zeros);

            for (int a = 0; a < selectedPoints.n_cols; a++)
            {
                volumeObj.col(a) = selectedPoints.col(a);
            }

            volumeObj.col(selectedPoints.n_cols) = notselectedpoints.col(j);

            double strengthh = ComputeStrength(volumeObj);

//            if (strengthh < strengthCut)
//            {
//
//                strengthh = 0;
//            }
            ocdCloud(selectedPoints.n_cols+1,j)=answer/this->constC * strengthh;
        }
    }
    double ComputeStrength(arma::mat & input)
    {

        arma::mat caylerMenger(input.n_cols+1, input.n_cols+1);
        caylerMenger(0,0) = 0;
        for (int i = 1; i < caylerMenger.n_cols; i++)
        {
            caylerMenger(caylerMenger.n_rows-1,i) = 1;
            caylerMenger(i,caylerMenger.n_cols-1) = 1;
        }
        caylerMenger(input.n_cols, input.n_cols) = 0;
        double halfsumofdistances = 0;
        for (int i = 0; i < caylerMenger.n_cols-1; i++)
        {
            for (int j = 0; j < caylerMenger.n_rows-1; j++)
            {
                double tmpdist = arma::norm(input.col(i)-input.col(j),2);
                caylerMenger(j,i) = tmpdist * tmpdist;
                halfsumofdistances = halfsumofdistances + tmpdist;
            }
        }
        halfsumofdistances  = halfsumofdistances / 4;


        double dt = arma::det(caylerMenger);
        double multiplier = (double) pow(-1, input.n_cols+1) / ((double) pow(2,input.n_cols) * (double) pow(tgamma(input.n_cols+1),2.0));
        double tmpdt = multiplier * dt;



        return tmpdt * tmpdt / (double) pow(halfsumofdistances,2*(input.n_cols-1)-1);

    }

    void ProduceDistanceMat(arma::mat & selectedPoints, arma::mat & distMat, std::vector<std::vector<bool>> & graph, std::vector<int> & combos, int numPoints)
    {
        distMat.resize(selectedPoints.n_cols + 1, selectedPoints.n_cols + 1);
        arma::vec origin(selectedPoints.n_rows, arma::fill::zeros);
        for (int i  = 0; i < selectedPoints.n_cols; i++)
        {
            double dist1 = arma::norm(origin-selectedPoints.col(i),2);
            if (graph[numPoints][combos[i]] == false)
            {
                dist1 = dist1 * -1;


            }
            distMat(selectedPoints.n_cols,i) = dist1;
            distMat(i,selectedPoints.n_cols) = dist1;
            for (int j = 0; j < selectedPoints.n_cols; j++)
            {
                double dist = arma::norm(selectedPoints.col(i)-selectedPoints.col(j),2);
                if (graph[combos[j]][combos[i]] == false)
                {
                    dist = dist* -1;
                }
                distMat(j,i) = dist;
            }
        }
    }


    //! inner_computator.initilize(cloudAunc, cloudA);

//    void initilize(std::pair<arma::mat, std::vector<std::vector<double>>> & cloudAunc, std::pair<arma::mat, std::vector<std::vector<double>>> & cloudA)
//    {
//
//
//
//
//    }




protected:

private:
    double constC;
};

#endif // SCD_H
