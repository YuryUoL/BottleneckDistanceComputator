#ifndef NCDGRABBER_H
#define NCDGRABBER_H


#include <iostream>
#include <fstream>
#include <vector>
#include <armadillo>
#include <sstream>

class NCDGrabber
{
public:
    NCDGrabber(std::string badboy)
    {

        rads.push_back(0);
        rads.push_back(0);
        ocd.push_back(std::vector<arma::mat>());
        dists.push_back(std::vector<arma::mat>());
        ocd.push_back(std::vector<arma::mat>());
        dists.push_back(std::vector<arma::mat>());
        constC = sqrt(3);
        aConst = 14;
        furthestPoint.push_back(0);
        furthestPoint.push_back(0);

        signs.push_back(std::vector<double>());
        signs.push_back(std::vector<double>());



    }
//std::pair<arma::mat, std::vector<std::vector<bool>>>

    double WhatWhat()
    {

    return -1;
    }

     void resInitilize(std::vector<std::vector<int>> & res, int n, int k, std::vector<int> current, boost::unordered::unordered_map<std::vector<int>,int> & umap, int & recnum)
    {
        std::unordered_set<int> duplicateCheck(current.begin(), current.end());


        recnum = recnum + 1;

        if (current.size() == k)
        {

            res.push_back(current);
            //  std::cout << "Res size is now: " << res.size() << std::endl;
            umap[current] = res.size()-1;
        }
        else
        {
            for (int i = 0; i < n; i++)
            {
                if(duplicateCheck.find(i) == duplicateCheck.end())
                {
                    std::vector<int> currentNew = current;
                    currentNew.push_back(i);
                    resInitilize(res,n,k,currentNew,umap, recnum);
                }
            }
        }
    }

    double GrabTheThings(std::pair<arma::mat, std::vector<std::vector<bool>>>& cloudAunc,  int k, std::string pathOut)
    {

        using Clock = std::chrono::steady_clock;
        using std::chrono::time_point;
        using std::chrono::duration_cast;
        //!using std::chrono::milliseconds;
        using std::chrono::microseconds;
        using namespace std::literals::chrono_literals;
//        using std::this_thread::sleep_for;
        using namespace std;
        boost::unordered::unordered_map<std::vector<int>,int> umap;
        std::vector<std::vector<int>> res;
        std::pair<arma::mat, std::vector<std::vector<bool>>> cloudA;
       // obj cloudB;
     //   invariantType tmp;
     //   inner_computator = tmp;

//        std::cout << "Startg 123" << std::endl;

        std::cout << "One " << std::endl;



       this->initilize(cloudAunc, cloudA,0);

         std::cout << "Two " << std::endl;

//       inner_computator.initilize(cloudBunc, cloudB,1);

       //  std::cout << "Three " << std::endl;

        int n = this->NumElem(cloudAunc);
        std::vector<int> current;
        int recnum = 0;
//       std::cout << "Resinitilize : " << n << " , " << k << std::endl;
       resInitilize(res, n, k, current, umap, recnum);

     std::cout << "Four " << std::endl;

    this->ConstructCloudInvariantCombo(cloudA, res,0);

            std::cout << "DDDD" << std::endl;

            std::cout << "Attempting to save to: " << pathOut << std::endl;

          //  std::cout << ocd[0].size() << std::endl;
         //   std::cout << dists[0].size() << std::endl;
        //    std::cout << signs[0].size() << std::endl;
          //  std::cout << res.size() << std::endl;


this->write_csv(pathOut, res, ocd[0] , dists[0], signs[0] );
                 //!  const std::vector<arma::mat>& NCD1,
                 //!  const std::vector<arma::mat>& NCD2,
                  //! const std::vector<double>& values)

    std::cout << "Five" << std::endl;

   std::cout << "Five" << std::endl;
   std::cout << "Five" << std::endl;
   std::cout << "Five" << std::endl;
   std::cout << "XXXX" << std::endl;

   return -1;

//        inner_computator.ConstructCloudInvariantCombo(cloudB, res,1);
}






    void write_csv(const std::string& filename,
                   const std::vector<std::vector<int>>& combos,
                   const std::vector<arma::mat>& NCD1,
                   const std::vector<arma::mat>& NCD2,
                   const std::vector<double>& values)
    {

        std::ofstream file(filename);
        if (!file.is_open())
        {
            std::cerr << "Failed to open " << filename << std::endl;
            return;
        }

        // Write header
        file << "combo,distribution,distance_matrix,sign\n";

        for (size_t i = 0; i < combos.size(); ++i)
        {
            // Convert combo to string
            std::ostringstream combo_stream;
            for (size_t j = 0; j < combos[i].size(); ++j)
            {
                combo_stream << combos[i][j];
                if (j != combos[i].size() - 1)
                    combo_stream << " ";
            }

            // Convert arma::mat to string
            auto mat_to_string = [](const arma::mat& mat) -> std::string
            {
                std::ostringstream ss;
                ss << "[";
                for (size_t r = 0; r < mat.n_rows; ++r)
                {
                    ss << "[";
                    for (size_t c = 0; c < mat.n_cols; ++c)
                    {
                        ss << mat(r, c);
                        if (c != mat.n_cols - 1)
                            ss << ",";
                    }
                    ss << "]";
                    if (r != mat.n_rows - 1)
                        ss << ",";
                }
                ss << "]";
                return ss.str();
            };

            std::string ncd1_str = mat_to_string(NCD1[i]);
            std::string ncd2_str = mat_to_string(NCD2[i]);

            // Write line to file
            file << "\"" << combo_stream.str() << "\","
                 << "\"" << ncd1_str << "\","
                 << "\"" << ncd2_str << "\","
                 << values[i] << "\n";
        }

        file.close();
        std::cout << "Data written to " << filename << std::endl;
    }


    int NumElem(std::pair<arma::mat, std::vector<std::vector<bool>>> & cloudAunc)
    {
        return cloudAunc.first.n_cols;

    }
    int Dimension(std::pair<arma::mat, std::vector<std::vector<bool>>> & cloudAunc)
    {

        return cloudAunc.first.n_rows;

    }

    void initilize(std::pair<arma::mat, std::vector<std::vector<bool>>> & cloudAunc, std::pair<arma::mat, std::vector<std::vector<bool>>> & cloudA, int index)
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

        rads[index] = radius;
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



    void ConstructCloudInvariantCombo(std::pair<arma::mat, std::vector<std::vector<bool>>> & cloudA, std::vector<std::vector<int>> & combos, int index)
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
            double sign = 0;
            if (cloudA.first.n_rows == combos[i].size())
            {
                sign = ComputeSign(selectedPoints);
            }

            //!std::get<0>(finalCollection) = sign;
            arma::mat distMat;
            //  std::cout << "Matrix selection succesful" << std::endl;
            //  off << "Combo: " << i << std::endl;
            this->ProduceSCDMatrixM(selectedPoints,notSelectedpoints,ocdCloud, cloudA.second, combos[i], counterCombos, cloudA.first.n_cols);
            //  std::cout << "Produce OCD succeful" << std::endl;
            this->ProduceDistanceMat(selectedPoints,distMat, cloudA.second, combos[i], cloudA.first.n_cols);


            //  std::cout << "Produce distMat succesful" << std::endl;
            signs[index].push_back(sign);
            ocd[index].push_back(ocdCloud);
            dists[index].push_back(distMat);
            //!    auto t = std::make_tuple(sign, distMat, ocdCloud);
            //!    finalCollection.push_back(t);

            //!    OCDCloudCollection.push_back(ocdCloud);
            //!   distMatrixCollection.push_back(distMat);
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

        ocdCloud.resize(selectedPoints.n_cols+1,notselectedpoints.n_cols);
        arma::vec origin(selectedPoints.n_rows, arma::fill::zeros);

        for (int j = 0; j < notselectedpoints.n_cols; j++)
        {
            for (int i = 0 ; i < selectedPoints.n_cols; i++)
            {

                double dist = arma::norm(notselectedpoints.col(j)-selectedPoints.col(i));
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

    double ComputeSign(arma::mat & selectedPoints)
    {
        double dett = arma::det(selectedPoints);
        double answer = 0;

        if (dett < 0)
        {
            answer = -1;
        }
        if (dett > 0)
        {
            answer = +1;
        }

        arma::mat volumeObj(selectedPoints.n_rows, selectedPoints.n_cols + 1, arma::fill::zeros);

        for (int a = 0; a < selectedPoints.n_cols; a++)
        {
            volumeObj.col(a) = selectedPoints.col(a);
        }

        //  volumeObj.col(selectedPoints.n_cols) = notselectedpoints.col(j);

        double strengthh = ComputeStrength(volumeObj);



        return answer * strengthh;

    }
    virtual ~NCDGrabber() {}

protected:

private:
    double aConst;
    double constC;
//    bdcomputatorType  bd_computator;
    std::vector<double> rads;
    std::vector<std::vector<arma::mat>> ocd;
    std::vector<std::vector<arma::mat>> dists;
    std::vector<std::vector<std::vector<int>>> storedMatchings;
    std::vector<std::vector<double>> storedDists;
    std::vector<int> furthestPoint;
    std::vector<int> finalMatch;
    std::vector<std::vector<double>> signs;

};

#endif // NCDGRABBER_H
