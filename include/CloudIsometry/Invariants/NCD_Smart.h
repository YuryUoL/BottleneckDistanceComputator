#ifndef NCD_SMART_H
#define NCD_SMART_H

template<typename bdcomputatorType>
class NCD_Smart
{
public:
    NCD_Smart()
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
    }
    virtual ~NCD_Smart() {}

    int NumElem(arma::mat & cloudAunc)
    {
        return cloudAunc.n_cols;

    }
    int Dimension(arma::mat & cloudAunc)
    {

        return cloudAunc.n_rows;

    }

    void initilize(arma::mat & cloudAunc, arma::mat & cloudA, int index)
    {

        arma::mat cloudOUT;
        double sizeofcloudin = cloudAunc.n_cols;

        arma::vec centroid(cloudAunc.n_rows);
        double radius = 0;
        for (int i = 0; i < cloudAunc.n_cols; i++)
        {
            arma::vec tmpZ =  cloudAunc.col(i);
            centroid = centroid + tmpZ;
        }
        double division = -1/sizeofcloudin;

        centroid = division * centroid;

        cloudOUT = cloudAunc;
        for (int i = 0; i < cloudOUT.n_cols; i++)
        {

            arma::vec blablaz = cloudAunc.col(i);
            arma::vec tmpd = blablaz + centroid;
            cloudOUT.col(i) = tmpd;
            double dudu = arma::norm(tmpd,2);

            if (dudu > radius)
            {
                radius = dudu;
                furthestPoint[index] = i;

            }
         //   radius = std::max(radius,dudu);

            //! radius = std::max(radius, sqrt(vectorTools::SquaredNorm(cloudOUT[i])));

        }

        cloudA = cloudOUT;
        rads[index] = radius;
        //    return radius;
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

    void ConstructCloudInvariantCombo(arma::mat & cloudA, std::vector<std::vector<int>> & combos,int index)
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
            this->selectMatrix(combos[i],cloudA,selectedPoints);
            this->theOthetMatrix(combos[i],counterCombos,cloudA,notSelectedpoints);
            arma::mat ocdCloud;
            arma::mat distMat;
            //  std::cout << "Matrix selection succesful" << std::endl;
            //  off << "Combo: " << i << std::endl;
            this->ProduceSCDMatrixM(selectedPoints,notSelectedpoints,ocdCloud, combos[i], counterCombos, cloudA.n_cols, index);
            //  std::cout << "Produce OCD succeful" << std::endl;
            this->ProduceDistanceMat(selectedPoints,distMat, combos[i], cloudA.n_cols, index);
            //  std::cout << "Produce distMat succesful" << std::endl;


            ocd[index].push_back(ocdCloud);
            dists[index].push_back(distMat);
        }

        storedMatchings = std::vector<std::vector<std::vector<int>>>(ocd[0].size(), std::vector<std::vector<int>>(ocd[1].size(), std::vector<int>()));
        storedDists = std::vector<std::vector<double>>(ocd[0].size(), std::vector<double>(ocd[1].size(), 0));

    }

    void ProduceSCDMatrixM(arma::mat & selectedPoints, arma::mat & notselectedpoints, arma::mat & ocdCloud, std::vector<int> & combos, std::vector<int> & counterCombos, int NumPoints, int index)
    {
        ocdCloud.resize(selectedPoints.n_cols+2,notselectedpoints.n_cols);
        arma::vec origin(selectedPoints.n_rows, arma::fill::zeros);

        for (int j = 0; j < notselectedpoints.n_cols; j++)
        {
            arma::mat finalDeterminant(notselectedpoints.n_rows, selectedPoints.n_rows);
            for (int i = 0 ; i < selectedPoints.n_cols; i++)
            {

                double dist = arma::norm(notselectedpoints.col(j)-selectedPoints.col(i));
                finalDeterminant.col(i) = notselectedpoints.col(j)-selectedPoints.col(i);

                ocdCloud(i,j) = (dist *dist) / (aConst * rads[index]);
            }

            double tmpdist = arma::norm(origin-notselectedpoints.col(j));

            ocdCloud(selectedPoints.n_cols,j) = (tmpdist * tmpdist) / (aConst * rads[index]);


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

    void ProduceDistanceMat(arma::mat & selectedPoints, arma::mat & distMat, std::vector<int> & combos, int numPoints, int index)
    {
        distMat.resize(selectedPoints.n_cols + 1, selectedPoints.n_cols + 1);
        arma::vec origin(selectedPoints.n_rows, arma::fill::zeros);
        for (int i  = 0; i < selectedPoints.n_cols; i++)
        {
            double dist1 = arma::norm(origin-selectedPoints.col(i),2);
            distMat(selectedPoints.n_cols,i) = dist1;
            distMat(i,selectedPoints.n_cols) = dist1;
            for (int j = 0; j < selectedPoints.n_cols; j++)
            {
                double dist = arma::norm(selectedPoints.col(i)-selectedPoints.col(j),2);

                distMat(j,i) = (dist * dist)/rads[index];
            }
        }
    }

    void clearAll()
    {
        ocd.clear();
        dists.clear();
        rads.clear();
        storedMatchings.clear();
            //!std::vector<std::vector<std::vector<int>>> storedMatchings;
        storedDists.clear();



    }

    double baseCase(int iii, int jjj)
    {

        double distocd;
        double permute = bd_computator.computeBetweenInvariantDistancePair(dists[0][iii], dists[1][jjj], ocd[0][iii], ocd[1][jjj], distocd, storedMatchings[iii][jjj]);
      //  std::cout << "For pair " << iii << " --- " << jjj << " Storing " << std::endl;
      //  std::cout <"pritning" << std::endl;
      //  printVector(storedMatchings[iii][jjj]);

        storedDists[iii][jjj] = distocd;
        double tmpdu = abs(rads[1] - rads[0]);
        return std::max(permute, tmpdu);

    }

    double score(std::vector<std::vector<double>> & distMat, std::vector<int> & A, std::vector<int> & B)
    {
        //std::vector<int> matcher;

        double dist = bd_computator.computeForGraph(distMat,finalMatch);
        return dist;

    }

        void printVector(std::vector<int> & in)
    {
        //  std::cout << "Size of vector: " << in.size() << std::endl;
        for (int i = 0; i < in.size(); i++)
        {

            std::cout << in[i] << ", ";
        }
        std::cout << std::endl;
    }


     void findMatching(std::vector<int> & matching, std::vector<std::vector<int>> & res, int n)
     {

        std::cout << "Radius A computed as: " << rads[0] << std::endl;
        std::cout << "Radius B computed as: "<< rads[1] << std::endl;


        int saveI = furthestPoint[0];
        int saveJ = finalMatch[saveI];



        std::cout << "First match is : " << saveI << " --- " << saveJ << std::endl;

        std::cout << "Between ocds: " << std::endl;

        std::cout << ocd[0][saveI] << std::endl;
        std::cout << ocd[1][saveJ] << std::endl;


        std::vector<int> TestMatching;
        bd_computator.computeBetweenCloudPair(ocd[0][saveI], ocd[1][saveJ]  , TestMatching);
        std::cout << "Testmatching: " << std::endl;
        printVector(TestMatching);


        std::set<int> sA(res[saveI].begin(), res[saveI].end());
        std::set<int> sB(res[saveJ].begin(), res[saveJ].end());

        std::vector<int> counterA;
        std::vector<int> counterB;

        for (int i = 0 ; i < n ; i++)
        {
            if (sA.find(i) == sA.end())
            {
                counterA.push_back(i);
            }
        }

         for (int i = 0 ; i < n ; i++)
        {
            if (sB.find(i) == sB.end())
            {
                counterB.push_back(i);
            }
        }

        std::vector<int> finalM(n, -1);
        finalM[res[saveI][0]] = res[saveJ][0];
        std::vector<int> corrMatch = storedMatchings[saveI][saveJ];

        std::cout << "Counter A : " << std::endl;
        printVector(counterA);
        std::cout << "Counter B : " << std::endl;
        printVector(counterB);
        std::cout << "Corrmatch: " << std::endl;
        printVector(corrMatch);


        for (int i = 0; i < corrMatch.size(); i++)
        {
     //   finalM[] = counterB[];
        finalM[counterA[i]] = counterB[corrMatch[i]];



        }




        matching = finalM;




     }

//    void findMatching(std::vector<int> & matching, std::vector<std::vector<int>> & res, int n)
//    {
//        int saveI;
//        int saveJ;
//        double minini = std::numeric_limits<double>::max();
//        for (int i = 0; i < storedDists.size(); i++)
//        {
//            for (int j = 0 ; j < storedDists[0].size(); j++)
//            {
//                if (storedDists[i][j] < minini)
//                {
//                    saveI = i;
//                    saveJ = j;
//                    minini = storedDists[i][j];
//
//                }
//            }
//        }
//
//
//        std::set<int> sA(res[saveI].begin(), res[saveI].end());
//        std::set<int> sB(res[saveJ].begin(), res[saveJ].end());
//
//        std::vector<int> counterA;
//        std::vector<int> counterB;
//
//        for (int i = 0 ; i < n ; i++)
//        {
//            if (sA.find(i) == sA.end())
//            {
//                counterA.push_back(i);
//            }
//        }
//
//         for (int i = 0 ; i < n ; i++)
//        {
//            if (sB.find(i) == sB.end())
//            {
//                counterB.push_back(i);
//            }
//        }
//
//        std::vector<int> finalM(n, 0);
//        finalM[res[saveI][0]] = res[saveJ][0];
//        for (int i = 0; i < storedMatchings[saveI][saveJ].size(); i++)
//        {
//        finalM[counterA[i]] = counterB[storedMatchings[saveI][saveJ][i]];
//
//
//
//        }
//
//
//
//
//        matching = finalM;
//
//    }



protected:

private:

    double aConst;
    double constC;
    bdcomputatorType  bd_computator;
    std::vector<double> rads;
    std::vector<std::vector<arma::mat>> ocd;
    std::vector<std::vector<arma::mat>> dists;
    std::vector<std::vector<std::vector<int>>> storedMatchings;
    std::vector<std::vector<double>> storedDists;
    std::vector<int> furthestPoint;
    std::vector<int> finalMatch;


};

#endif // NCD_SMART_H
