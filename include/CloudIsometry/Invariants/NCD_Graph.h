#ifndef NCD_GRAPH_H
#define NCD_GRAPH_H

template<typename bdcomputatorType>
class NCD_Graph
{
public:
    NCD_Graph()
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

    int NumElem(std::vector<std::vector<double>> & cloudAunc)
    {
        return cloudAunc.size();

    }


    void initilize(std::vector<std::vector<double>> & cloudAunc, std::vector<std::vector<double>> & cloudA, int index)
    {
        std::vector<std::vector<double>> k = cloudAunc;
        cloudA = k;

        //! Do nothing
    }

    void initilize(std::vector<std::vector<double>> & cloudAunc, std::vector<std::vector<double>> & cloudA, int index, std::ofstream & tp)
    {
        std::vector<std::vector<double>> k = cloudAunc;
        cloudA = k;

        //! Do nothing
    }

    void findCounterCombos(std::vector<int> & indices, std::vector<int> & counterCombos, int n)
    {
        std::set<int> s(indices.begin(), indices.end());
        for (int i = 0 ; i < n ; i++)
        {
            if (s.find(i) == s.end())
            {
                counterCombos.push_back(i);
                //  output.col(j) = input.col(i);
                //  j++;
            }
        }


    }

    void ProduceDistanceMat(std::vector<std::vector<double>> & dists, arma::mat & distMat, std::vector<int> & combos, int numPoints, int index)
    {
      //  std::cout << "Entered produce distance mat method" << std::endl;
     //   std::cout << " Dists size : " << dists.size() << std::endl;
     //   std::cout << "Dists: (" << dists.size() << " , " << dists[0].size() << ")" << std::endl;
        distMat.resize(combos.size(), combos.size() );
        //  arma::vec origin(selectedPoints.n_rows, arma::fill::zeros);
        for (int i  = 0; i < combos.size(); i++)
        {
            for (int j = 0; j < combos.size(); j++)
            {
             //   std::cout << "i: " << combos[i] << std::endl;
             //   std::cout << "j: " << combos[j] << std::endl;

                double duduz = dists[combos[i]][combos[j]];
                distMat(j,i) = duduz;
              // std::cout << "Succesful: " << duduz << std::endl;
                //! (dist * dist)/rads[index];
            }
        }
    //    std::cout << "Finnished succefully" << std::endl;
    }
    template<typename T>
      void printVector(std::vector<T> & in, std::ofstream & debugS)
    {
        //  std::cout << "Size of vector: " << in.size() << std::endl;
        for (int i = 0; i < in.size(); i++)
        {

            debugS << in[i] << ", ";
        }
        debugS << std::endl;
    }

    void ConstructCloudInvariantCombo(std::vector<std::vector<double>> & cloudA, std::vector<std::vector<int>> & combos,int index, std::ofstream & debugS)
    {
        for (int i = 0; i < combos.size(); i++)
        {
           // std::cout << "Staritng combo " << i << std::endl;

            //   arma::mat selectedPoints;
            //   arma::mat notSelectedpoints;

            //CrystData::ConstructSRCloud(cloud,basisCollection,combos[i],SrCloud);
            // CrystData::ConstructDistanceMatrix(basisCollection,DistMat);
            std::vector<int> counterCombos;
            //   std::cout << "Selecting matrices" << std::endl;
            this->findCounterCombos(combos[i], counterCombos, cloudA.size());
            //this->theOthetMatrix(combos[i],counterCombos,cloudA,notSelectedpoints);
            arma::mat ocdCloud;
            arma::mat distMat;
            //  std::cout << "Matrix selection succesful" << std::endl;
          //  std::cout << "Combo: " << i << std::endl;
            this->ProduceSCDMatrixM(cloudA, ocdCloud, combos[i], counterCombos, cloudA.size(), index);
          //    std::cout << "Produce OCD succeful" << std::endl;
            this->ProduceDistanceMat(cloudA, distMat, combos[i], cloudA.size(), index);
          //   std::cout << "Produce distMat succesful" << std::endl;


            ocd[index].push_back(ocdCloud);
            dists[index].push_back(distMat);
            printVector(combos[i] , debugS);
            debugS << ocdCloud << std::endl;
            debugS << distMat << std::endl;


            //debugS << << std::endl;

           // std::cout <<"Insert succesful" << std::endl;
        }

        storedMatchings = std::vector<std::vector<std::vector<int>>>(ocd[0].size(), std::vector<std::vector<int>>(ocd[1].size(), std::vector<int>()));
        storedDists = std::vector<std::vector<double>>(ocd[0].size(), std::vector<double>(ocd[1].size(), 0));

    }

    void ConstructCloudInvariantCombo(std::vector<std::vector<double>> & cloudA, std::vector<std::vector<int>> & combos,int index)
    {
        for (int i = 0; i < combos.size(); i++)
        {
           // std::cout << "Staritng combo " << i << std::endl;

            //   arma::mat selectedPoints;
            //   arma::mat notSelectedpoints;

            //CrystData::ConstructSRCloud(cloud,basisCollection,combos[i],SrCloud);
            // CrystData::ConstructDistanceMatrix(basisCollection,DistMat);
            std::vector<int> counterCombos;
            //   std::cout << "Selecting matrices" << std::endl;
            this->findCounterCombos(combos[i], counterCombos, cloudA.size());
            //this->theOthetMatrix(combos[i],counterCombos,cloudA,notSelectedpoints);
            arma::mat ocdCloud;
            arma::mat distMat;
            //  std::cout << "Matrix selection succesful" << std::endl;
          //  std::cout << "Combo: " << i << std::endl;
            this->ProduceSCDMatrixM(cloudA, ocdCloud, combos[i], counterCombos, cloudA.size(), index);
          //    std::cout << "Produce OCD succeful" << std::endl;
            this->ProduceDistanceMat(cloudA, distMat, combos[i], cloudA.size(), index);
          //   std::cout << "Produce distMat succesful" << std::endl;


            ocd[index].push_back(ocdCloud);
            dists[index].push_back(distMat);

           // std::cout <<"Insert succesful" << std::endl;
        }

        storedMatchings = std::vector<std::vector<std::vector<int>>>(ocd[0].size(), std::vector<std::vector<int>>(ocd[1].size(), std::vector<int>()));
        storedDists = std::vector<std::vector<double>>(ocd[0].size(), std::vector<double>(ocd[1].size(), 0));

    }

    void ProduceSCDMatrixM(std::vector<std::vector<double>> & dists, arma::mat & ocdCloud, std::vector<int> & combos, std::vector<int> & counterCombos, int NumPoints, int index)
    {

        ocdCloud.resize(combos.size(),counterCombos.size());
        for (int j  = 0; j < counterCombos.size(); j++)
        {
            for (int i = 0; i < combos.size(); i++)
            {
                ocdCloud(i,j) = dists[combos[i]][counterCombos[j]];
                //! (dist * dist)/rads[index];
            }
        }

    }

     void clearAll()
    {
        ocd.clear();
        dists.clear();
        rads.clear();
     //   storedMatchings.clear();
        //!std::vector<std::vector<std::vector<int>>> storedMatchings;
      //  storedDists.clear();



    }

    double baseCase(int iii, int jjj,std::vector<int> & A, std::vector<int> & B)
    {

        double distocd;
        std::vector<int> matchd;
        //! storedMatchings[iii][jjj]

        bdcomputatorType bd_computator_tmp;

    //    std::cout << "Dists[0] : " << dists[0].size() << std::endl;
     //   std::cout << "Ocd[0]: "<< ocd[0].size() << std::endl;
      //  std::cout << "Dists[1] : " << dists[1].size() << std::endl;
      //  std::cout << "Ocd[1]: "<< ocd[1].size() << std::endl;

        double permute = bd_computator_tmp.computeBetweenInvariantDistancePair(dists[0][iii], dists[1][jjj], ocd[0][iii], ocd[1][jjj], distocd, matchd);
        //  std::cout << "For pair " << iii << " --- " << jjj << " Storing " << std::endl;
        //  std::cout <"pritning" << std::endl;
        //  printVector(storedMatchings[iii][jjj]);

       // storedDists[iii][jjj] = distocd;
        double tmpdu = abs(rads[1] - rads[0]);
        return std::max(permute, tmpdu);

    }

      double baseCase(int iii, int jjj,std::vector<int> & A, std::vector<int> & B, std::ofstream & debugS)
    {

        double distocd;
        std::vector<int> matchd;
        //! storedMatchings[iii][jjj]

        bdcomputatorType bd_computator_tmp;

    //    std::cout << "Dists[0] : " << dists[0].size() << std::endl;
     //   std::cout << "Ocd[0]: "<< ocd[0].size() << std::endl;
      //  std::cout << "Dists[1] : " << dists[1].size() << std::endl;
      //  std::cout << "Ocd[1]: "<< ocd[1].size() << std::endl;

        double permute = bd_computator_tmp.computeBetweenInvariantDistancePair(dists[0][iii], dists[1][jjj], ocd[0][iii], ocd[1][jjj], distocd, matchd);
        //  std::cout << "For pair " << iii << " --- " << jjj << " Storing " << std::endl;
        //  std::cout <"pritning" << std::endl;
        //  printVector(storedMatchings[iii][jjj]);

       // storedDists[iii][jjj] = distocd;
        double tmpdu = abs(rads[1] - rads[0]);
        double ff = std::max(permute, tmpdu);

        debugS << "Basecase " << iii << ", " << jjj << " : " << ff << std::endl;

        return ff;

    }

    double score(std::vector<std::vector<double>> & distMat, std::ofstream & debugS, std::vector<int> & A, std::vector<int> & B)
    {

        double dist = bd_computator.computeForGraph(distMat,finalMatch);
        debugS << "Score: " << dist << std::endl;
        return dist;

    }

    double score(std::vector<std::vector<double>> & distMat , std::vector<int> & A, std::vector<int> & B)
    {
        //std::vector<int> matcher;

        double dist = bd_computator.computeForGraph(distMat,finalMatch);

  //      debugS << "Score: " << dist << std::endl;

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
        matching = this->finalMatch;
    }


//    void findMatching(std::vector<int> & matching, std::vector<std::vector<int>> & res, int n)
//    {
//
//     //   std::cout << "Radius A computed as: " << rads[0] << std::endl;
//       // std::cout << "Radius B computed as: "<< rads[1] << std::endl;
//
//
//        int saveI = furthestPoint[0];
//        int saveJ = finalMatch[saveI];
//
//
//
////        std::cout << "First match is : " << saveI << " --- " << saveJ << std::endl;
//
//  //      std::cout << "Between ocds: " << std::endl;
//
//    //    std::cout << ocd[0][saveI] << std::endl;
//      //  std::cout << ocd[1][saveJ] << std::endl;
//
//
//        std::vector<int>corrMatch;
//        bd_computator.computeBetweenCloudPair(ocd[0][saveI], ocd[1][saveJ], corrMatch);
//       // std::cout << "Testmatching: " << std::endl;
//       // printVector(TestMatching);
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
//        for (int i = 0 ; i < n ; i++)
//        {
//            if (sB.find(i) == sB.end())
//            {
//                counterB.push_back(i);
//            }
//        }
//
//        std::vector<int> finalM(n, 0);
//        finalM[res[saveI][0]] = res[saveJ][0];
//        //!std::vector<int> corrMatch = storedMatchings[saveI][saveJ];
//
//     //   std::cout << "Counter A : " << std::endl;
//    //    printVector(counterA);
//   //     std::cout << "Counter B : " << std::endl;
//    //    printVector(counterB);
//    //    std::cout << "Corrmatch: " << std::endl;
//   //     printVector(corrMatch);
//
//
//        for (int i = 0; i < corrMatch.size(); i++)
//        {
//            //   finalM[] = counterB[];
//            finalM[counterA[i]] = counterB[corrMatch[i]];
//
//        }
//
//
//
//
//        matching = finalM;
//
//
//
//
//    }
//
//    void findMatching(std::vector<int> & matching, std::vector<std::vector<int>> & res, int n)
//    {
//
//
//      //!  int saveI = furthestPoint[0];
//        int saveI = 0;
//        int saveJ = finalMatch[saveI];
//
//        selectAindx = saveI;
//        selectBindx = saveJ;
//
//        std::vector<int>corrMatch;
//        bd_computator.computeBetweenCloudPair(ocd[0][selectAindx], ocd[1][selectBindx], corrMatch);
//        // std::cout << "Testmatching: " << std::endl;
//        // printVector(TestMatching);
//
//
//        std::set<int> sA(res[selectAindx].begin(), res[selectAindx].end());
//        std::set<int> sB(res[selectBindx].begin(), res[selectBindx].end());
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
//        for (int i = 0 ; i < n ; i++)
//        {
//            if (sB.find(i) == sB.end())
//            {
//                counterB.push_back(i);
//            }
//        }
//
//        std::vector<int> finalM(n, 0);
//        for (int i = 0; i < res[selectAindx].size(); i++)
//        {
//            finalM[res[selectAindx][i]] = res[selectBindx][i];
//        }
//        //!std::vector<int> corrMatch = storedMatchings[saveI][saveJ];
//
//        //   std::cout << "Counter A : " << std::endl;
//        //    printVector(counterA);
//        //     std::cout << "Counter B : " << std::endl;
//        //    printVector(counterB);
//        //    std::cout << "Corrmatch: " << std::endl;
//        //     printVector(corrMatch);
//
//
//        for (int i = 0; i < corrMatch.size(); i++)
//        {
//            //   finalM[] = counterB[];
//            finalM[counterA[i]] = counterB[corrMatch[i]];
//
//        }
//
//
//
//
//        matching = finalM;
//
//
//        //!  matching = finalMatch;
//
//    }



    virtual ~NCD_Graph() {}

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
    std::vector<double> signs;
};

#endif // NCD_GRAPH_H
