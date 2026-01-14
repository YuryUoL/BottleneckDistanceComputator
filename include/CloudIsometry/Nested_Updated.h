#ifndef NESTED_UPDATED_H
#define NESTED_UPDATED_H

template<typename invariantType, typename obj>
class Nested_Updated
{
public:
    Nested_Updated() {}
    //:debug(false) {}
    // Nested_Updated(std::string debugPath):dp(debugPath),debug(true) {}
    virtual ~Nested_Updated() {}

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

//    void printVector(std::vector<int> & in)
//    {
//        //  std::cout << "Size of vector: " << in.size() << std::endl;
//        for (int i = 0; i < in.size(); i++)
//        {
//
//            std::cout << in[i] << ", ";
//        }
//        std::cout << std::endl;
//    }


    void printVector(std::vector<int> & in)
    {
        //  std::cout << "Size of vector: " << in.size() << std::endl;
        for (int i = 0; i < in.size(); i++)
        {

            std::cout << in[i] << ", ";
        }
        std::cout << std::endl;
    }



    void printVector(std::vector<int> & in, std::ofstream & debugS)
    {
        //  std::cout << "Size of vector: " << in.size() << std::endl;
        for (int i = 0; i < in.size(); i++)
        {

            debugS<< in[i] << ", ";
        }
        debugS << std::endl;
    }

    template<typename tt>
    void printMatrix(std::vector<std::vector<tt>> & in, std::ofstream & debugS)
    {
        //  std::cout << "Size of vector: " << in.size() << std::endl;
        for (int i = 0; i < in.size(); i++)
        {
            for (int j = 0; j < in[0].size(); j++)
            {
                if (j < in[0].size() - 1)
                {
                    debugS<< in[i][j] << ", ";
                }
                else
                {
                    debugS<< in[i][j];
                }

            }
            debugS << std::endl;
        }

    }

    void printMatrix(arma::mat & in, std::ofstream & debugS)
    {
        debugS << in << std::endl;
        //  std::cout << "Size of vector: " << in.size() << std::endl;
//        for (int i = 0; i < in.size(); i++)
//        {
//            for (int j = 0; j < in[0].size(); j++)
//            {
//                debugS<< in[i][j] << ", ";
//
//            }
//            debugS << std::endl;
//        }

    }

    void printHandling(std::vector<int> & in, std::vector<int> & two, std::ofstream & debugS)
    {

        debugS << "Handling: ";

        for (int i = 0; i < in.size(); i++)
        {

            debugS<< in[i] << ", ";
        }
        debugS << "--- ";
        for (int i = 0; i < two.size(); i++)
        {

            debugS<< two[i] << ", ";
        }
        debugS << std::endl;

    }

    double recursiveLayeredDistanceComputator(std::vector<int> A, std::vector<int> B, int k, int n,
            boost::unordered::unordered_map<std::vector<int>,int> & umap, std::ofstream & debugS)
    {




        std::unordered_set<int> AduplicateCheck(A.begin(), A.end());
        std::unordered_set<int> BduplicateCheck(B.begin(), B.end());
        // std::cout << "Running " << current.size() << std::endl;
        if (A.size() == k)
        {
            //   std::cout << "BaseCase" << std::endl;
            //    std::cout << "Retrieving distance" << std::endl;
            //   std::vector<int> vecA;
            ////  vecA.insert(vecA.end(), A.begin(), A.end());
            //  std::vector<int> vecB;
            //  vecB.insert(vecB.end(), B.begin(), B.end());

            int iii = umap[A];
            int jjj = umap[B];



            //arma::mat firstMat = distA[iii];
            //  arma::mat secondMat= distB[jjj];

            //   double linftydist = 0;
            //  double bddist = 0;

            //   std::cout << "BaseCase activated" << iii << " |  " << jjj  << std::endl;
            //!  printVector(A, debugS);
            //!   printVector(B, debugS);
            printHandling(A, B, debugS);
            double permute = inner_computator.baseCase(iii,jjj,A,B, debugS);
            //   std::cout << "Out of basecase" << std::endl;
            //    std::cout << "Permute computed as : " << permute << std::endl;
            //      double permute = bd_computator.computeBetweenInvariantDistancePair(distA[iii], distB[jjj], ocdA[iii], ocdB[jjj]);
            //     outs << "Combo " << iii << ", " << jjj << " : " << permute << std::endl;
            return permute;
            // double linftydist = inner_computator.computeBetweenDistanceClouds(distA[iii], distB[jjj]);
            //  = inner_computator.computeBetweenInvariants(ocdA[iii], ocdB[jjj]);

            //    double permute = std::max(linftydist, bddist);



            //  return permute;

            // double damm = retrieveDistance(A,B, distMatrix, umap);
            // return damm;
        }
        else
        {
            //  std::cout << "Score" << std::endl;
            std::vector<std::vector<double>> distMat(n - A.size(), std::vector<double>(n - B.size()));
            int kx = 0;
            int ky = 0;
            for (int i = 0; i < n; i++)
            {
                if (AduplicateCheck.find(i) == AduplicateCheck.end())
                {
                    for (int j = 0; j < n; j++)
                    {
                        if (BduplicateCheck.find(j) == BduplicateCheck.end())
                        {
                            std::vector <int> newA = A;
                            std::vector <int> newB = B;
                            newA.push_back(i);
                            newB.push_back(j);
                            //! ,outs
                            //        std::cout << "Visiting " << kx << " , " << ky << std::endl;
                            double kak = recursiveLayeredDistanceComputator(newA, newB, k, n,umap, debugS);

                            //   std::cout << "DistMaT size: " << distMat.size() << std::endl;
                            //    std::cout << "Allocated : " << kx << " , " << ky << " dist: " << kak << std::endl;

                            distMat[kx][ky] = kak;
                            ky = ky + 1;
                        }

                    }
                    ky = 0;
                    kx = kx + 1;

                }
            }

//!            matching = std::vector<int>(n-A.size());

            //  std::cout << "Score computing" << std::endl;
            printHandling(A, B, debugS);
            printMatrix(distMat, debugS);

            double dist = inner_computator.score(distMat, debugS, A, B);
            //  std::cout << "Score computed" << std::endl;


            //      outs << "Distnace computation was : " << dist << std::endl;
            //  outs << "Larger combo: " << iii << " , " << jjj << " : " << dist << std::endl;




            return dist;
        }



    }

    double recursiveLayeredDistanceComputator(std::vector<int> A, std::vector<int> B, int k, int n,
            boost::unordered::unordered_map<std::vector<int>,int> & umap)
    {


        //  std::cout << "Comparing" << A[0] << " to " << B[0] << std::endl;
        // std::cout << "Vector A: " << std::endl;
        //  std::cout << "Printing vectors A,B: " << std::endl;
        //   printVector(A);
        //     printVector(B);
        //    printVector(A, outs);
        //   outs << "Vector B: " << std::endl;
        //   printVector(B,  outs);
        //std::cout << "Enter Correct Method [x]" << std::endl;

//        std::cout << "Set A is: ";
//        for (int i = 0; i < A.size(); i++)
//        {
//            std::cout << A[i] << ", ";
//        }
//        std::cout << std::endl;

        std::unordered_set<int> AduplicateCheck(A.begin(), A.end());
        std::unordered_set<int> BduplicateCheck(B.begin(), B.end());
        // std::cout << "Running " << current.size() << std::endl;
        if (A.size() == k)
        {
            //   std::cout << "BaseCase" << std::endl;
            //    std::cout << "Retrieving distance" << std::endl;
            //   std::vector<int> vecA;
            ////  vecA.insert(vecA.end(), A.begin(), A.end());
            //  std::vector<int> vecB;
            //  vecB.insert(vecB.end(), B.begin(), B.end());

            int iii = umap[A];
            int jjj = umap[B];



            //arma::mat firstMat = distA[iii];
            //  arma::mat secondMat= distB[jjj];

            //   double linftydist = 0;
            //  double bddist = 0;

            //   std::cout << "BaseCase activated" << iii << " |  " << jjj  << std::endl;
            double permute = inner_computator.baseCase(iii,jjj,A,B);
            //   std::cout << "Out of basecase" << std::endl;
            //    std::cout << "Permute computed as : " << permute << std::endl;
            //      double permute = bd_computator.computeBetweenInvariantDistancePair(distA[iii], distB[jjj], ocdA[iii], ocdB[jjj]);
            //     outs << "Combo " << iii << ", " << jjj << " : " << permute << std::endl;
            return permute;
            // double linftydist = inner_computator.computeBetweenDistanceClouds(distA[iii], distB[jjj]);
            //  = inner_computator.computeBetweenInvariants(ocdA[iii], ocdB[jjj]);

            //    double permute = std::max(linftydist, bddist);



            //  return permute;

            // double damm = retrieveDistance(A,B, distMatrix, umap);
            // return damm;
        }
        else
        {
            //  std::cout << "Score" << std::endl;
            std::vector<std::vector<double>> distMat(n - A.size(), std::vector<double>(n - B.size()));
            int kx = 0;
            int ky = 0;
            for (int i = 0; i < n; i++)
            {
                if (AduplicateCheck.find(i) == AduplicateCheck.end())
                {
                    for (int j = 0; j < n; j++)
                    {
                        if (BduplicateCheck.find(j) == BduplicateCheck.end())
                        {
                            std::vector <int> newA = A;
                            std::vector <int> newB = B;
                            newA.push_back(i);
                            newB.push_back(j);
                            //! ,outs
                            //        std::cout << "Visiting " << kx << " , " << ky << std::endl;
                            double kak = recursiveLayeredDistanceComputator(newA, newB, k, n,umap);

                            //   std::cout << "DistMaT size: " << distMat.size() << std::endl;
                            //    std::cout << "Allocated : " << kx << " , " << ky << " dist: " << kak << std::endl;

                            distMat[kx][ky] = kak;
                            ky = ky + 1;
                        }

                    }
                    ky = 0;
                    kx = kx + 1;

                }
            }

//!            matching = std::vector<int>(n-A.size());

            //  std::cout << "Score computing" << std::endl;
            double dist = inner_computator.score(distMat, A, B);
            //  std::cout << "Score computed" << std::endl;


            //      outs << "Distnace computation was : " << dist << std::endl;
            //  outs << "Larger combo: " << iii << " , " << jjj << " : " << dist << std::endl;




            return dist;
        }
    }

    double Compute_Isometry(obj & cloudAunc, obj & cloudBunc, std::vector<int> & matching,long & totalTimeForSMI, long & totalTimeForBD)
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
        obj cloudA;
        obj cloudB;
        invariantType tmp;
        inner_computator = tmp;
        int k = inner_computator.Dimension(cloudAunc)-1;

        return Compute_Isometry(cloudAunc,cloudBunc, matching,totalTimeForSMI,totalTimeForBD, k);

    }

    double Compute_Isometry(obj & cloudAunc, obj & cloudBunc, std::vector<int> & matching,long & totalTimeForSMI, long & totalTimeForBD, int k)
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
        obj cloudA;
        obj cloudB;
     //   invariantType tmp;
     //   inner_computator = tmp;

//        std::cout << "Startg 123" << std::endl;

        //std::cout << "One " << std::endl;

        inner_computator.initilize(cloudAunc, cloudA,0);

      //   std::cout << "Two " << std::endl;

        inner_computator.initilize(cloudBunc, cloudB,1);

       //  std::cout << "Three " << std::endl;

        int n = inner_computator.NumElem(cloudAunc);
        std::vector<int> current;
        int recnum = 0;
      //  std::cout << "Resinitilize : " << n << " , " << k << std::endl;
        resInitilize(res, n, k, current, umap, recnum);
        time_point<Clock> start = Clock::now();

     //   std::cout << "Four " << std::endl;

        inner_computator.ConstructCloudInvariantCombo(cloudA, res,0);

      //   std::cout << "Five" << std::endl;

        inner_computator.ConstructCloudInvariantCombo(cloudB, res,1);

     //    std::cout << "Six " << std::endl;

        time_point<Clock> endd = Clock::now();
        //! milliseconds diff = duration_cast<milliseconds>(endd - start);
        microseconds diff = duration_cast<microseconds>(endd - start);
        totalTimeForSMI = diff.count();
        time_point<Clock> startTwo = Clock::now();
        std::vector <int> AA;
        std::vector <int> BB;
        double dumdum = recursiveLayeredDistanceComputator(AA, BB, k,n, umap);

        // std::cout << "Seven " << std::endl;

        inner_computator.findMatching(matching, res, n);

        // std::cout << "Eight " << std::endl;

        time_point<Clock> enddTwo = Clock::now();
        //!milliseconds diffTwo = duration_cast<milliseconds>(enddTwo - startTwo);
        microseconds diffTwo = duration_cast<microseconds>(enddTwo - startTwo);
        totalTimeForBD = diffTwo.count();

    //    std::cout << "Preparing exit " << std::endl;

        return dumdum;


    }

    double Compute_Isometry(obj & cloudAunc, obj & cloudBunc, std::vector<int> & matching,long & totalTimeForSMI, long & totalTimeForBD, std::ofstream & debugS)
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
        obj cloudA;
        obj cloudB;
        invariantType tmp;
        inner_computator = tmp;
        //inner_computator = invariantType();
        int k = inner_computator.Dimension(cloudAunc)-1;

        return Compute_Isometry(cloudAunc,cloudBunc, matching,totalTimeForSMI,totalTimeForBD, k, debugS);

    }



    double Compute_Isometry(obj & cloudAunc, obj & cloudBunc, std::vector<int> & matching,long & totalTimeForSMI, long & totalTimeForBD, int k, std::ofstream & debugS)
    {

        using Clock = std::chrono::steady_clock;
        using std::chrono::time_point;
        using std::chrono::duration_cast;
        using std::chrono::milliseconds;
        using namespace std::literals::chrono_literals;
        using namespace std;
        boost::unordered::unordered_map<std::vector<int>,int> umap;
        std::vector<std::vector<int>> res;
        obj cloudA;
        obj cloudB;
        invariantType tmp;
        inner_computator = tmp;
        //   void printMatrix(std::vector<std::vector<tt>> & in, std::ofstream & debugS)

        debugS << "Cloud A:" << std::endl;
        //  debugS << cloudAunc << std::endl;
        printMatrix( cloudAunc, debugS);

        debugS << "Cloud B:" << std::endl;
        // debugS << cloudBunc << std::endl;
        printMatrix( cloudBunc, debugS);


        inner_computator.initilize(cloudAunc, cloudA,0,debugS);

        inner_computator.initilize(cloudBunc, cloudB,1,debugS);
        int n = inner_computator.NumElem(cloudAunc);
        //int k = inner_computator.Dimension(cloudAunc)-1;
        std::vector<int> current;
        int recnum = 0;
        resInitilize(res, n, k, current, umap, recnum);

        time_point<Clock> start = Clock::now();
        debugS << "Initilizing A:" << std::endl;
        inner_computator.ConstructCloudInvariantCombo(cloudA, res,0, debugS);
        debugS << "Initilizing B" << std::endl;
        inner_computator.ConstructCloudInvariantCombo(cloudB, res,1, debugS);


        time_point<Clock> endd = Clock::now();
        milliseconds diff = duration_cast<milliseconds>(endd - start);
        totalTimeForSMI = diff.count();
        time_point<Clock> startTwo = Clock::now();
        std::vector <int> AA;
        std::vector <int> BB;
        double dumdum = recursiveLayeredDistanceComputator(AA, BB, k,n, umap, debugS);
        inner_computator.findMatching(matching, res, n);
        time_point<Clock> enddTwo = Clock::now();
        milliseconds diffTwo = duration_cast<milliseconds>(enddTwo - startTwo);
        totalTimeForBD = diffTwo.count();
        debugS << " Super final answer is : " << dumdum << std::endl;
        return dumdum;


    }


protected:

private:
    invariantType  inner_computator;
    std::string dp;
    bool debug;
    //  bdcomputatorType  bd_computator;
};

#endif // NESTED_UPDATED_H
