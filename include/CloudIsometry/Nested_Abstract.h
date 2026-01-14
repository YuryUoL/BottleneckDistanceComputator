#ifndef NESTED_ABSTRACT_H
#define NESTED_ABSTRACT_H

#include <boost/unordered_map.hpp>
#include <unordered_map>
#include <set>

template<typename invariantType, typename innerType, typename bdcomputatorType, typename obj>
class Nested_Abstract
{
public:
    Nested_Abstract() {}
    Nested_Abstract(invariantType & t, bdcomputatorType & d):inner_computator(t), bd_computator(d) {}

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


    double recursiveLayeredDistanceComputator(std::vector<innerType> & DataA, std::vector<innerType> & DataB,
            std::vector<int> & matching, std::vector<int> A, std::vector<int> B, int k, int n,
            boost::unordered::unordered_map<std::vector<int>,int> & umap)
    {
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
            //  std::cout << "Retrieving distance" << std::endl;
            std::vector<int> vecA;
            ////  vecA.insert(vecA.end(), A.begin(), A.end());
            //  std::vector<int> vecB;
            //  vecB.insert(vecB.end(), B.begin(), B.end());

            int iii = umap[A];
            int jjj = umap[B];

            //arma::mat firstMat = distA[iii];
            //  arma::mat secondMat= distB[jjj];

            //   double linftydist = 0;
            //  double bddist = 0;

            double permute = inner_computator.computeSpecific(DataA[iii], DataB[jjj], bd_computator);

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
                            double kak = recursiveLayeredDistanceComputator(DataA, DataB, matching, newA, newB, k, n,umap);
                            distMat[kx][ky] = kak;
                            ky = ky + 1;
                        }

                    }
                    ky = 0;
                    kx = kx + 1;

                }
            }

            matching = std::vector<int>(n-A.size());

            double dist = bd_computator.computeForGraph(distMat,matching);


            return dist;
        }



    }



    double Compute_Isometry(obj & cloudAunc, obj & cloudBunc, std::vector<int> & matching,long & totalTimeForSMI, long & totalTimeForBD, int k)
    {

        using Clock = std::chrono::steady_clock;
        using std::chrono::time_point;
        using std::chrono::duration_cast;
        using std::chrono::milliseconds;
        using namespace std::literals::chrono_literals;
//        using std::this_thread::sleep_for;
        using namespace std;

        boost::unordered::unordered_map<std::vector<int>,int> umap;
        std::vector<std::vector<int>> res;
        obj cloudA;
        obj cloudB;

        // std::cout << "No fragmentation x0 " << std::endl;


        inner_computator.initilize(cloudAunc, cloudA);
        inner_computator.initilize(cloudBunc, cloudB);

        //  std::cout << "No fragmentation x1 " << std::endl;


        int n = inner_computator.NumElem(cloudAunc);
      //  int k = inner_computator.Dimension(cloudAunc)-1;
        std::vector<int> current;
        int recnum = 0;
        resInitilize(res, n, k, current, umap, recnum);

        //  std::cout << "No fragmentation x2" << std::endl;



      //!  std::vector<arma::mat> distA;
      //!  std::vector<arma::mat> distB;
      //!  std::vector<arma::mat> ocdA;
      //!  std::vector<arma::mat> ocdB;
        time_point<Clock> start = Clock::now();
        std::vector<innerType> aData;
        std::vector<innerType> bData;
        inner_computator.ConstructCloudInvariantCombo(cloudA, res, aData);
        inner_computator.ConstructCloudInvariantCombo(cloudB, res, bData);

        //   std::cout << "Invariants computed " << std::endl;

        // std::cout << "No fragmentation x3" << std::endl;


        time_point<Clock> endd = Clock::now();
        milliseconds diff = duration_cast<milliseconds>(endd - start);
        totalTimeForSMI = diff.count();
        time_point<Clock> startTwo = Clock::now();
        std::vector <int> AA;
        std::vector <int> BB;



        //  std::cout << distA[0].n_cols << std::endl;
        // std::cout << distB[0].n_cols << std::endl;

        // std::cout << ocdA[0].n_cols << std::endl;

        //  std::cout << ocdB[0].n_cols << std::endl;



        //   std::cout << "Starting recursive distance computations " << std::endl;

        double dumdum = recursiveLayeredDistanceComputator(aData, bData, matching,AA, BB, k,n, umap);

        // std::cout << "Recursive distance computations have been finnished " << std::endl;
        // std::cout << "No fragmentation x4" << std::endl;

        time_point<Clock> enddTwo = Clock::now();
        milliseconds diffTwo = duration_cast<milliseconds>(enddTwo - startTwo);
        totalTimeForBD = diffTwo.count();

        return dumdum;


    }

    virtual ~Nested_Abstract() {}

protected:

private:
    invariantType  inner_computator;
    bdcomputatorType  bd_computator;
};

#endif // NESTED_ABSTRACT_H
