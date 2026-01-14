#ifndef ORIGINAL_BASE_H
#define ORIGINAL_BASE_H

#include <boost/unordered_map.hpp>
#include <unordered_map>
#include <set>
#include <chrono>
#include "CloudIsometry/GeneralInvariantMethods.h"

template<typename invariantType, typename bdcomputatorType, typename obj>
class Original_Base
{



public:
    Original_Base():inner_computator(invariantType()),bd_computator(bdcomputatorType()) {}
    virtual ~Original_Base() {}

    void find_comb(std::vector<std::vector<int>>& res,std::vector<int> temp, int ind,int st, int end,int k)
    {
        if(ind == k)
        {
            res.push_back(temp);
            return;
        }
        for(int i = st; i<end && end-i+1 >= k-ind; i++)
        {
            temp[ind] = i;
            find_comb(res,temp,ind+1,i+1,end,k);
        }
    }

    void permuteMatrix(arma::mat & input, arma::mat & output, std::vector<int> & permutation)
    {


        output.resize(input.n_rows,input.n_cols);


//        for (int i = 0; i < permutation.size(); i++)
//        {
//            output.row(permutation[i]) = input.row(i);
//
//
//        }

        for (int i = 0; i < input.n_rows; i++)
        {
            for (int j = 0; j < input.n_cols; j++)
            {
                output(permutation[i],permutation[j]) = input(i,j);
            }

        }
    }

    void PermuteRowAndReorder(arma::mat & input, arma::mat & output, std::vector<int> & permutation, int sign)
    {

        //!arma::mat tmpStep;
        output.resize(input.n_rows,input.n_cols);
        for (int i = 0; i < permutation.size(); i++)
        {
            output.row(permutation[i]) = input.row(i);
        }
        for(int i = 0; i < input.n_cols; i++)
        {
            output.row(output.n_rows-1)(i) = sign * input.row(output.n_rows-1)(i);
        }



        // for (int i = 0; i < tmpStep.n_cols; i++)
        //  {
        //      tmpStep(tmpStep.n_rows-1, i) = sign * tmpStep(tmpStep.n_rows-1, i);
        //   }
        //  this->LexiographicOrdering(tmpStep, output);
    }


    double Compute_Isometry(obj & cloudAunc, obj & cloudBunc, std::vector<int> & matching,long & totalTimeForSMI, long & totalTimeForBD)
    {
        using Clock = std::chrono::steady_clock;
        using std::chrono::time_point;
        using std::chrono::duration_cast;
        using std::chrono::milliseconds;
        using namespace std::literals::chrono_literals;
//using std::this_thread::sleep_for;
        using namespace std;



        boost::unordered::unordered_map<std::vector<int>,int> umap;
        std::vector<std::vector<int>> res;
        obj cloudA;
        obj cloudB;
        inner_computator.initilize(cloudAunc, cloudA);
        inner_computator.initilize(cloudBunc, cloudB);
        int n = inner_computator.NumElem(cloudAunc);
        int k = inner_computator.Dimension(cloudAunc) - 1;
        //!std::vector<int> current;
        int recnum = 0;
        std::vector<int> temp(k,0);
        // CrystData::find_comb(res,temp,0,0,n,k);
        find_comb(res, temp,0,0,n,k);
        std::vector<arma::mat> distA;
        std::vector<arma::mat> distB;
        std::vector<arma::mat> ocdA;
        std::vector<arma::mat> ocdB;
        time_point<Clock> start = Clock::now();
        inner_computator.ConstructCloudInvariantCombo(cloudA, res, ocdA, distA);
        inner_computator.ConstructCloudInvariantCombo(cloudB, res, ocdB, distB);
        time_point<Clock> endd = Clock::now();
        milliseconds diff = duration_cast<milliseconds>(endd - start);
        totalTimeForSMI = diff.count();
        time_point<Clock> startTwo = Clock::now();

        std::vector<std::vector<double>> distMatrix(res.size(),std::vector<double>(res.size()));
        for (int i = 0; i < res.size(); i++)
        {
            for (int j = 0; j < res.size(); j++)
            {


                arma::mat normalizeddistA = distA[i];
                arma::mat normalizeddistB = distB[j];
                arma::mat normalizedocdA = ocdA[i];
                arma::mat normalizedocdB = ocdB[j];




                std::vector<int> selectB;

                for (int ij = 0; ij < k; ij++ )
                {
                    selectB.push_back(ij);
                }
                double dist = std::numeric_limits<double>::max();
                int perm = 0;
                do
                {
                    std::vector<int> permuteTmp = selectB;
                    int sign = GeneralInvariantMethods::permutationSign(permuteTmp);


                    arma::mat normalizeddistBperm;
                    arma::mat normalizeddocdBperm;

                    permuteTmp.push_back(selectB.size());



                    permuteMatrix(normalizeddistB, normalizeddistBperm, permuteTmp);



                    permuteTmp.push_back(selectB.size()+1);
                    PermuteRowAndReorder(normalizedocdB, normalizeddocdBperm, permuteTmp, sign);




                    double permute = bd_computator.computeBetweenInvariantDistancePair(normalizeddistA, normalizeddistBperm, normalizedocdA, normalizeddocdBperm);


                    dist = std::min(dist,permute);


                    perm = perm + 1;
                }
                while (next_permutation(selectB.begin(), selectB.end()));

                distMatrix[i][j] = dist;


            }
        }
        double dist = bd_computator.computeForGraph(distMatrix,matching);

        time_point<Clock> enddTwo = Clock::now();
        milliseconds diffTwo = duration_cast<milliseconds>(enddTwo - startTwo);
        totalTimeForBD = diffTwo.count();
        return dist;



    }

protected:

private:
    invariantType inner_computator;
    bdcomputatorType bd_computator;
};

#endif // ORIGINAL_BASE_H
