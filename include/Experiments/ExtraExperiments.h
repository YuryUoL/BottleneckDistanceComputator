#ifndef EXTRAEXPERIMENTS_H
#define EXTRAEXPERIMENTS_H

#include "CloudIsometry/Invariants/vectorTools.h"
#include "CloudIsometry/Invariants/SCD_Updated.h"
#include "CloudIsometry/Nested_Abstract.h"
#include <math.h>

#include "CloudIsometry/Invariants/NCD_Smart.h"
#include "CloudIsometry/Invariants/NCD_Smart_Update.h"
#include "CloudIsometry/Invariants/PBI_Update.h"
#include "CloudIsometry/Nested_Updated.h"
#include "Bottleneckdistance/Misc/GraphTreeControllerVector.h"


class ExtraExperiments
{
public:
    ExtraExperiments() {}
    virtual ~ExtraExperiments() {}

    static void randomGraphGenerator(std::vector<std::vector<bool>> & matt)
    {
        for (int i = 0 ; i < matt.size(); i++)
        {
            for (int j = 0; j < matt[i].size(); j++)
            {
                int a=rand()%2;
                if (a > 0)
                {
                    matt[i][j] = 1;

                }
                else
                {
                    matt[i][j] = 0;

                }

            }

        }



    }

    static double fRand(double fMin, double fMax)
    {
        double f = (double)rand() / RAND_MAX;
        return fMin + f * (fMax - fMin);
    }


    static void generateRandomRotationMatrix2D(arma::mat & rotMat)
    {
        double phi = ExtraExperiments::fRand(0,2*M_PI);
        rotMat = arma::mat(2,2);
        rotMat(0,0) = cos(phi);
        rotMat(1,0) = -sin(phi);
        rotMat(0,1) = sin(phi);
        rotMat(1,1) = cos(phi);


    }

    static void generateRandomRotationMatrix3D(arma::mat & fMat)
    {

        fMat = arma::mat(3,3);
        double phi = ExtraExperiments::fRand(0,2*M_PI);
        double a = ExtraExperiments::fRand(0,2*M_PI);
        double b = ExtraExperiments::fRand(0,2*M_PI);
        arma::mat rotMat(3,3);
        rotMat(0,0) = cos(phi);
        rotMat(1,0) = -sin(phi);
        rotMat(0,1) = sin(phi);
        rotMat(1,1) = cos(phi);
        rotMat(2,2) = 1;

        arma::mat rotMatd(3,3);

        rotMatd(0,0) = cos(b);
        rotMatd(0,2) = sin(b);
        rotMatd(1,1) = 1;
        rotMatd(2,0) = -1*sin(b);
        rotMatd(2,2) = cos(b);

        arma::mat rotMatz(3,3);
        rotMatz(0,0) = 1;
        rotMatz(1,1) = cos(a);
        rotMatz(1,2) = -1*sin(a);
        rotMatz(2,1) = sin(a);
        rotMatz(2,2) = cos(a);

        fMat = rotMat;






    }

    static void generateRandomVector(arma::vec & p, int dim, double epp)
    {
        p = arma::vec(dim);
        //   std::cout << p << std::endl;
        //   std::cout << p.size() << std::endl;
        for (int i = 0; i < dim; i++)
        {
            //   std::cout << "In iteration i " << i << std::endl;
            double tmp = ExtraExperiments::fRand(-epp, epp);
            p[i] = tmp;
        }

        //   std::cout << "Method finnished succesfully " << std::endl;



    }

//    static void RotateCloud(arma::mat & rotatable, arma::mat & rotMat)
//    {
//
//
//    }

    static void RearrangeCloud(arma::mat & in, arma::mat & out, std::vector<int> & mapp)
    {
        out = arma::mat(in.n_rows, in.n_cols);
        for (int i = 0; i < mapp.size(); i++)
        {
            out.col(mapp[i]) = in.col(i);
        }


    }

    static void printVector(std::vector<int> & in)
    {
        //  std::cout << "Size of vector: " << in.size() << std::endl;
        for (int i = 0; i < in.size(); i++)
        {

            std::cout << in[i] << ", ";
        }
        std::cout << std::endl;
    }

    static void printVector(std::vector<double> & in)
    {
        //  std::cout << "Size of vector: " << in.size() << std::endl;
        for (int i = 0; i < in.size(); i++)
        {

            std::cout << in[i] << ", ";
        }
        std::cout << std::endl;
    }

    static void printMatrix( std::vector<std::vector<double>> & distMat)
    {

        for (int i = 0; i < distMat.size(); i++)
        {
            ExtraExperiments::printVector(distMat[i]);

        }

    }

    static void computeDistanceMatrix(arma::mat & A, arma::mat & B, std::vector<std::vector<double>> & distMat)
    {
        distMat = std::vector<std::vector<double>>(A.n_cols, std::vector<double>(B.n_cols));
        for (int i = 0; i < A.n_cols; i++)
        {
            for (int j = 0; j < B.n_cols; j++)
            {
                distMat[i][j] = arma::norm(A.col(i)-B.col(j), 2);

            }

        }


    }

    static void SampleTest(std::string Apath, std::string Bpath)
    {

        typedef BottleneckCalculator<Bipartie_Graph,GraphTreeController> bipartiebd;
        typedef bd_general_interface<bipartiebd, bipartiebd> naiveBdCalc;
        typedef bd_general_interface<BoostBottleneck,BoostBottleneck> boostBdCalc;
        arma::mat A;
        arma::mat B;
        A.load(Apath, arma::arma_ascii );
        B.load(Bpath, arma::arma_ascii );

        Nested_Base<SCD,  boostBdCalc, std::pair<arma::mat, std::vector<std::vector<bool>>>> boostSCD;

        std::vector<std::vector<bool>> graphA(A.n_cols+1, std::vector<bool>(A.n_cols+1, true));
        std::vector<std::vector<bool>> graphB(B.n_cols+1, std::vector<bool>(B.n_cols+1, true));
        long totalTimeForSMIstd;
        long totalTimeForBD;
        std::vector<int> matching;

        auto firstP = std::make_pair(A, graphA);
        auto secondP = std::make_pair(B, graphB);



        double res = boostSCD.Compute_Isometry(firstP,secondP, matching, totalTimeForSMIstd, totalTimeForBD);


        std::cout << "Cloud A:" << std::endl;
        std::cout << A << std::endl;
        std::cout << "Cloud B: " << std::endl;
        std::cout << B << std::endl;

        std::cout << " Result was : " << res << std::endl;
        std::cout << "Matching was: " << std::endl;
        ExtraExperiments::printVector(matching);


    }


    static void PerturbSingleNCDTest(int dim, int numPoints,  double curEpsion, int repeat, int numEpsilon, double & res, double & ratiod, double & realDist, double & ratiodistbd,
                                     double & totalTimeForSMIstdA, double & totalTimeForBDA, double boxS)
    {
        long totalTimeForSMIstd;
        long totalTimeForBD;
        int iter = numPoints;
        std::cout << "Iteration: " << iter << "curEpsilon: " << numEpsilon << " repeat: " << repeat << " started" << std::endl;
        typedef BottleneckCalculator<Bipartie_Graph,GraphTreeControllerVector> bipartiebd;
        typedef bd_general_interface<BoostBottleneck,BoostBottleneck> boostBdCalc;
        typedef bd_general_interface<bipartiebd, bipartiebd> naiveBdCalc;
        Nested_Updated<PBI_Update<naiveBdCalc>, arma::mat> nestedSCD;

        std::vector<std::vector<double>> baseCloud;
        ExtraExperiments::generateRandomCloudInBoxNd(baseCloud, numPoints, boxS, dim);

        std::vector<std::string> simpleStringVector;
        std::vector<std::vector<double>> perturbedCloud;
        double perturbError = ExtraExperiments::perturbCloud(baseCloud,perturbedCloud, curEpsion);

        arma::mat A;
        arma::mat B;

        vectorTools::vec_to_mat(baseCloud, A);
        vectorTools::vec_to_mat(perturbedCloud, B);

        std::cout << "Original A" << std::endl;
        std::cout << A << std::endl;
        std::cout << "Original B" << std::endl;
        std::cout << B << std::endl;



        std::vector<int> justPermute;
        for (int i = 0; i < numPoints; i++)
        {
            justPermute.push_back(i);

        }
        auto rng = std::default_random_engine {};

        //!   std::shuffle(std::begin(justPermute), std::end(justPermute), rng);


        arma::mat BB;
        ExtraExperiments::RearrangeCloud(B, BB, justPermute);

        EuclideanBDCompute<bipartiebd> eucBd;

        std::vector<int> matchingOne;
        double purebd = eucBd.compute(A,BB, matchingOne);

        arma::mat rotMat;

        if (dim == 2)
        {

            generateRandomRotationMatrix2D(rotMat);

        }
        else
        {

            generateRandomRotationMatrix3D(rotMat);


        }

        //  std::cout << "Rotation matrix: " << std::endl;

        //  std::cout << rotMat << std::endl;

        BB = rotMat * BB;

        arma::vec tt;
        generateRandomVector(tt, dim, 1.0);


        ShiftBy(BB,tt);



        // std::cout << "Cloud BB:" << std::endl;
        //  std::cout << BB << std::endl;



        //  std::cout << "Everything initilized " << std::endl;

        std::vector<int> matching;

        //   std::cout << A << std::endl;

        //    std::cout << BB << std::endl;
        std::ofstream debugS("/home/yury/LocalProjects/March19Project/debug/debugSingle.txt");
        res = nestedSCD.Compute_Isometry(A,BB, matching, totalTimeForSMIstd, totalTimeForBD, debugS);
        debugS.close();
        ratiod = res / perturbError;

        //   std::cout << "Obtained res" << std::endl;
        //std::vector<int> matchingZ;
        //double resd = nestedSCD.Compute_Isometry(A,BB, matchingZ, totalTimeForSMIstd, totalTimeForBD);



        arma::mat AA;

        std::cout << "Original matching" << std::endl;

        ExtraExperiments::printVector(justPermute);


        std::cout << "Matching found" << std::endl;

        //  std::cout << "A dimensinos : " << A.n_rows << " , " << A.n_cols << std::endl;
        //  std::cout << "Matching size: " << matching.size() << std::endl;
        ExtraExperiments::printVector(matching);

        ExtraExperiments::RearrangeCloud(A, AA, matching);




        arma::mat Rv;
        arma::vec tttv;


        findRotation(AA, BB,  Rv, tttv);

        arma::mat AAA = Rv*AA;

        ExtraExperiments::ShiftBy(AAA,tttv);

        realDist = 0;


        for (int i = 0; i < matching.size(); i++)
        {
            //   std::cout << "Comparing " << std::endl;
            //     std::cout << A.col(i) << std::endl;
            //   std::cout << BB.col(matching[i]) << std::endl;
            double distP = arma::norm(AAA.col(i)-BB.col(i), 2);
            //  double distP2 = arma::norm(A.col(i)-BB.col(matchingZ[i]), 2);

            //   std::cout << i << " --- " << matching[i] << " : "<< distP << std::endl;
            //   if (distP > curEpsion)
            //   {
            //       correctMatching = false;

            //  }
            realDist = std::max(realDist, distP);
            //  realDistTwo = std::max(realDistTwo, distP2);

            //  std::cout << "Matched : " << i << " --- " << justPermute[i] << std::endl;

        }

        std::cout << "RealDist: " << realDist << std::endl;
        std::cout << "Res: " << res << std::endl;
        std::cout << "Perturb error: " << perturbError << std::endl;


        ratiodistbd = realDist / res;

        //    std::cout << "Cloud AAA:" << std::endl;
        //     std::cout << AAA << std::endl;
        //         std::cout << "Cloud BB:" << std::endl;
//         std::cout << BB << std::endl;



        totalTimeForSMIstdA = (double) totalTimeForSMIstd;
        totalTimeForBDA = (double) totalTimeForBD;

        // std::cout << "Iteration: " << iter << " repeat: " << repeat << " ended" << std::endl;
        std::cout << "Iteration: " << iter << "curEpsilon: " << numEpsilon << " repeat: " << repeat << " ended" << std::endl;
        std::cout << "It returned : " << res << " , " << ratiod << " , " << realDist << ", " << ratiodistbd << " ," << totalTimeForSMIstdA << " , " << totalTimeForBDA << std::endl;



    }


    static void PerturbSingleNCDSmartComplicatedZ(int dim, int numPoints,  double curEpsion, int repeat, int numEpsilon, double & res, double & ratiod, double & realDist, double & ratiodistbd,
            double & totalTimeForSMIstdA, double & totalTimeForBDA)
    {
        long totalTimeForSMIstd;
        long totalTimeForBD;
        int iter = numPoints;
        std::cout << "Iteration: " << iter << "curEpsilon: " << numEpsilon << " repeat: " << repeat << " started" << std::endl;
        typedef BottleneckCalculator<Bipartie_Graph,GraphTreeControllerVector> bipartiebd;
        typedef bd_general_interface<BoostBottleneck,BoostBottleneck> boostBdCalc;
        typedef bd_general_interface<bipartiebd, bipartiebd> naiveBdCalc;

        Nested_Updated<PBI_Update<naiveBdCalc>, arma::mat> nestedSCD;

        std::vector<std::vector<double>> baseCloud;
        ExtraExperiments::generateRandomCloudInBoxNd(baseCloud, numPoints, 1.0, dim);

        std::vector<std::string> simpleStringVector;
        std::vector<std::vector<double>> perturbedCloud;
        double perturbError = ExtraExperiments::perturbCloud(baseCloud,perturbedCloud, curEpsion);

        arma::mat A;
        arma::mat B;




        vectorTools::vec_to_mat(baseCloud, A);
        vectorTools::vec_to_mat(perturbedCloud, B);

        //    std::cout << "Original A" << std::endl;
        //      std::cout << A << std::endl;
        //     std::cout << "Original B" << std::endl;
        //     std::cout << B << std::endl;



        std::vector<int> justPermute;
        for (int i = 0; i < numPoints; i++)
        {
            justPermute.push_back(i);

        }
        auto rng = std::default_random_engine {};

        std::shuffle(std::begin(justPermute), std::end(justPermute), rng);


        arma::mat BB;
        ExtraExperiments::RearrangeCloud(B, BB, justPermute);

        EuclideanBDCompute<bipartiebd> eucBd;

        std::vector<int> matchingOne;
        double purebd = eucBd.compute(A,BB, matchingOne);

        arma::mat rotMat;

        if (dim == 2)
        {

            generateRandomRotationMatrix2D(rotMat);

        }
        else
        {

            generateRandomRotationMatrix3D(rotMat);


        }

        //  std::cout << "Rotation matrix: " << std::endl;

        //  std::cout << rotMat << std::endl;

        BB = rotMat * BB;

        arma::vec tt;
        generateRandomVector(tt, dim, 1.0);


        ShiftBy(BB,tt);



        // std::cout << "Cloud BB:" << std::endl;
        //  std::cout << BB << std::endl;



        //  std::cout << "Everything initilized " << std::endl;

        std::vector<int> matching;

        //   std::cout << A << std::endl;

        //    std::cout << BB << std::endl;

        res = nestedSCD.Compute_Isometry(A,BB, matching, totalTimeForSMIstd, totalTimeForBD);
        ratiod = res / perturbError;

        //  std::cout << "Obtained res" << std::endl;
        //std::vector<int> matchingZ;
        //double resd = nestedSCD.Compute_Isometry(A,BB, matchingZ, totalTimeForSMIstd, totalTimeForBD);



        arma::mat AA;

        //std::cout << "EXperiments finnished " << std::endl;

        //  std::cout << "A dimensinos : " << A.n_rows << " , " << A.n_cols << std::endl;
        //  std::cout << "Matching size: " << matching.size() << std::endl;
        //  ExtraExperiments::printVector(matching);

        ExtraExperiments::RearrangeCloud(A, AA, matching);




        arma::mat Rv;
        arma::vec tttv;


        findRotation(AA, BB,  Rv, tttv);

        arma::mat AAA = Rv*AA;

        ExtraExperiments::ShiftBy(AAA,tttv);

        realDist = 0;


        for (int i = 0; i < matching.size(); i++)
        {
            //   std::cout << "Comparing " << std::endl;
            //     std::cout << A.col(i) << std::endl;
            //   std::cout << BB.col(matching[i]) << std::endl;
            double distP = arma::norm(AAA.col(i)-BB.col(i), 2);
            //  double distP2 = arma::norm(A.col(i)-BB.col(matchingZ[i]), 2);

            //   std::cout << i << " --- " << matching[i] << " : "<< distP << std::endl;
            //   if (distP > curEpsion)
            //   {
            //       correctMatching = false;

            //  }
            realDist = std::max(realDist, distP);
            //  realDistTwo = std::max(realDistTwo, distP2);

            //  std::cout << "Matched : " << i << " --- " << justPermute[i] << std::endl;

        }
        ratiodistbd = realDist / res;

        //    std::cout << "Cloud AAA:" << std::endl;
        //     std::cout << AAA << std::endl;
        //         std::cout << "Cloud BB:" << std::endl;
//         std::cout << BB << std::endl;



        totalTimeForSMIstdA = (double) totalTimeForSMIstd;
        totalTimeForBDA = (double) totalTimeForBD;

        // std::cout << "Iteration: " << iter << " repeat: " << repeat << " ended" << std::endl;
        std::cout << "Iteration: " << iter << "curEpsilon: " << numEpsilon << " repeat: " << repeat << " ended" << std::endl;
        std::cout << "It returned : " << res << " , " << ratiod << " , " << realDist << ", " << ratiodistbd << " ," << totalTimeForSMIstdA << " , " << totalTimeForBDA << std::endl;



    }

    static void printCloud(std::string fname, arma::mat & veccd)
    {
        std::ofstream out(fname);
        for (int i = 0; i < veccd.n_cols; i++)
        {
            for (int j = 0; j < veccd.n_rows; j++)
            {
                out << veccd.col(i)(j);
                if (j < veccd.n_rows - 1)
                {
                    out << ",";
                }



            }
            out << std::endl;



        }
        //    for (int i = 0; i < veccd.size(); i++)
        //   {
        //     out << veccd[i];
        //       if (i < veccd.size()-1)
        //     {
        //          out << ",";
        //     }
        // }
        out.close();


    }

    static void PerturbSingleNDPSmartComplicatedZ(int dim, int numPoints,  double curEpsion, int repeat, int numEpsilon, double & res, double & ratiod, double & realDist, double & ratiodistbd,
            double & totalTimeForSMIstdA, double & totalTimeForBDA, std::string pathOutA, std::string pathOutB)
    {
        long totalTimeForSMIstd;
        long totalTimeForBD;
        int iter = numPoints;
        std::cout << "Iteration: " << iter << "curEpsilon: " << numEpsilon << " repeat: " << repeat << " started" << std::endl;
        typedef BottleneckCalculator<Bipartie_Graph,GraphTreeControllerVector> bipartiebd;
        typedef bd_general_interface<BoostBottleneck,BoostBottleneck> boostBdCalc;
        typedef bd_general_interface<bipartiebd, bipartiebd> naiveBdCalc;

        Nested_Updated<PBI_Update<naiveBdCalc>, arma::mat> nestedSCD;

        std::vector<std::vector<double>> baseCloud;
        ExtraExperiments::generateRandomCloudInBoxNd(baseCloud, numPoints, 1.0, dim);

        std::vector<std::string> simpleStringVector;
        std::vector<std::vector<double>> perturbedCloud;
        double perturbError = ExtraExperiments::perturbCloud(baseCloud,perturbedCloud, curEpsion);

        arma::mat A;
        arma::mat B;




        vectorTools::vec_to_mat(baseCloud, A);
        vectorTools::vec_to_mat(perturbedCloud, B);

        //    std::cout << "Original A" << std::endl;
        //      std::cout << A << std::endl;
        //     std::cout << "Original B" << std::endl;
        //     std::cout << B << std::endl;



        std::vector<int> justPermute;
        for (int i = 0; i < numPoints; i++)
        {
            justPermute.push_back(i);

        }
        auto rng = std::default_random_engine {};

        std::shuffle(std::begin(justPermute), std::end(justPermute), rng);


        arma::mat BB;
        ExtraExperiments::RearrangeCloud(B, BB, justPermute);

        EuclideanBDCompute<bipartiebd> eucBd;

        std::vector<int> matchingOne;
        double purebd = eucBd.compute(A,BB, matchingOne);

        arma::mat rotMat;

        if (dim == 2)
        {

            generateRandomRotationMatrix2D(rotMat);

        }
        else
        {

            generateRandomRotationMatrix3D(rotMat);


        }

        //  std::cout << "Rotation matrix: " << std::endl;

        //  std::cout << rotMat << std::endl;

        BB = rotMat * BB;

        arma::vec tt;
        generateRandomVector(tt, dim, 1.0);


        ShiftBy(BB,tt);



        // std::cout << "Cloud BB:" << std::endl;
        //  std::cout << BB << std::endl;



        //  std::cout << "Everything initilized " << std::endl;

        std::vector<int> matching;

        //   std::cout << A << std::endl;

        //    std::cout << BB << std::endl;

        ExtraExperiments::printCloud(pathOutA, A);
        ExtraExperiments::printCloud(pathOutB, BB);


        res = nestedSCD.Compute_Isometry(A,BB, matching, totalTimeForSMIstd, totalTimeForBD);
        ratiod = res / perturbError;

        //  std::cout << "Obtained res" << std::endl;
        //std::vector<int> matchingZ;
        //double resd = nestedSCD.Compute_Isometry(A,BB, matchingZ, totalTimeForSMIstd, totalTimeForBD);



        arma::mat AA;

        //std::cout << "EXperiments finnished " << std::endl;

        //  std::cout << "A dimensinos : " << A.n_rows << " , " << A.n_cols << std::endl;
        //  std::cout << "Matching size: " << matching.size() << std::endl;
        //  ExtraExperiments::printVector(matching);

        ExtraExperiments::RearrangeCloud(A, AA, matching);




        arma::mat Rv;
        arma::vec tttv;


        findRotation(AA, BB,  Rv, tttv);

        arma::mat AAA = Rv*AA;

        ExtraExperiments::ShiftBy(AAA,tttv);

        realDist = 0;


        for (int i = 0; i < matching.size(); i++)
        {
            //   std::cout << "Comparing " << std::endl;
            //     std::cout << A.col(i) << std::endl;
            //   std::cout << BB.col(matching[i]) << std::endl;
            double distP = arma::norm(AAA.col(i)-BB.col(i), 2);
            //  double distP2 = arma::norm(A.col(i)-BB.col(matchingZ[i]), 2);

            //   std::cout << i << " --- " << matching[i] << " : "<< distP << std::endl;
            //   if (distP > curEpsion)
            //   {
            //       correctMatching = false;

            //  }
            realDist = std::max(realDist, distP);
            //  realDistTwo = std::max(realDistTwo, distP2);

            //  std::cout << "Matched : " << i << " --- " << justPermute[i] << std::endl;

        }
        ratiodistbd = realDist / res;

        //    std::cout << "Cloud AAA:" << std::endl;
        //     std::cout << AAA << std::endl;
        //         std::cout << "Cloud BB:" << std::endl;
//         std::cout << BB << std::endl;



        totalTimeForSMIstdA = (double) totalTimeForSMIstd;
        totalTimeForBDA = (double) totalTimeForBD;

        // std::cout << "Iteration: " << iter << " repeat: " << repeat << " ended" << std::endl;
        std::cout << "Iteration: " << iter << "curEpsilon: " << numEpsilon << " repeat: " << repeat << " ended" << std::endl;
        std::cout << "It returned : " << res << " , " << ratiod << " , " << realDist << ", " << ratiodistbd << " ," << totalTimeForSMIstdA << " , " << totalTimeForBDA << std::endl;



    }



    static void saveToCsv(std::string fname, std::vector<double> & veccd)
    {
        std::ofstream out(fname);
        for (int i = 0; i < veccd.size(); i++)
        {
            out << veccd[i];
            if (i < veccd.size()-1)
            {
                out << ",";
            }
        }
        out.close();


    }

    static void NCDZSuperFast(std::string path, int bigI, int bigJ, int dim, int numPoints,  double curEpsion, int repeat, int numEpsilon)
    {

        std::vector<double> maxx(6, std::numeric_limits<double>::min());
        std::vector<double> minn(6, std::numeric_limits<double>::max());
        std::vector<double> avg (6, 0 );
        for (int i = 0; i < repeat; i++)
        {
            std::vector<double> params(6,0);
            ExtraExperiments::PerturbSingleNCDSmartComplicatedZ(dim, numPoints,  curEpsion, i, numEpsilon, params[0], params[1],params[2],params[3], params[4], params[5]);
            for (int j = 0; j < 6; j++)
            {
                maxx[j] = std::max(maxx[j], params[j]);
                minn[j] = std::min(minn[j], params[j]);
                avg[j] =  avg[j] + params[j] / repeat;

            }
        }
        for (int i = 0; i < 6; i++)
        {
            saveToCsv(path + "Avg/" + std::to_string(bigI) + "_" +std::to_string(bigJ) + ".csv", avg);
            saveToCsv(path + "Max/" + std::to_string(bigI) + "_" + std::to_string(bigJ) + ".csv", maxx);
            saveToCsv(path + "Min/" + std::to_string(bigI) + "_" + std::to_string(bigJ) + ".csv", minn);


        }


    }

     static void MatchingExperimentsTime(int threads, int dim, int numPoints, int stepPoints, double epsilon, int repeat, std::string outputCore, std::string outputDist)
    {

        int stepsEps = stepPoints;
        std::vector<std::vector<double>> maxx(1, std::vector<double>(stepsEps, std::numeric_limits<double>::min()));
        std::vector<std::vector<double>> minn(1, std::vector<double>(stepsEps, std::numeric_limits<double>::max()));
        std::vector<std::vector<double>> avg(1, std::vector<double>(stepsEps, 0));
        std::vector<std::vector<std::vector<double>>> params(stepsEps,std::vector<std::vector<double>>(repeat, std::vector<double>(6, 0)));
        //  std::vector<double> minn(6, std::numeric_limits<double>::max());
        //   std::vector<double> avg (6, 0 );

        //     std::vector<int> numPointsVector;
        //     BS::thread_pool pool(threads);
        //     for (int i = 0; i < stepsPts; i++)
        //     {
        //       numPointsVector.push_back((i+1)*stepSizePts);
        //  double lel = (double) (i+1)*stepSizePts;
//            ptsVec.push_back(lel);
        //  }

        std::vector<double> epsPointsVector;
        for (int i = 0; i < stepPoints; i++)
        {
            epsPointsVector.push_back((i+1)*numPoints);
        }

        BS::thread_pool pool(threads);

        //    for (int k = 2; k < stepsEps; k++)
        //     {

        //    int i = 4;
        // //        ExtraExperiments::NCDZSuperFast(outputPath, i, k, dim,  numPointsVector[i], epsPointsVector[k], repeat, k);

        //   pool.push_task(ExtraExperiments::NCDZSuperFast,outputPath, k, i, dim,  numPointsVector[i], epsPointsVector[k], repeat, k);

        //      }




        for (int k = 0; k < stepsEps; k++)
        {
            for (int b = 0; b < repeat; b++)
            {
                //!     static void NCDZSuperFast(std::string path, int bigI, int bigJ, int dim, int numPoints,  double curEpsion, int repeat, int numEpsilon)

           //     int numPoints = numPointsVector[i];
                double curEpsion =  epsPointsVector[k];
                int numEpsilon = k;
                pool.push_task(ExtraExperiments::PerturbSingleNDPSmartComplicatedZ,dim, curEpsion,  epsilon, b, numEpsilon,
                               std::ref(params[k][b][0]), std::ref(params[k][b][1]),std::ref(params[k][b][2]),
                               std::ref(params[k][b][3]), std::ref(params[k][b][4]), std::ref(params[k][b][5]),
                               outputCore  + "/" +  std::to_string(k) + "/A/out_" + std::to_string(b) + ".csv"  ,
                               outputCore + "/" +  std::to_string(k) + "/B/out_" + std::to_string(b) + ".csv"  );

                //         ExtraExperiments::PerturbSingleNCDSmartComplicatedZ(dim, numPoints,  curEpsion, b, numEpsilon,
                //        (params[i][k][0]),(params[i][k][1]),(params[i][k][2]),
                //         (params[i][k][3]), (params[i][k][4]), (params[i][k][5]));

            }


            //    ExtraExperiments::NCDZSuperFast(outputPath, i, k, dim,  numPointsVector[i], epsPointsVector[k], repeat, k);



        }



        pool.wait_for_tasks();



        for (int k = 0; k < stepsEps; k++)
            {

                for (int b = 0; b < repeat; b++)
                {

                    maxx[0][k] = std::max(maxx[0][k], params[k][b][4] + params[k][b][5] );
                    minn[0][k] = std::min(minn[0][k], params[k][b][4] + params[k][b][5] );
                    avg[0][k] =  avg[0][k] + (params[k][b][4] + params[k][b][5]) / repeat;
                    }


            }

         //   std::string outputDist, std::string outputFinalMatch

            PrintExperimentsModi(maxx[0], minn[0], avg[0], epsPointsVector, outputDist );
        //    PrintExperimentsModi(maxx[1], minn[1], avg[1], epsPointsVector, outputFinalMatch );

          //!  PrintExperimentsModi(std::vector<double> & maxx, std::vector<double> & minn, std::vector<double> & avg, std::vector<double> & epsilons, std::string outputpath )

     //   std::string path = outputPath;

//        for (int i = 0; i < stepsPts; i++)
//        {
//
//            for (int k = 0; k < stepsEps; k++)
//            {
//                int bigI = i;
//                int bigJ = k;
//
//                //     for (int i = 0; i < 6; i++)
//                //    {
//                saveToCsv(path + "Avg/" + std::to_string(bigI) + "_" +std::to_string(bigJ) + ".csv", avg[i][k]);
//                saveToCsv(path + "Max/" + std::to_string(bigI) + "_" + std::to_string(bigJ) + ".csv", maxx[i][k]);
//                saveToCsv(path + "Min/" + std::to_string(bigI) + "_" + std::to_string(bigJ) + ".csv", minn[i][k]);
//
//
//                //     }
//            }
//        }




    }

    static void MatchingExperimentsMatch(int threads, int dim, int numPoints, int stepsEps, double stepSizeEps, int repeat, std::string outputCore, std::string outputDist, std::string outputFinalMatch)
    {

        std::vector<std::vector<double>> maxx(2, std::vector<double>(stepsEps, std::numeric_limits<double>::min()));
        std::vector<std::vector<double>> minn(2, std::vector<double>(stepsEps, std::numeric_limits<double>::max()));
        std::vector<std::vector<double>> avg(2, std::vector<double>(stepsEps, 0));
        std::vector<std::vector<std::vector<double>>> params(stepsEps,std::vector<std::vector<double>>(repeat, std::vector<double>(6, 0)));
        //  std::vector<double> minn(6, std::numeric_limits<double>::max());
        //   std::vector<double> avg (6, 0 );

        //     std::vector<int> numPointsVector;
        //     BS::thread_pool pool(threads);
        //     for (int i = 0; i < stepsPts; i++)
        //     {
        //       numPointsVector.push_back((i+1)*stepSizePts);
        //  double lel = (double) (i+1)*stepSizePts;
//            ptsVec.push_back(lel);
        //  }

        std::vector<double> epsPointsVector;
        for (int i = 0; i < stepsEps; i++)
        {
            epsPointsVector.push_back((i+1)*stepSizeEps);
        }

        BS::thread_pool pool(threads);

        //    for (int k = 2; k < stepsEps; k++)
        //     {

        //    int i = 4;
        // //        ExtraExperiments::NCDZSuperFast(outputPath, i, k, dim,  numPointsVector[i], epsPointsVector[k], repeat, k);

        //   pool.push_task(ExtraExperiments::NCDZSuperFast,outputPath, k, i, dim,  numPointsVector[i], epsPointsVector[k], repeat, k);

        //      }




        for (int k = 0; k < stepsEps; k++)
        {
            for (int b = 0; b < repeat; b++)
            {
                //!     static void NCDZSuperFast(std::string path, int bigI, int bigJ, int dim, int numPoints,  double curEpsion, int repeat, int numEpsilon)

           //     int numPoints = numPointsVector[i];
                double curEpsion =  epsPointsVector[k];
                int numEpsilon = k;
                pool.push_task(ExtraExperiments::PerturbSingleNDPSmartComplicatedZ,dim, numPoints,  curEpsion, b, numEpsilon,
                               std::ref(params[k][b][0]), std::ref(params[k][b][1]),std::ref(params[k][b][2]),
                               std::ref(params[k][b][3]), std::ref(params[k][b][4]), std::ref(params[k][b][5]),
                               outputCore  + "/" +  std::to_string(k) + "/A/out_" + std::to_string(b) + ".csv"  ,
                               outputCore + "/" +  std::to_string(k) + "/B/out_" + std::to_string(b) + ".csv"  );

                //         ExtraExperiments::PerturbSingleNCDSmartComplicatedZ(dim, numPoints,  curEpsion, b, numEpsilon,
                //        (params[i][k][0]),(params[i][k][1]),(params[i][k][2]),
                //         (params[i][k][3]), (params[i][k][4]), (params[i][k][5]));

            }


            //    ExtraExperiments::NCDZSuperFast(outputPath, i, k, dim,  numPointsVector[i], epsPointsVector[k], repeat, k);



        }



        pool.wait_for_tasks();



        for (int k = 0; k < stepsEps; k++)
            {

                for (int b = 0; b < repeat; b++)
                {

                    maxx[0][k] = std::max(maxx[0][k], params[k][b][0]);
                    minn[0][k] = std::min(minn[0][k], params[k][b][0]);
                    avg[0][k] =  avg[0][k] + params[k][b][0] / repeat;

                    maxx[1][k] = std::max(maxx[1][k], params[k][b][2]);
                    minn[1][k] = std::min(minn[1][k], params[k][b][2]);
                    avg[1][k] =  avg[1][k] + params[k][b][2] / repeat;

                    }


            }

         //   std::string outputDist, std::string outputFinalMatch

            PrintExperimentsModi(maxx[0], minn[0], avg[0], epsPointsVector, outputDist );
            PrintExperimentsModi(maxx[1], minn[1], avg[1], epsPointsVector, outputFinalMatch );

          //!  PrintExperimentsModi(std::vector<double> & maxx, std::vector<double> & minn, std::vector<double> & avg, std::vector<double> & epsilons, std::string outputpath )

     //   std::string path = outputPath;

//        for (int i = 0; i < stepsPts; i++)
//        {
//
//            for (int k = 0; k < stepsEps; k++)
//            {
//                int bigI = i;
//                int bigJ = k;
//
//                //     for (int i = 0; i < 6; i++)
//                //    {
//                saveToCsv(path + "Avg/" + std::to_string(bigI) + "_" +std::to_string(bigJ) + ".csv", avg[i][k]);
//                saveToCsv(path + "Max/" + std::to_string(bigI) + "_" + std::to_string(bigJ) + ".csv", maxx[i][k]);
//                saveToCsv(path + "Min/" + std::to_string(bigI) + "_" + std::to_string(bigJ) + ".csv", minn[i][k]);
//
//
//                //     }
//            }
//        }




    }




    static void MatchingExperiments2DBBA(int threads, int dim, int stepsPts, int stepSizePts, int stepsEps, double stepSizeEps, int repeat, std::string outputPath)
    {

        std::vector<std::vector<std::vector<double>>> maxx(stepsPts, std::vector<std::vector<double>>(stepsEps, std::vector<double>(6, std::numeric_limits<double>::min())));
        std::vector<std::vector<std::vector<double>>> minn(stepsPts, std::vector<std::vector<double>>(stepsEps, std::vector<double>(6, std::numeric_limits<double>::max())));
        std::vector<std::vector<std::vector<double>>> avg(stepsPts, std::vector<std::vector<double>>(stepsEps, std::vector<double>(6, 0)));
        std::vector<std::vector<std::vector<double>>> params(stepsPts, std::vector<std::vector<double>>(stepsEps, std::vector<double>(6, 0)));
        //  std::vector<double> minn(6, std::numeric_limits<double>::max());
        //   std::vector<double> avg (6, 0 );

        std::vector<int> numPointsVector;
        BS::thread_pool pool(threads);
        for (int i = 0; i < stepsPts; i++)
        {
            numPointsVector.push_back((i+1)*stepSizePts);
            //  double lel = (double) (i+1)*stepSizePts;
//            ptsVec.push_back(lel);
        }

        std::vector<double> epsPointsVector;
        for (int i = 0; i < stepsEps; i++)
        {
            epsPointsVector.push_back((i+1)*stepSizeEps);
        }

//        BS::thread_pool pool(threads);

        //    for (int k = 2; k < stepsEps; k++)
        //     {

        //    int i = 4;
        // //        ExtraExperiments::NCDZSuperFast(outputPath, i, k, dim,  numPointsVector[i], epsPointsVector[k], repeat, k);

        //   pool.push_task(ExtraExperiments::NCDZSuperFast,outputPath, k, i, dim,  numPointsVector[i], epsPointsVector[k], repeat, k);

        //      }


        for (int i = 0; i < stepsPts; i++)
        {

            for (int k = 0; k < stepsEps; k++)
            {
                for (int b = 0; b < repeat; b++)
                {
                    //!     static void NCDZSuperFast(std::string path, int bigI, int bigJ, int dim, int numPoints,  double curEpsion, int repeat, int numEpsilon)

                    int numPoints = numPointsVector[i];
                    double curEpsion =  epsPointsVector[k];
                    int numEpsilon = k;
                    pool.push_task(ExtraExperiments::PerturbSingleNCDSmartComplicatedZ,dim, numPoints,  curEpsion, b, numEpsilon,
                                   std::ref(params[i][k][0]), std::ref(params[i][k][1]),std::ref(params[i][k][2]),
                                   std::ref(params[i][k][3]), std::ref(params[i][k][4]), std::ref(params[i][k][5]));

                    //         ExtraExperiments::PerturbSingleNCDSmartComplicatedZ(dim, numPoints,  curEpsion, b, numEpsilon,
                    //        (params[i][k][0]),(params[i][k][1]),(params[i][k][2]),
                    //         (params[i][k][3]), (params[i][k][4]), (params[i][k][5]));

                }


                //    ExtraExperiments::NCDZSuperFast(outputPath, i, k, dim,  numPointsVector[i], epsPointsVector[k], repeat, k);



            }


        }
        pool.wait_for_tasks();

        for (int i = 0; i < stepsPts; i++)
        {

            for (int k = 0; k < stepsEps; k++)
            {

                for (int b = 0; b < 6; b++)
                {
                    maxx[i][k][b] = std::max(maxx[i][k][b], params[i][k][b]);
                    minn[i][k][b] = std::min(minn[i][k][b], params[i][k][b]);
                    avg[i][k][b] =  avg[i][k][b] + params[i][k][b] / repeat;

                }
            }
        }
        std::string path = outputPath;

        for (int i = 0; i < stepsPts; i++)
        {

            for (int k = 0; k < stepsEps; k++)
            {
                int bigI = i;
                int bigJ = k;

                //     for (int i = 0; i < 6; i++)
                //    {
                saveToCsv(path + "Avg/" + std::to_string(bigI) + "_" +std::to_string(bigJ) + ".csv", avg[i][k]);
                saveToCsv(path + "Max/" + std::to_string(bigI) + "_" + std::to_string(bigJ) + ".csv", maxx[i][k]);
                saveToCsv(path + "Min/" + std::to_string(bigI) + "_" + std::to_string(bigJ) + ".csv", minn[i][k]);


                //     }
            }
        }




    }


    static void MatchingExperiments2DBB(int threads, int dim, int stepsPts, int stepSizePts, int stepsEps, double stepSizeEps, int repeat, std::string outputPath)
    {

        std::vector<int> numPointsVector;
        BS::thread_pool pool(threads);
        for (int i = 0; i < stepsPts; i++)
        {
            numPointsVector.push_back((i+1)*stepSizePts);
            //  double lel = (double) (i+1)*stepSizePts;
//            ptsVec.push_back(lel);
        }

        std::vector<double> epsPointsVector;
        for (int i = 0; i < stepsEps; i++)
        {
            epsPointsVector.push_back((i+1)*stepSizeEps);
        }

//        BS::thread_pool pool(threads);

        //    for (int k = 2; k < stepsEps; k++)
        //     {

        //    int i = 4;
        // //        ExtraExperiments::NCDZSuperFast(outputPath, i, k, dim,  numPointsVector[i], epsPointsVector[k], repeat, k);

        //   pool.push_task(ExtraExperiments::NCDZSuperFast,outputPath, k, i, dim,  numPointsVector[i], epsPointsVector[k], repeat, k);

        //      }


        for (int i = 0; i < stepsPts; i++)
        {

            for (int k = 0; k < stepsEps; k++)
            {


                //    ExtraExperiments::NCDZSuperFast(outputPath, i, k, dim,  numPointsVector[i], epsPointsVector[k], repeat, k);

                pool.push_task(ExtraExperiments::NCDZSuperFast,outputPath, i, k, dim,  numPointsVector[i], epsPointsVector[k], repeat, k);

            }


        }


        pool.wait_for_tasks();


    }






    static void PerturbSingleNCDSmartComplicated(int dim, int numPoints,  double curEpsion, int repeat, int numEpsilon, double & res, double & ratiod, double & realDist, double & ratiodistbd,
            long & totalTimeForSMIstd, long & totalTimeForBD)
    {
        int iter = numEpsilon;
        std::cout << "Iteration: " << iter << " repeat: " << repeat << " started" << std::endl;
        typedef BottleneckCalculator<Bipartie_Graph,GraphTreeController> bipartiebd;
        typedef bd_general_interface<BoostBottleneck,BoostBottleneck> boostBdCalc;
        typedef bd_general_interface<bipartiebd, bipartiebd> naiveBdCalc;

        Nested_Updated<NCD_Smart_Update<naiveBdCalc>, arma::mat> nestedSCD;

        std::vector<std::vector<double>> baseCloud;
        ExtraExperiments::generateRandomCloudInBoxNd(baseCloud, numPoints, 1.0, dim);

        std::vector<std::string> simpleStringVector;
        std::vector<std::vector<double>> perturbedCloud;
        double perturbError = ExtraExperiments::perturbCloud(baseCloud,perturbedCloud, curEpsion);

        arma::mat A;
        arma::mat B;




        vectorTools::vec_to_mat(baseCloud, A);
        vectorTools::vec_to_mat(perturbedCloud, B);

        //    std::cout << "Original A" << std::endl;
        //      std::cout << A << std::endl;
        //     std::cout << "Original B" << std::endl;
        //     std::cout << B << std::endl;



        std::vector<int> justPermute;
        for (int i = 0; i < numPoints; i++)
        {
            justPermute.push_back(i);

        }
        auto rng = std::default_random_engine {};

        std::shuffle(std::begin(justPermute), std::end(justPermute), rng);


        arma::mat BB;
        ExtraExperiments::RearrangeCloud(B, BB, justPermute);

        EuclideanBDCompute<bipartiebd> eucBd;

        std::vector<int> matchingOne;


        double purebd = eucBd.compute(A,BB, matchingOne);

        arma::mat rotMat;

        generateRandomRotationMatrix2D(rotMat);

        BB = rotMat * BB;

        arma::vec tt;
        generateRandomVector(tt, 2, 1.0);


        ShiftBy(BB,tt);



        // std::cout << "Cloud BB:" << std::endl;
        //  std::cout << BB << std::endl;



        std::vector<int> matching;
        res = nestedSCD.Compute_Isometry(A,BB, matching, totalTimeForSMIstd, totalTimeForBD);
        ratiod = res / perturbError;
        //std::vector<int> matchingZ;
        //double resd = nestedSCD.Compute_Isometry(A,BB, matchingZ, totalTimeForSMIstd, totalTimeForBD);



        arma::mat AA;

        //std::cout << "EXperiments finnished " << std::endl;

        //  std::cout << "A dimensinos : " << A.n_rows << " , " << A.n_cols << std::endl;
        //  std::cout << "Matching size: " << matching.size() << std::endl;
        //  ExtraExperiments::printVector(matching);

        ExtraExperiments::RearrangeCloud(A, AA, matching);




        arma::mat Rv;
        arma::vec tttv;


        findRotation(AA, BB,  Rv, tttv);

        arma::mat AAA = Rv*AA;

        ExtraExperiments::ShiftBy(AAA,tttv);

        realDist = 0;


        for (int i = 0; i < matching.size(); i++)
        {
            //   std::cout << "Comparing " << std::endl;
            //     std::cout << A.col(i) << std::endl;
            //   std::cout << BB.col(matching[i]) << std::endl;
            double distP = arma::norm(AAA.col(i)-BB.col(i), 2);
            //  double distP2 = arma::norm(A.col(i)-BB.col(matchingZ[i]), 2);

            //   std::cout << i << " --- " << matching[i] << " : "<< distP << std::endl;
            //   if (distP > curEpsion)
            //   {
            //       correctMatching = false;

            //  }
            realDist = std::max(realDist, distP);
            //  realDistTwo = std::max(realDistTwo, distP2);

            //  std::cout << "Matched : " << i << " --- " << justPermute[i] << std::endl;

        }
        ratiodistbd = realDist / res;

        //    std::cout << "Cloud AAA:" << std::endl;
        //     std::cout << AAA << std::endl;
        //         std::cout << "Cloud BB:" << std::endl;
//         std::cout << BB << std::endl;




        std::cout << "Iteration: " << iter << " repeat: " << repeat << " ended" << std::endl;



    }

    static void PerturbSingleNCDSmartSimple(int dim, int numPoints,  double curEpsion, int repeat, int numEpsilon, double & res, double & ratiod, double & realDist, double & ratiodistbd,
                                            long & totalTimeForSMIstd, long & totalTimeForBD)
    {
        int iter = numEpsilon;
        std::cout << "Iteration: " << iter << " repeat: " << repeat << " started" << std::endl;
        typedef BottleneckCalculator<Bipartie_Graph,GraphTreeController> bipartiebd;
        typedef bd_general_interface<BoostBottleneck,BoostBottleneck> boostBdCalc;
        typedef bd_general_interface<bipartiebd, bipartiebd> naiveBdCalc;

        Nested_Updated<NCD_Smart_Update<naiveBdCalc>, arma::mat> nestedSCD;

        std::vector<std::vector<double>> baseCloud;
        ExtraExperiments::generateRandomCloudInBoxNd(baseCloud, numPoints, 1.0, dim);

        std::vector<std::string> simpleStringVector;
        std::vector<std::vector<double>> perturbedCloud;
        double perturbError = ExtraExperiments::perturbCloud(baseCloud,perturbedCloud, curEpsion);

        arma::mat A;
        arma::mat B;




        vectorTools::vec_to_mat(baseCloud, A);
        vectorTools::vec_to_mat(perturbedCloud, B);

        //   std::cout << "Original A" << std::endl;
        std::cout << A << std::endl;
        // std::cout << "Original B" << std::endl;
        std::cout << B << std::endl;



        std::vector<int> justPermute;
        for (int i = 0; i < numPoints; i++)
        {
            justPermute.push_back(i);

        }
        auto rng = std::default_random_engine {};

        std::shuffle(std::begin(justPermute), std::end(justPermute), rng);


        arma::mat BB;
        ExtraExperiments::RearrangeCloud(B, BB, justPermute);

        //      std::cout << "Cloud BB:" << std::endl;
        //   std::cout << BB << std::endl;

        EuclideanBDCompute<bipartiebd> eucBd;

        std::vector<int> matchingOne;
        double purebd = eucBd.compute(A,BB, matchingOne);

        std::vector<int> matching;
        res = nestedSCD.Compute_Isometry(A,BB, matching, totalTimeForSMIstd, totalTimeForBD);
        ratiod = res / perturbError;
        //std::vector<int> matchingZ;
        //double resd = nestedSCD.Compute_Isometry(A,BB, matchingZ, totalTimeForSMIstd, totalTimeForBD);

        realDist = 0;


        for (int i = 0; i < matching.size(); i++)
        {
            //   std::cout << "Comparing " << std::endl;
            //     std::cout << A.col(i) << std::endl;
            //   std::cout << BB.col(matching[i]) << std::endl;
            double distP = arma::norm(A.col(i)-BB.col(matching[i]), 2);
            //  double distP2 = arma::norm(A.col(i)-BB.col(matchingZ[i]), 2);

            //   std::cout << i << " --- " << matching[i] << " : "<< distP << std::endl;
            //   if (distP > curEpsion)
            //   {
            //       correctMatching = false;

            //  }
            realDist = std::max(realDist, distP);
            //  realDistTwo = std::max(realDistTwo, distP2);

            //  std::cout << "Matched : " << i << " --- " << justPermute[i] << std::endl;

        }
        ratiodistbd = realDist / res;

        //   arma::mat AA;

        //   std::cout << "EXperiments finnished " << std::endl;

        //  std::cout << "A dimensinos : " << A.n_rows << " , " << A.n_cols << std::endl;
        //  std::cout << "Matching size: " << matching.size() << std::endl;
        //  ExtraExperiments::printVector(matching);

        //  ExtraExperiments::RearrangeCloud(A, AA, matching);



        //   std::cout << "Cloud AA:" << std::endl;
        //   std::cout << AA << std::endl;

        //    std::cout << "Iteration: " << iter << " repeat: " << repeat << " ended" << std::endl;



    }



    static void PerturbSingleNCDSmart(int dim, int numPoints,  double curEpsion, int repeat, int numEpsilon)
    {
        typedef BottleneckCalculator<Bipartie_Graph,GraphTreeController> bipartiebd;
        typedef bd_general_interface<BoostBottleneck,BoostBottleneck> boostBdCalc;
        typedef bd_general_interface<bipartiebd, bipartiebd> naiveBdCalc;

        Nested_Updated<NCD_Smart_Update<naiveBdCalc>, arma::mat> nestedSCD;
        Nested_Updated<NCD_Smart_Update<boostBdCalc>, arma::mat> boostSCD;

        std::vector<std::vector<double>> baseCloud;
        ExtraExperiments::generateRandomCloudInBoxNd(baseCloud, numPoints, 1.0, dim);


        std::vector<std::string> simpleStringVector;
        std::vector<std::vector<double>> perturbedCloud;
        double perturbError = ExtraExperiments::perturbCloud(baseCloud,perturbedCloud, curEpsion);

        arma::mat A;
        arma::mat B;



        vectorTools::vec_to_mat(baseCloud, A);
        vectorTools::vec_to_mat(perturbedCloud, B);

        auto rng = std::default_random_engine {};

        std::vector<int> justPermute;
        for (int i = 0; i < numPoints; i++)
        {
            justPermute.push_back(i);

        }
        std::shuffle(std::begin(justPermute), std::end(justPermute), rng);

        arma::mat BB;
        ExtraExperiments::RearrangeCloud(B, BB, justPermute);

        std::cout << "just permute: " << std::endl;
        ExtraExperiments::printVector(justPermute);

        long totalTimeForSMIstd;
        long totalTimeForBD;

        std::cout << A << std::endl;
        std::cout << BB << std::endl;

        std::vector<int> matchA;
        std::vector<int> matchB;
        bipartiebd uw;
        uw.bottleneck_between_clouds(A,BB,matchA);
        BoostBottleneck zd;
        zd.bottleneck_between_clouds(A,BB,matchB);
        std::cout << "normal bruteforce matching: " << std::endl;
        ExtraExperiments::printVector(matchA);
        std::cout << "boost bruteforce matching:  " << std::endl;
        ExtraExperiments::printVector(matchB);


        std::string folder = "/home/yury/LocalProjects/dbg/MatricesTwo/";
        std::string saveFileA = "A" + std::to_string(numEpsilon) +  ".mat";
        std::string saveFileB = "B" + std::to_string(numEpsilon) + ".mat";
        std::string saveFileBB = "BB" + std::to_string(numEpsilon) + ".mat";

        A.save(folder + saveFileA, arma::arma_ascii);
        B.save(folder + saveFileB, arma::arma_ascii);
        BB.save(folder + saveFileBB, arma::arma_ascii);


        std::vector<std::vector<double>> distMat;

        ExtraExperiments::computeDistanceMatrix(A, B, distMat);
        ExtraExperiments::printMatrix(distMat);



        EuclideanBDCompute<bipartiebd> eucBd;

        std::vector<int> matchingOne;
        double purebd = eucBd.compute(A,BB, matchingOne);


        std::cout << "Cloud A: " << std::endl;

        std::cout << A << std::endl;

        std::cout << "Cloud B: " << std::endl;

        std::cout << B << std::endl;

        std::cout << "Euclidean BD: " << purebd << std::endl;


        std::vector<int> matching;
        std::cout << "Launching isometry computations" << std::endl;
        double res = boostSCD.Compute_Isometry(A,BB, matching, totalTimeForSMIstd, totalTimeForBD);
        std::vector<int> matchingZ;
        double resd = nestedSCD.Compute_Isometry(A,BB, matchingZ, totalTimeForSMIstd, totalTimeForBD);


        std::cout << "Computations launched " << std::endl;

        std::cout << "Boost obtained: " << std::endl;
        ExtraExperiments::printVector(matching);

        std::cout << "Res is : " << res << " | ratio: " << res/purebd << std::endl;

        double realDist = 0;
        double realDistTwo = 0;

        for (int i = 0; i < matching.size(); i++)
        {
            //   std::cout << "Comparing " << std::endl;
            //     std::cout << A.col(i) << std::endl;
            //   std::cout << BB.col(matching[i]) << std::endl;
            double distP = arma::norm(A.col(i)-BB.col(matching[i]), 2);
            double distP2 = arma::norm(A.col(i)-BB.col(matchingZ[i]), 2);

            std::cout << i << " --- " << matching[i] << " : "<< distP << std::endl;
            //   if (distP > curEpsion)
            //   {
            //       correctMatching = false;

            //  }
            realDist = std::max(realDist, distP);
            realDistTwo = std::max(realDistTwo, distP2);

            //  std::cout << "Matched : " << i << " --- " << justPermute[i] << std::endl;

        }
        std::cout << "Realdist was : " << realDist << std::endl;
        std::cout << "RealdistTwo was : " << realDistTwo << std::endl;



    }

    static void PerturbExperimentsSingleNew(int dim, int numPoints,  double curEpsion,
                                            int repeat, int numEpsilon, double & res, double & ratiod, bool & correctMatching, double & realDist)
    {

        std::cout << "Epsilon: " << numEpsilon << " repeat: " << repeat << " started" << std::endl;

        typedef BottleneckCalculator<Bipartie_Graph,GraphTreeController> bipartiebd;
        typedef bd_general_interface<bipartiebd, bipartiebd> naiveBdCalc;
        typedef bd_general_interface<BoostBottleneck,BoostBottleneck> boostBdCalc;

        EuclideanBDCompute<bipartiebd> eucBd;
        //!  Nested_Base<SCD, naiveBdCalc, std::pair<arma::mat, std::vector<std::vector<bool>>>> nestedNaive;
        //Nested_Base<PBI, naiveBdCalc, arma::mat> nestedNaive;
        Nested_Base<SCD,  naiveBdCalc, std::pair<arma::mat, std::vector<std::vector<bool>>>> nestedSCD;
        Nested_Base<SCD,  boostBdCalc, std::pair<arma::mat, std::vector<std::vector<bool>>>> boostSCD;

        std::vector<std::vector<double>> baseCloud;
        ExtraExperiments::generateRandomCloudInBoxNd(baseCloud, numPoints, 1.0, dim);


        std::vector<std::string> simpleStringVector;
        std::vector<std::vector<double>> perturbedCloud;
        ExtraExperiments::perturbCloud(baseCloud,perturbedCloud, curEpsion);
        auto rng = std::default_random_engine {};

        std::vector<int> justPermute;
        for (int i = 0; i < numPoints; i++)
        {
            justPermute.push_back(i);

        }
        std::shuffle(std::begin(justPermute), std::end(justPermute), rng);

        arma::mat A;
        arma::mat B;


        vectorTools::vec_to_mat(baseCloud, A);
        vectorTools::vec_to_mat(perturbedCloud, B);
        arma::mat BB;
        ExtraExperiments::RearrangeCloud(B, BB, justPermute);
        arma::mat C = BB;

        arma::mat rotMat;
        ExtraExperiments::generateRandomRotationMatrix2D(rotMat);
        BB = rotMat * BB;



        //    std::cout << "Matrix A" << std::endl;
        //    std::cout << A << std::endl;

        //    std::cout << "Matrix B " << std::endl;
        //    std::cout << B << std::endl;

        //    std::cout << "Matrix BB " << std::endl;
        //   std::cout << BB << std::endl;


        std::vector<std::vector<bool>> graphA(A.n_cols+1, std::vector<bool>(A.n_cols+1, true));
        std::vector<std::vector<bool>> graphB(B.n_cols+1, std::vector<bool>(B.n_cols+1, true));



        //     ExtraExperiments::randomGraphGenerator(graphA);
        //    graphB = graphA;
//!        ExtraExperiments::randomGraphGenerator(graphB);

        long totalTimeForSMIstd;
        long totalTimeForBD;
        std::vector<int> matching;

        auto firstP = std::make_pair(A, graphA);
        auto secondP = std::make_pair(BB, graphB);

        res = boostSCD.Compute_Isometry(firstP,secondP, matching, totalTimeForSMIstd, totalTimeForBD);


//        long totalTimeForSMIstd;
        //    long totalTimeForBD;
        std::vector<int> matchingtwo;

        auto firstP2 = std::make_pair(A, graphA);
        auto secondP2 = std::make_pair(C, graphB);

        double resd = boostSCD.Compute_Isometry(firstP2,secondP2, matchingtwo, totalTimeForSMIstd, totalTimeForBD);

        std::vector<int> matchingthree;

        auto firstP3 = std::make_pair(A, graphA);
        auto secondP3 = std::make_pair(B, graphB);

        double resd2 = nestedSCD.Compute_Isometry(firstP3,secondP3, matchingthree, totalTimeForSMIstd, totalTimeForBD);

        std::vector<int> matchingfour;

        auto firstP4 = std::make_pair(A, graphA);
        auto secondP4 = std::make_pair(B, graphB);

        double resd4 = boostSCD.Compute_Isometry(firstP4,secondP4, matchingfour, totalTimeForSMIstd, totalTimeForBD);

//        std::vector<int> matchingthree;
//
//        auto firstP3 = std::make_pair(A, graphA);
//        auto secondP3 = std::make_pair(C, graphB);
//
//        double ruru = boostSCD.Compute_Isometry(firstP3,secondP3, matchingthree, totalTimeForSMIstd, totalTimeForBD);


        std::vector<int> matchingOne;
        double purebd = eucBd.compute(A,BB, matchingOne);

        //res = nestedNaive.Compute_Isometry(A,B, matching, totalTimeForSMIstd, totalTimeForBD);
        //  std::cout << "Starting the program" << std::endl;
        //   std::cout << "A cols: " << A.n_cols << " rows: " << A.n_rows << std::endl;
        //   std::cout << "B cols: " << B.n_cols << " rows: " << B.n_rows << std::endl;
        //  std::cout << "B cols: " << B.n_cols << " rows: " << B.n_rows <<< std::endl;

        //   std::vector<int> matching2;
        //  auto firstP2 = std::make_pair(A, graphA);
        //  auto secondP2 = std::make_pair(C, graphB);

        //  double plainDist = nestedSCD.Compute_Isometry(firstP2,secondP2, matching2, totalTimeForSMIstd, totalTimeForBD);


//        std::vector<int> matching3;
//        auto firstP3 = std::make_pair(B, graphA);
//        auto secondP3 = std::make_pair(BB, graphB);
//
//        double stupidPairCompare = nestedSCD.Compute_Isometry(firstP3,secondP3, matching3, totalTimeForSMIstd, totalTimeForBD);
//
//        std::cout << "Stupid pair compare: " << stupidPairCompare << std::endl;
//
//         for (int i = 0; i < matching3.size(); i++)
//        {
//
//            std::cout << "Matched : " << i << " --- " << matching3[i] << std::endl;
//
//
//        }




        //  std::cout << "Program finnished" << std::endl;
        ratiod = res / purebd;

        //!  res =
        //!  res = CrystData::computeDistanceBetweenTwoClouds(baseCloud,perturbedCloud,matchingA,matchingB,simpleStringVector);

        //! Debug INformation:

        // std::cout << "Plain dist is: " << plainDist << std::endl;

//       std::cout << " Official results into: " << res << " Bd dist: " << purebd << " Ratio: " << ratiod << std::endl;
//
//       std::cout << "Original matching was:" << std::endl;
//

        //std::cout << A << std::endl;
        //   std::cout << C << std::endl;
        //  std::cout << "Curepsilon: " << curEpsion << std::endl;

        correctMatching = true;
        realDist = -9000;

        std::cout << A << std::endl;
        std::cout << C << std::endl;

        for (int i = 0; i < justPermute.size(); i++)
        {
            std::cout << "Comparing" << std::endl;
            std::cout << A.col(i) << std::endl;
            std::cout << C.col(matching[i]) << std::endl;
            double distP = arma::norm(A.col(i)-C.col(matching[i]), "inf");

            std::cout << i << " --- " << matching[i] << "|  should be: " << justPermute[i] << " | matchingtwo is: " << matchingtwo[i] << " : " << distP << std::endl;
            if (distP > curEpsion)
            {
                correctMatching = false;

            }
            realDist = std::max(realDist, distP);

            //  std::cout << "Matched : " << i << " --- " << justPermute[i] << std::endl;

        }
//
//        std::cout << "Following matching was obtained: " << std::endl;
//
//        for (int i = 0; i < matching.size(); i++)
//        {
//
//            std::cout << "Matched : " << i << " --- " << matching[i] << std::endl;
//
//
//        }

        std::cout << "Real matching should be:" << std::endl;
        ExtraExperiments::printVector(justPermute);
        std::cout << "Boost obtained (with rotation): " << std::endl;
        ExtraExperiments::printVector(matching);
        std::cout << "Boost obtained: (without rotation)" << std::endl;
        ExtraExperiments::printVector(matchingtwo);
        std::cout << "No permutation original: " << resd2 << std::endl;
        ExtraExperiments::printVector(matchingthree);
        std::cout << "No permutation boost:" << resd4 << std::endl;
        ExtraExperiments::printVector(matchingfour);

        std::cout << "Original clouds : " << std::endl;

        std::cout << A << std::endl;

        std::cout << B << std::endl;

        std::vector<std::vector<double>> distMat;

        ExtraExperiments::computeDistanceMatrix(A, B, distMat);
        ExtraExperiments::printMatrix(distMat);


        std::string folder = "/home/yury/LocalProjects/dbg/Matrices/";
        std::string saveFileA = "A" + std::to_string(numEpsilon) +  ".mat";
        std::string saveFileB = "B" + std::to_string(numEpsilon) + ".mat";

        A.save(folder + saveFileA, arma::arma_ascii);
        B.save(folder + saveFileB, arma::arma_ascii);


        //! ExtraExperiments::printVector(matchingthree);
        std::cout << "Resd is (no rotation) : " << resd << std::endl;
        std::cout << "Epsilon: " << numEpsilon << " repeat: " << repeat << " ended" << std::endl;

    }

    static void printMatrix(std::vector<std::vector<double>> & values, std::vector<double> & xCoord, std::vector<double> & yCoord, std::string output)
    {
        std::cout << "xcoordsize" << xCoord.size() << std::endl;
        std::cout << "ycoordsize" << yCoord.size() << std::endl;
        std::ofstream out(output);
        out << "0,";
        for (int i = 0; i < xCoord.size(); i++)
        {
            out << xCoord[i];
            if (i < xCoord.size()-1)
            {
                out << ",";
            }

        }
        out << std::endl;
        for (int i = 0; i < values.size(); i++)
        {
            out << yCoord[i] << ",";

            for (int j = 0; j < values[i].size(); j++)
            {

                out << values[i][j];
                if (j < values[i].size()-1)
                {
                    out << ",";
                }
            }
            out << std::endl;
        }



    }

    static void CgenerateCloudsArray(std::vector<double> & inbuds,std::vector<std::vector<double>> & cloud,std::vector<std::vector<double>> & cloudPair)
    {
        double a = inbuds[0];
        double b = inbuds[1];
        double c = inbuds[2];
        double d = inbuds[3];
        cloud.push_back(std::vector<double> {4*a,0});
        cloud.push_back(std::vector<double> {b,c});
        cloud.push_back(std::vector<double> {-b,-c});
        cloud.push_back(std::vector<double> {0,4*d});

        cloudPair.push_back(std::vector<double> {4*a,0});
        cloudPair.push_back(std::vector<double> {b,c});
        cloudPair.push_back(std::vector<double> {-b,-c});
        cloudPair.push_back(std::vector<double> {0,-4*d});

    }

    static void PrintExperiments(std::vector<double> & finalResult, std::vector<double> & maxx, std::vector<double> & minn, std::vector<double> & avg, std::vector<double> & epsilons, std::string outputpath )
    {
        std::ofstream out(outputpath);
        out << "e,result,maxRatio,minRatio,avgRatio" << std::endl;
        for (int i = 0; i < finalResult.size(); i++)
        {
            out << epsilons[i] << ", " << finalResult[i] << ", " << maxx[i] << ", " << minn[i] << ", " << avg[i] << std::endl;
        }
        out.close();


    }

      static void PrintExperimentsModi(std::vector<double> & maxx, std::vector<double> & minn, std::vector<double> & avg, std::vector<double> & epsilons, std::string outputpath )
    {
        std::ofstream out(outputpath);
        out << "e,maxRatio,minRatio,avgRatio" << std::endl;
        for (int i = 0; i < maxx.size(); i++)
        {
            out << epsilons[i] << ", " << maxx[i] << ", " << minn[i] << ", " << avg[i] << std::endl;
        }
        out.close();


    }

    static double shiftSingle(std::vector<double> & vv, double epsilon)
    {
        double dist = std::numeric_limits<double>::max();
        std::vector<double> perturbVector(vv.size());
        while(dist > epsilon)
        {
            perturbVector = std::vector<double>(vv.size());
            for (int i = 0; i < vv.size(); i++)
            {
                double epp = -epsilon + 2*((double) rand() / (RAND_MAX)) * epsilon;
                perturbVector[i] = perturbVector[i] + epp;
            }
            dist = sqrt(vectorTools::SquaredNorm(perturbVector));
        }
        vectorTools::vectorAddition(vv,perturbVector);
        return dist;
    }

    static double perturbCloud(std::vector<std::vector<double>> & cloudIN, std::vector<std::vector<double>> & cloudOUT, double epsilon)
    {
        double dd = std::numeric_limits<double>::min();
        for (int i  = 0; i < cloudIN.size(); i++)
        {
            std::vector<double> tmpElem = cloudIN[i];
            double dum = shiftSingle(tmpElem, epsilon);
            dd = std::max(dd, dum);
            cloudOUT.push_back(tmpElem);
        }
        return dd;
    }

    static void PerturbExperimentsSingle(int dim, int numPoints,  double curEpsion,
                                         int repeat, int numEpsilon, double & res, double & ratiod)
    {

        std::cout << "Epsilon: " << numEpsilon << " repeat: " << repeat << " started" << std::endl;

        typedef BottleneckCalculator<Bipartie_Graph,GraphTreeController> bipartiebd;
        typedef bd_general_interface<bipartiebd, bipartiebd> naiveBdCalc;

        EuclideanBDCompute<bipartiebd> eucBd;
        //!  Nested_Base<SCD, naiveBdCalc, std::pair<arma::mat, std::vector<std::vector<bool>>>> nestedNaive;
        Nested_Base<PBI, naiveBdCalc, arma::mat> nestedNaive;
        Nested_Abstract<SCD_Updated, std::tuple<double,arma::mat, arma::mat>, naiveBdCalc, std::pair<arma::mat, std::vector<std::vector<bool>>>> nestedSCD;

        std::vector<std::vector<double>> baseCloud;
        ExtraExperiments::generateRandomCloudInBoxNd(baseCloud, numPoints, 1.0, dim);


        std::vector<std::string> simpleStringVector;
        std::vector<std::vector<double>> perturbedCloud;
        ExtraExperiments::perturbCloud(baseCloud,perturbedCloud, curEpsion);
        //! auto rng = std::default_random_engine {};
        //!std::shuffle(std::begin(perturbedCloud), std::end(perturbedCloud), rng);
        arma::mat A;
        arma::mat B;

        long totalTimeForSMIstd;
        long totalTimeForBD;
        std::vector<int> matching;
        vectorTools::vec_to_mat(baseCloud, A);
        vectorTools::vec_to_mat(perturbedCloud, B);

        std::vector<std::vector<bool>> graphA(A.n_cols+1, std::vector<bool>(A.n_cols+1, true));

        std::vector<std::vector<bool>> graphB(B.n_cols+1, std::vector<bool>(B.n_cols+1, true));

        ExtraExperiments::randomGraphGenerator(graphA);
        graphB = graphA;
//!        ExtraExperiments::randomGraphGenerator(graphB);

        auto firstP = std::make_pair(A, graphA);
        auto secondP = std::make_pair(B, graphB);

        std::vector<int> matchingOne;
        double purebd = eucBd.compute(A,B, matchingOne);

        //res = nestedNaive.Compute_Isometry(A,B, matching, totalTimeForSMIstd, totalTimeForBD);
        //  std::cout << "Starting the program" << std::endl;
        //   std::cout << "A cols: " << A.n_cols << " rows: " << A.n_rows << std::endl;
        //   std::cout << "B cols: " << B.n_cols << " rows: " << B.n_rows << std::endl;
        //  std::cout << "B cols: " << B.n_cols << " rows: " << B.n_rows <<< std::endl;
        res = nestedSCD.Compute_Isometry(firstP,secondP, matching, totalTimeForSMIstd, totalTimeForBD, 2);
        //  std::cout << "Program finnished" << std::endl;
        ratiod = res / purebd;

        //!  res =
        //!  res = CrystData::computeDistanceBetweenTwoClouds(baseCloud,perturbedCloud,matchingA,matchingB,simpleStringVector);
        std::cout << "Epsilon: " << numEpsilon << " repeat: " << repeat << " ended" << std::endl;
    }

    static void generateRandomCloudInBoxNd(std::vector<std::vector<double>> & cloud, int n, double boxSize, int dimm)
    {
        for (int j = 0; j < n; j++)
        {
            std::vector<double> elem;
            for (int i = 0; i < dimm; i++)
            {
                double tmpValue = ((double) rand() / (RAND_MAX)) * boxSize;
                elem.push_back(tmpValue);
            }
            cloud.push_back(elem);
        }
    }

    static void LinkageExperiments(int threads, int dim, int numPoints, double maxEpsilon, int numEpsilon, int repeat, std::string outputPath)
    {

        BS::thread_pool pool(threads);
        std::vector<double> epsilons;

        std::vector<std::string> simpleStringVector;



        double epsilonstep = maxEpsilon / (double) (numEpsilon);
        std::vector<std::vector<double>> resultsMassive(numEpsilon, std::vector<double>(repeat));
        std::vector<std::vector<double>> ratioMassive(numEpsilon, std::vector<double>(repeat));
        std::vector<std::thread> vec_thr;

        for (int i = 0; i < numEpsilon; i++)
        {
            double curEpsilon = epsilonstep * (i+1);
            epsilons.push_back(curEpsilon);
            for (int j = 0; j < repeat; j++)
            {
                //    std::cout << "Combo: " << i << ", " << j << std::endl;

                pool.push_task(ExtraExperiments::PerturbExperimentsSingle,dim, numPoints, curEpsilon, j, i, std::ref(resultsMassive[i][j]), std::ref(ratioMassive[i][j]));

            }

        }

        pool.wait_for_tasks();
        // for(auto& t: vec_thr)
        //    t.join();

        std::vector<double> finalResults;
        std::vector<double> ratioMax;
        std::vector<double> ratioMin;
        std::vector<double> ratioAvg;

        for (int i = 0; i < numEpsilon; i++)
        {

            double sumdistance = 0;
            double maxxRatio = -9999999;
            double minnRatio = 99999999;
            double avgRatio = 0;
            for (int j = 0; j < repeat; j++)
            {
                sumdistance = sumdistance + resultsMassive[i][j];
                avgRatio = avgRatio + ratioMassive[i][j];
                maxxRatio = std::max(maxxRatio, ratioMassive[i][j]);
                minnRatio = std::min(minnRatio, ratioMassive[i][j]);
            }
            sumdistance = sumdistance / (double) repeat;
            avgRatio = avgRatio /  (double) repeat;
            finalResults.push_back(sumdistance);
            ratioMax.push_back(maxxRatio);
            ratioMin.push_back(minnRatio);
            ratioAvg.push_back(avgRatio);
        }

        ExtraExperiments::PrintExperiments(finalResults, ratioMax, ratioMin, ratioAvg, epsilons, outputPath);
    }


    static void PerturbationExperiments(int threads, int dim, int numPoints, double maxEpsilon, int numEpsilon, int repeat, std::string outputPath)
    {

        BS::thread_pool pool(threads);
        std::vector<double> epsilons;

        std::vector<std::string> simpleStringVector;



        double epsilonstep = maxEpsilon / (double) (numEpsilon);
        std::vector<std::vector<double>> resultsMassive(numEpsilon, std::vector<double>(repeat));
        std::vector<std::vector<double>> ratioMassive(numEpsilon, std::vector<double>(repeat));
        std::vector<std::thread> vec_thr;

        for (int i = 0; i < numEpsilon; i++)
        {
            double curEpsilon = epsilonstep * (i+1);
            epsilons.push_back(curEpsilon);
            for (int j = 0; j < repeat; j++)
            {
                //    std::cout << "Combo: " << i << ", " << j << std::endl;

                pool.push_task(ExtraExperiments::PerturbExperimentsSingle,dim, numPoints, curEpsilon, j, i, std::ref(resultsMassive[i][j]), std::ref(ratioMassive[i][j]));

            }

        }

        pool.wait_for_tasks();
        // for(auto& t: vec_thr)
        //    t.join();

        std::vector<double> finalResults;
        std::vector<double> ratioMax;
        std::vector<double> ratioMin;
        std::vector<double> ratioAvg;

        for (int i = 0; i < numEpsilon; i++)
        {

            double sumdistance = 0;
            double maxxRatio = -9999999;
            double minnRatio = 99999999;
            double avgRatio = 0;
            for (int j = 0; j < repeat; j++)
            {
                sumdistance = sumdistance + resultsMassive[i][j];
                avgRatio = avgRatio + ratioMassive[i][j];
                maxxRatio = std::max(maxxRatio, ratioMassive[i][j]);
                minnRatio = std::min(minnRatio, ratioMassive[i][j]);
            }
            sumdistance = sumdistance / (double) repeat;
            avgRatio = avgRatio /  (double) repeat;
            finalResults.push_back(sumdistance);
            ratioMax.push_back(maxxRatio);
            ratioMin.push_back(minnRatio);
            ratioAvg.push_back(avgRatio);
        }

        ExtraExperiments::PrintExperiments(finalResults, ratioMax, ratioMin, ratioAvg, epsilons, outputPath);
    }

    static void experimentsConducting(int n, std::string folder)
    {

        typedef BottleneckCalculator<Bipartie_Graph,GraphTreeController> bipartiebd;
        typedef bd_general_interface<bipartiebd, bipartiebd> naiveBdCalc;

        EuclideanBDCompute<bipartiebd> eucBd;
        //! Nested_Updated<PBI_Update, naiveBdCalc, arma::mat> nestedNaive;

        Nested_Updated<PBI_Update<naiveBdCalc>, arma::mat> nestedNaive;

        //!   bd = nestedNaive.Compute_Isometry(matrices[u], matrices[v], matching, totalTimeForSMIstd, totalTimeForBD);

        std::cout << "Method started " << std::endl;


        std::vector<int> mappingA = std::vector<int> {0,1,2,3};
        std::vector<int> mappingB = std::vector<int> {0,1,2,3};
        std::vector<int> one = std::vector<int> {0,1};
        std::vector<int> two = std::vector<int> {0,2};
        std::vector<int> three = std::vector<int> {0,3};
        std::vector<int> four = std::vector<int> {1,2};
        std::vector<int> five = std::vector<int> {1,3};
        std::vector<int> six = std::vector<int> {2,3};
        std::vector<std::vector<int>> inclusion;
        std::vector<std::vector<int>> exclusion;
        inclusion.push_back(one);
        inclusion.push_back(two);
        inclusion.push_back(three);
        inclusion.push_back(four);
        inclusion.push_back(five);
        inclusion.push_back(six);
        exclusion.push_back(six);
        exclusion.push_back(five);
        exclusion.push_back(four);
        exclusion.push_back(three);
        exclusion.push_back(two);
        exclusion.push_back(one);
        double divider = (double) n;
        double stepsizeA = 1.0 / divider;
        double stepsizeB = 4.0 / divider;
        double stepsizeC = 4.0 / divider;
        double stepsizeD = 1.0 / divider;
        std::vector<double> dividers = std::vector<double> {stepsizeA,stepsizeB,stepsizeC,stepsizeD};
        std::vector<double> constants = std::vector<double> {0.5,2.0,2.0,0.5};
        std::vector<std::string> names{"a","b","c","d"};

        std::cout << "Finnished preparade " << std::endl;

        for (int i = 0; i < inclusion.size(); i++)
        {
            std::vector<double> valuesX;
            std::vector<double> valuesY;
            std::vector<std::vector<double>> matrix(n, std::vector<double>(n));
            std::string fileoutput = folder + "(" + names[inclusion[i][0]] + ", " + names[inclusion[i][1]] + ")" + "(" + names[exclusion[i][0]] + ":" + std::to_string(constants[exclusion[i][0]]) + ")"
                                     + "(" + names[exclusion[i][1]] + ":" + std::to_string(constants[exclusion[i][1]]) + ")"   + ".csv";
            for (int j = 1; j <= n; j++)
            {
                valuesX.push_back(j*dividers[inclusion[i][0]]);
                valuesY.push_back(j*dividers[inclusion[i][1]]);
                for (int k = 1; k <= n; k++)
                {
                    std::cout << "Begginig " << i << " , " << j << ", " << k << std::endl;
                    std::vector<std::vector<double>> twoDcloudA;
                    std::vector<std::vector<double>> twoDcloudB;
                    std::vector<double> inpudd(4);
                    inpudd[inclusion[i][0]] = j*dividers[inclusion[i][0]];
                    inpudd[inclusion[i][1]] = k*dividers[inclusion[i][1]];
                    inpudd[exclusion[i][0]] = constants[exclusion[i][0]];
                    inpudd[exclusion[i][1]] = constants[exclusion[i][1]];
                    ExtraExperiments::CgenerateCloudsArray(inpudd,twoDcloudA,twoDcloudB);
                    //   std::cout << "Cloud generated" << std::endl;
                    //    CrystData::PrintVector<double>(inpudd);
                    // //   CrystData::PrintMatrix(twoDcloudA);
                    //    CrystData::PrintMatrix(twoDcloudB);
                    std::vector<std::string> outputpairs;
                    //! double answer = CrystData::computeDistanceBetweenTwoClouds(twoDcloudA,twoDcloudB,mappingA,mappingB,outputpairs);

                    arma::mat A;
                    arma::mat B;
                    long totalTimeForSMIstd;
                    long totalTimeForBD;
                    std::vector<int> matching;
                    vectorTools::vec_to_mat(twoDcloudA, A);
                    vectorTools::vec_to_mat(twoDcloudB, B);


                    double bd = nestedNaive.Compute_Isometry(A,B, matching, totalTimeForSMIstd, totalTimeForBD);

                    matrix[j-1][k-1] = bd;
                }
            }
            ExtraExperiments::printMatrix(matrix,valuesX,valuesY,fileoutput);
        }



    }

    static void TimeExperimentsSingle(int dim, int numPoints,  double curEpsion,
                                      int repeat, int iter, long & totalTimeForSMIstd, long & totalTimeForBD)
    {

        std::cout << "Iteration: " << iter << " repeat: " << repeat << " started" << std::endl;

        typedef BottleneckCalculator<Bipartie_Graph,GraphTreeController> bipartiebd;
        typedef bd_general_interface<bipartiebd, bipartiebd> naiveBdCalc;


        EuclideanBDCompute<bipartiebd> eucBd;
        //!  Nested_Base<SCD, naiveBdCalc, std::pair<arma::mat, std::vector<std::vector<bool>>>> nestedNaive;
        Nested_Base<PBI, naiveBdCalc, arma::mat> nestedNaive;
        Nested_Abstract<SCD_Updated, std::tuple<double,arma::mat, arma::mat>, naiveBdCalc, std::pair<arma::mat, std::vector<std::vector<bool>>>> nestedSCD;


        std::vector<std::vector<double>> baseCloud;
        ExtraExperiments::generateRandomCloudInBoxNd(baseCloud, numPoints, 1.0, dim);


        std::vector<std::string> simpleStringVector;
        std::vector<std::vector<double>> perturbedCloud;
        ExtraExperiments::perturbCloud(baseCloud,perturbedCloud, curEpsion);
        // auto rng = std::default_random_engine {};
        //std::shuffle(std::begin(perturbedCloud), std::end(perturbedCloud), rng);
        arma::mat A;
        arma::mat B;
        // long totalTimeForSMIstd;
        // long totalTimeForBD;
        std::vector<int> matching;
        vectorTools::vec_to_mat(baseCloud, A);
        vectorTools::vec_to_mat(perturbedCloud, B);
        std::vector<int> matchingOne;

        std::vector<std::vector<bool>> graphA(A.n_cols+1, std::vector<bool>(A.n_cols+1, true));
        std::vector<std::vector<bool>> graphB(B.n_cols+1, std::vector<bool>(B.n_cols+1, true));

        ExtraExperiments::randomGraphGenerator(graphA);
        graphB = graphA;
        //  ExtraExperiments::randomGraphGenerator(graphB);

        auto firstP = std::make_pair(A, graphA);
        auto secondP = std::make_pair(B, graphB);

        double res = nestedSCD.Compute_Isometry(firstP,secondP, matching, totalTimeForSMIstd, totalTimeForBD, 2);


        // double purebd = eucBd.compute(A,B, matchingOne);

        //! double res = nestedNaive.Compute_Isometry(A,B, matching, totalTimeForSMIstd, totalTimeForBD);
        //  ratiod = res / purebd;

        //!  res =
        //!  res = CrystData::computeDistanceBetweenTwoClouds(baseCloud,perturbedCloud,matchingA,matchingB,simpleStringVector);
        std::cout << "Iteration: " << iter << " repeat: " << repeat << " ended" << std::endl;
    }




    static void TimeExperiments(int threads, int dim, int steps, int stepSize, int repeat, double curEpsilon, std::string outputPath)
    {

        BS::thread_pool pool(threads);

        std::vector<std::vector<long>> timeInvariant(steps, std::vector<long>(repeat));
        std::vector<std::vector<long>> timeBd(steps, std::vector<long>(repeat));


        std::vector<int> numPointsVector;

        for (int i = 0; i < steps; i++)
        {
            numPointsVector.push_back((i+1)*stepSize);
        }


        for (int i = 0; i < steps; i++)
        {


            for (int j = 0; j < repeat; j++)
            {

                //! TimeExperimentsSingle(int dim, int numPoints,  double curEpsion,int repeat, int iter, long & totalTimeForSMIstd, long & totalTimeForBD)
                pool.push_task(ExtraExperiments::TimeExperimentsSingle,dim, (i+1)*stepSize, curEpsilon, j, i, std::ref(timeInvariant[i][j]), std::ref(timeBd[i][j]));



            }


        }
        pool.wait_for_tasks();


//             avgSMI = (long) avgSMI / repeat;
//            avgSBM = (long) avgSBM / repeat;
//            tSMI.push_back(avgSMI);
//            tSBM.push_back(avgSBM);

        std::vector<double> tSMI;
        std::vector<double> tSBM;
        std::vector<double> tSum;

        for (int i = 0; i < steps; i++)
        {

            double sumSmi = 0;
            double sumTSBM = 0;
            for (int j = 0; j < repeat; j++)
            {
                sumSmi = sumSmi + timeInvariant[i][j];
                sumTSBM = sumTSBM + timeBd[i][j];

            }
            sumSmi = sumSmi / (double) repeat;
            sumTSBM = sumTSBM / (double) repeat;
            tSMI.push_back(sumSmi );
            tSBM.push_back(sumTSBM);
            double sumdx = sumSmi + sumTSBM;
            tSum.push_back(sumdx);
            //   sumdistance = sumdistance / (double) repeat;
            //    avgRatio = avgRatio /  (double) repeat;
            //    finalResults.push_back(sumdistance);
            //    ratioMax.push_back(maxxRatio);
            //    ratioMin.push_back(minnRatio);
            //   ratioAvg.push_back(avgRatio);
        }


        ExtraExperiments::PrintExperimentsTime(numPointsVector, tSMI, tSBM, tSum, outputPath);
    }

    static void PrintExperimentsTime(std::vector<int> & finalResult, std::vector<double> & tSMI, std::vector<double> & tSBM, std::vector<double> & tSum, std::string outputpath )
    {
        std::ofstream out(outputpath);
        out << "p, tSMI, tSBM, tSum " << std::endl;
        for (int i = 0; i < finalResult.size(); i++)
        {
            out << finalResult[i] << ", " << tSMI[i] << ", " << tSBM[i] << ", " << tSum[i] << std::endl;
        }
        out.close();


    }

    static void MatchingExperiments(int threads, int dim, int steps, int stepSize, int repeat, double curEpsilon, std::string outputPath)
    {

        BS::thread_pool pool(threads);

        std::vector<std::vector<long>> timeInvariant(steps, std::vector<long>(repeat));
        std::vector<std::vector<long>> timeBd(steps, std::vector<long>(repeat));

        std::vector<std::vector<std::vector<double>>> DataVec(6,std::vector<std::vector<double>>(steps, std::vector<double>(repeat)));
        //   std::vector<std::vector<double>> resratiovec(steps, std::vector<double>(repeat));
        //    std::vector<std::vector<double>> distmatchvec(steps, std::vector<double>(repeat));
        //   std::vector<std::vector<double>> ratiodistbdvec(steps, std::vector<double>(repeat));


        std::vector<int> numPointsVector;

        for (int i = 0; i < steps; i++)
        {
            numPointsVector.push_back((i+1)*stepSize);
        }


        for (int i = 0; i < steps; i++)
        {


            for (int j = 0; j < repeat; j++)
            {

                //! TimeExperimentsSingle(int dim, int numPoints,  double curEpsion,int repeat, int iter, long & totalTimeForSMIstd, long & totalTimeForBD)
//                pool.push_task(ExtraExperiments::PerturbSingleNCDSmartComplicatedZ,dim, (i+1)*stepSize, curEpsilon, j, i,
//                               std::ref(DataVec[0][i][j]),  std::ref(DataVec[1][i][j]),  std::ref(DataVec[2][i][j]),  std::ref(DataVec[3][i][j]),
//                               std::ref(timeInvariant[i][j]), std::ref(timeBd[i][j]));

                pool.push_task(ExtraExperiments::PerturbSingleNCDSmartComplicatedZ,dim, (i+1)*stepSize, curEpsilon, j, i,
                               std::ref(DataVec[0][i][j]),  std::ref(DataVec[1][i][j]),  std::ref(DataVec[2][i][j]),  std::ref(DataVec[3][i][j]),
                               std::ref(DataVec[4][i][j]), std::ref(DataVec[5][i][j]));

            }


        }
        pool.wait_for_tasks();
        std::cout << "Experiments finnisehd " << std::endl;


//             avgSMI = (long) avgSMI / repeat;
//            avgSBM = (long) avgSBM / repeat;
//            tSMI.push_back(avgSMI);
//            tSBM.push_back(avgSBM);

        std::vector<double> tSMI;
        std::vector<double> tSBM;
        std::vector<double> tSum;
        std::vector<std::vector<double>> extraData(6, std::vector<double>());
        std::vector<std::vector<double>> extraDataMax(6, std::vector<double>());
        std::vector<std::vector<double>> extraDataMin(6, std::vector<double>());
        //  std::vector<double> resAvg;
        //  std::vector<double> rasratioAvg;
        //  std::vector<double> distmatchvecAvg;
        //  std::vector<double> ratiodistbdvecAvg;

        for (int i = 0; i < steps; i++)
        {

            double sumSmi = 0;
            double sumTSBM = 0;

            std::vector<double> sumData(6,0);
            std::vector<double> maxData(6,std::numeric_limits<double>::min());
            std::vector<double> minData(6,std::numeric_limits<double>::max());

            //  double sumRes;
            // double sumResRatio;
            //  double sumdistmatch;
            //  double sumratiodistbd;

            for (int j = 0; j < repeat; j++)
            {
                sumSmi = sumSmi + timeInvariant[i][j];
                sumTSBM = sumTSBM + timeBd[i][j];
                for (int k = 0; k < 6; k++)
                {
                    maxData[k] = std::max(maxData[k], DataVec[k][i][j]);
                    minData[k] = std::min(minData[k], DataVec[k][i][j]);
                    sumData[k] = sumData[k] + DataVec[k][i][j];

                }


            }
            for (int k = 0; k < 6; k++)
            {

                sumData[k] = sumData[k] / (double) repeat;

            }
//            for (int k = 0; k < 6; k++)
//            {
//
//             //   extraData[k].push_back(sumData[k]);
//                extraDataMax[k].push_back(maxData[k]);
//                extraDataMin[k].push_back(minData[k]);
//
//            }
            //   sumdistance = sumdistance / (double) repeat;
            //    avgRatio = avgRatio /  (double) repeat;
            //    finalResults.push_back(sumdistance);
            //    ratioMax.push_back(maxxRatio);
            //    ratioMin.push_back(minnRatio);
            //   ratioAvg.push_back(avgRatio);
        }


        //ExtraExperiments::PrintExperimentsMatchTwo(extraDataMax, extraDataMin, extraData, numPointsVector, tSMI, tSBM, tSum, outputPath);
        ExtraExperiments::PrintExperimentsMatchTwo(extraDataMax, extraDataMin, extraData, numPointsVector, outputPath);
    }

    static void MatchingExperiments2D(int threads, int dim, int stepsPts, int stepSizePts, int stepsEps, double stepSizeEps, int repeat, std::string outputPath)
    {

        BS::thread_pool pool(threads);


        std::vector<std::vector<std::vector<std::vector<double>>>> DataVec(6,std::vector<std::vector<std::vector<double>>>(stepsEps, std::vector<std::vector<double>>(stepsPts, std::vector<double>(repeat))));
        //   std::vector<std::vector<double>> resratiovec(steps, std::vector<double>(repeat));
        //    std::vector<std::vector<double>> distmatchvec(steps, std::vector<double>(repeat));
        //   std::vector<std::vector<double>> ratiodistbdvec(steps, std::vector<double>(repeat));


        std::vector<int> numPointsVector;
        std::vector<double> ptsVec;

        for (int i = 0; i < stepsPts; i++)
        {
            numPointsVector.push_back((i+1)*stepSizePts);
            double lel = (double) (i+1)*stepSizePts;
            ptsVec.push_back(lel);
        }

        std::vector<double> epsPointsVector;
        for (int i = 0; i < stepsEps; i++)
        {
            epsPointsVector.push_back((i+1)*stepSizeEps);


        }

        std::cout << "Sizes: " << DataVec.size() << " , " << DataVec[0].size() << ", " << DataVec[0][0].size() << " , " << DataVec[0][0][0].size() << std::endl;
        std::cout << "epsPointsVector: " << epsPointsVector.size() << std::endl;
        std::cout << "Ptspointvector: " << numPointsVector.size() << std::endl;

        for (int k = 0; k < stepsEps; k++)
        {

            for (int i = 0; i < stepsPts; i++)
            {


                for (int j = 0; j < repeat; j++)
                {

                    //! TimeExperimentsSingle(int dim, int numPoints,  double curEpsion,int repeat, int iter, long & totalTimeForSMIstd, long & totalTimeForBD)
                    pool.push_task(ExtraExperiments::PerturbSingleNCDSmartComplicatedZ,dim,  numPointsVector[i], epsPointsVector[k], j, k,
                                   std::ref(DataVec[0][k][i][j]),  std::ref(DataVec[1][k][i][j]),  std::ref(DataVec[2][k][i][j]),  std::ref(DataVec[3][k][i][j]),
                                   std::ref(DataVec[4][k][i][j]),std::ref(DataVec[5][k][i][j]));



                }


            }
        }

        std::cout << "Finnished computing everything " << std::endl;
        pool.wait_for_tasks();
        std::cout << "Finnished computing everything x2 " << std::endl;


//             avgSMI = (long) avgSMI / repeat;
//            avgSBM = (long) avgSBM / repeat;
//            tSMI.push_back(avgSMI);
//            tSBM.push_back(avgSBM);

        // std::vector<double> tSMI;
        // std::vector<double> tSBM;
        //  std::vector<double> tSum;
        std::vector<std::vector<std::vector<std::vector<double>>>> extraData(3, std::vector<std::vector<std::vector<double>>>(7, std::vector<std::vector<double>>(stepsEps,std::vector<double>(stepsPts,0))));
        //  std::vector<std::vector<std::vector<double>>>  extraDataMax(7, std::vector<std::vector<double>>(stepsEps,std::vector<double>(stepsPts,0)));
        //   std::vector<std::vector<std::vector<double>>>  extraDataMin(7, std::vector<std::vector<double>>(stepsEps,std::vector<double>(stepsPts,0)));
        //  std::vector<double> resAvg;
        //  std::vector<double> rasratioAvg;
        //  std::vector<double> distmatchvecAvg;
        //  std::vector<double> ratiodistbdvecAvg;

        std::vector<std::string> attributeProperty{"max","min","avg"};

        std::vector<std::string> attributeNames{"NBM(cloud,pertrubation)", "NBM / pertubation_error", "Euclid_match_dist", "Euclid_match_dist / NBM", "NCD_Time_ms","NBM_Time_ms","total_time_ms"};

        for (int h = 0; h < stepsEps; h++)
        {

            for (int i = 0; i < stepsPts; i++)
            {

                double sumSmi = 0;
                double sumTSBM = 0;

                std::vector<double> sumData(7,0);
                std::vector<double> maxData(7,std::numeric_limits<double>::min());
                std::vector<double> minData(7,std::numeric_limits<double>::max());

                //  double sumRes;
                // double sumResRatio;
                //  double sumdistmatch;
                //  double sumratiodistbd;

                for (int j = 0; j < repeat; j++)
                {
                    for (int k = 0; k < 6; k++)
                    {
                        maxData[k] = std::max(maxData[k], DataVec[k][h][i][j]);
                        minData[k] = std::min(minData[k], DataVec[k][h][i][j]);
                        sumData[k] = sumData[k] + DataVec[k][h][i][j];

                    }
                    maxData[6] = std::max(maxData[6], DataVec[4][h][i][j] + DataVec[5][h][i][j]);
                    minData[6] = std::min(minData[6],DataVec[4][h][i][j] + DataVec[5][h][i][j]);
                    sumData[6] = sumData[6] + DataVec[4][h][i][j] + DataVec[5][h][i][j];



                }
                // sumSmi = sumSmi / (double) repeat;
                //   sumTSBM = sumTSBM / (double) repeat;
                for (int k = 0; k < 7; k++)
                {

                    sumData[k] = sumData[k] / (double) repeat;

                }
                for (int k = 0; k < 7; k++)
                {

                    extraData[0][k][h][i] = sumData[k];
                    extraData[1][k][h][i] =(maxData[k]);
                    extraData[2][k][h][i] =(minData[k]);

                }
                //   sumdistance = sumdistance / (double) repeat;
                //    avgRatio = avgRatio /  (double) repeat;
                //    finalResults.push_back(sumdistance);
                //    ratioMax.push_back(maxxRatio);
                //    ratioMin.push_back(minnRatio);
                //   ratioAvg.push_back(avgRatio);
            }
        }

        for (int pp = 0; pp < 3; pp++)
        {
            for (int qq = 0; qq < 7; qq++)
            {
                std::string ffpath = outputPath + attributeProperty[pp] + "_" + attributeNames[qq] + ".csv";
                ExtraExperiments::printMatrix(extraData[pp][qq],epsPointsVector,ptsVec,ffpath);
                // ExtraExperiments::PrintExperimentsMatch(extraDataMax, extraDataMin, extraData, numPointsVector, tSMI, tSBM, tSum, outputPath)

            }
        }

        //! Print averages first:



        //!    ExtraExperiments::PrintExperimentsMatch(extraDataMax, extraDataMin, extraData, numPointsVector, tSMI, tSBM, tSum, outputPath);
    }

    static void PrintExperimentsMatch(std::vector<std::vector<double>> & extraDataMax, std::vector<std::vector<double>> & extraDataMin,
                                      std::vector<std::vector<double>> & extraData, std::vector<int> & finalResult,
                                      std::vector<double> & tSMI, std::vector<double> & tSBM, std::vector<double> & tSum, std::string outputpath )
    {
        std::ofstream out(outputpath);
        out << "p,tSMI,tSBM,tSum,";
        out << "Avg-NBM,Avg-NBM-Ratio,Avg-MatchDist,Avg-MatchDistRatio,";
        out <<"Min-NBM,Min-NBM-Ratio,Min-MatchDist,Min-MatchDistRatio,";
        out << "Max-NBM,Max-NBM-Ratio,Max-MatchDist,Max-MatchDistRatio";
        out << std::endl;
        for (int i = 0; i < finalResult.size(); i++)
        {
            out << finalResult[i] << "," << tSMI[i] << "," << tSBM[i] << "," << tSum[i];
            for (int k = 0; k < 4; k++)
            {
                out << "," << extraData[k][i];
            }
            for (int k = 0; k < 4; k++)
            {
                out << "," << extraDataMin[k][i];
            }
            for (int k = 0; k < 4; k++)
            {
                out << "," << extraDataMax[k][i];
            }
            out << std::endl;
        }
        out.close();


    }

    static void PrintExperimentsMatchTwo(std::vector<std::vector<double>> & extraDataMax, std::vector<std::vector<double>> & extraDataMin,
                                         std::vector<std::vector<double>> & extraData, std::vector<int> & finalResult,std::string outputpath )
    {
        std::ofstream out(outputpath);
        out << "p";
        out << "Avg-NBM,Avg-NBM-Ratio,Avg-MatchDist,Avg-MatchDistRatio,Avg-invariant,Avg-bd";
        out <<"Min-NBM,Min-NBM-Ratio,Min-MatchDist,Min-MatchDistRatio,Min-invariant,Min-bd";
        out << ",Max-NBM,Max-NBM-Ratio,Max-MatchDist,Max-MatchDistRatio,Max-invariant,Max-bd";
        out << std::endl;
        for (int i = 0; i < finalResult.size(); i++)
        {
            out << finalResult[i];
            //out << finalResult[i] << "," << tSMI[i] << "," << tSBM[i] << "," << tSum[i];
            for (int k = 0; k < 6; k++)
            {
                out << "," << extraData[k][i];
            }
            for (int k = 0; k < 6; k++)
            {
                out << "," << extraDataMin[k][i];
            }
            for (int k = 0; k < 6; k++)
            {
                out << "," << extraDataMax[k][i];
            }
            out << std::endl;
        }
        out.close();


    }

    static void findCentroid(arma::mat & A, arma::vec & p)
    {
        arma::vec t(A.n_rows);
        for (int i = 0; i < A.n_cols; i++)
        {
            t = t + A.col(i);
        }
        t = t / A.n_cols;
        p = t;

    }
    static void ShiftBy(arma::mat & A, arma::vec & p)
    {
        for (int i = 0; i < A.n_cols; i++)
        {
            arma::vec tt = A.col(i) + p;
            A.col(i) = tt;
        }


    }



    static void findRotation(arma::mat & A, arma::mat & B, arma::mat & R, arma::vec & t)
    {
        arma::vec centroid_A;
        arma::vec centroid_B;
        findCentroid(A, centroid_A);
        findCentroid(B, centroid_B);
        arma::vec neg_centroid_A = (-1)*centroid_A;
        arma::vec neg_centroid_B = (-1)*centroid_B;
        arma::mat Am = A;
        arma::mat Bm = B;
        ShiftBy(Am, neg_centroid_A);
        ShiftBy(Bm, neg_centroid_B);
        arma::mat H = Am * Bm.t();
        arma::mat U;
        arma::vec S;
        arma::mat Vt;
        arma::svd(U,S,Vt, H);
        R = Vt * U.t();
        double dett = arma::det(R);
        //   std::cout << "DEterminant found as : " << dett << std::endl;
        arma::mat UnityMatrix = arma::mat(A.n_rows,A.n_rows,arma::fill::eye);

        //  std::cout << "Unity matrix " << std::endl;

        //std::cout << UnityMatrix << std::endl;
        std::cout << "Determinant detected as " << dett << std::endl;
        if (dett < 0)
        {


            UnityMatrix(A.n_rows-1, A.n_rows-1) = -1;
            R = Vt * UnityMatrix * U.t();

            //! Vt.t()

        }
        t = -R * centroid_A + centroid_B;


    }

    static void rotationTest(int dim, int numPoints)
    {
        std::cout << "Starting rotation test " << std::endl;
        std::vector<std::vector<double>> baseCloud;
        ExtraExperiments::generateRandomCloudInBoxNd(baseCloud, numPoints, 1.0, dim);
        arma::mat A;
        arma::mat B;
        arma::mat rotMat;
        vectorTools::vec_to_mat(baseCloud, A);
        generateRandomRotationMatrix2D(rotMat);
        B = rotMat * A;
        arma::vec tt;
        generateRandomVector(tt, 2, 1.0);

        std::cout << "Shift vector tt: " << std::endl;
        std::cout << tt << std::endl;
        ShiftBy(B,tt);

        std::cout << "Cloud A:" << std::endl;
        std::cout << A << std::endl;
        std::cout << "Cloud B" << std::endl;
        std::cout << B << std::endl;
        std::cout << "rotMat original:" << std::endl;
        std::cout << rotMat << std::endl;


        arma::mat R;
        arma::vec t;
        findRotation( A, B,  R,  t);

        std::cout << "Matrix R: " << std::endl;
        std::cout << R << std::endl;
        std::cout << "vector T:" << std::endl;
        std::cout << t << std::endl;


    }



protected:

private:
};

#endif // EXTRAEXPERIMENTS_H
