#include <iostream>

#include "Bottleneckdistance/Algorithms/BoostBottleneck.h"
#include "Bottleneckdistance/Algorithms/Classic/BottleneckCalculator.h"
#include "Bottleneckdistance/Misc/Bipartie_Graph.h"
#include "Bottleneckdistance/Misc/GraphTreeController.h"
#include "Bottleneckdistance/Misc/rtree_graph.h"
#include "Bottleneckdistance/Misc/rtreeController.h"
#include "Bottleneckdistance/bd_general_interface.h"
#include "Bottleneckdistance/Misc/EuclideanBDCompute.h"
#include "CloudIsometry/Original_Base.h"
#include "CloudIsometry/Nested_Base.h"
#include "CloudIsometry/Invariants/SCD.h"
#include "CloudIsometry/Invariants/PBI.h"
#include "CloudIsometry/Invariants/PBI_Original.h"
#include "Experiments/ExperimentRunner.h"
#include "Experiments/ExtraExperiments.h"
#include "Experiments/GraphExperiments.h"



using namespace std;

static void Vec2Mat( std::vector<std::vector<double>> & in, arma::mat & out)
{
    out.resize(in[0].size(), in.size());
    for (int i = 0; i < in.size() ; i++)
    {
        out.col(i) = arma::vec(in[i]);
    }


}

static void PermuteAbit(arma::mat & in, double epsilon)
{
    for (int i = 0; i < in.n_cols; i++)
    {
        for (int j = 0; j < in.n_rows; j++)
        {
            double epp = -epsilon + 2*((double) rand() / (RAND_MAX)) * epsilon;
            in.col(i)(j) = in.col(i)(j) + epp;


        }
    }
}

static void PermuteAbit(std::vector<std::vector<double>> & in, double epsilon)
{
    for (int i = 0; i < in.size(); i++)
    {
        for (int j = 0; j < in[0].size(); j++)
        {
            double epp = -epsilon + 2*((double) rand() / (RAND_MAX)) * epsilon;
            in[i][j] = in[i][j] + epp;


        }
    }
}





static void AdequareTestingSCDVer2()
{
    //  LayeredSCD pp(0.43,1);

    typedef BottleneckCalculator<Bipartie_Graph,GraphTreeController> bipartiebd;
    typedef bd_general_interface<bipartiebd, bipartiebd> naiveBdCalc;
    typedef bd_general_interface<BoostBottleneck,BoostBottleneck> boostBdCalc;

    EuclideanBDCompute<bipartiebd> eucBd;
    BottleneckCalculator<Bipartie_Graph,GraphTreeController> naiveBd;


  //  Original_Base<SCD, naiveBdCalc , std::pair<arma::mat, std::vector<std::vector<bool>>>> origNaive;
  //  Original_Base<SCD, boostBdCalc , std::pair<arma::mat, std::vector<std::vector<bool>>>> origBoost;
  //!  Nested_Base<SCD, naiveBdCalc , std::pair<arma::mat, std::vector<std::vector<bool>>>> nestedNaive;

    Nested_Base<PBI, naiveBdCalc , arma::mat> nestedNaive;
    Nested_Base<PBI_Original, naiveBdCalc , std::vector<std::vector<double>>> nestedClassic;

  //  Nested_Base<SCD, boostBdCalc , std::pair<arma::mat, std::vector<std::vector<bool>>>> nestedBoost;

    std::vector<std::string> failures;

    //LayeredSCD zz(0.43,2);
    //LayeredSCD pp(0.43,1);
    //LayeredSCD qq(0.43,0);
    int jmax = 15;
    int maxxx = 20;
    for(int i = 5; i < maxxx; i++)
    {
        int failure = 0;

        std::vector<std::vector<bool>> graphA(i+1, std::vector<bool>(i+1, true));
        std::vector<std::vector<bool>> graphB = graphA;
        for (int j = 0; j < jmax ; j++)
        {
            std::vector<int> matching;
            std::vector<std::vector<int>> res;
            long totalTimeForSMI;
            long totalTimeForBD;
            std::vector<std::vector<double>> twoDcloudA;

            for (int a = 0; a < i; a++)
            {
                double rx = ((double) rand() / (RAND_MAX));
                double ry = ((double) rand() / (RAND_MAX));
                double rz = ((double) rand() / (RAND_MAX));
                twoDcloudA.push_back(std::vector<double> {rx,ry,rz});

            }
            std::vector<std::vector<double>> twoDcloudB = twoDcloudA;

            PermuteAbit(twoDcloudB, 0.2);
            auto rng = std::default_random_engine {};

            arma::mat cloudA;
            arma::mat cloudB;
            arma::mat cloudC;

            Vec2Mat( twoDcloudA, cloudA);
            Vec2Mat( twoDcloudB, cloudB);

            std::shuffle(std::begin(twoDcloudB), std::end(twoDcloudB), rng);

            Vec2Mat(twoDcloudB, cloudC);


           // PermuteAbit(cloudB, 0.1);


            std::string folder = "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/ArmaFix7/";
            std::string saveFileA = "A" + std::to_string(i) + "_" + std::to_string(j) + ".mat";
            std::string saveFileB = "B" + std::to_string(i) + "_" + std::to_string(j) + ".mat";

          //!  cloudA.save(folder + saveFileA, arma::arma_ascii);
          //!  cloudB.save(folder + saveFileB, arma::arma_ascii);


            //     std::cout << "Cloud A: " << std::endl;
            //   std::cout << cloudA << std::endl;
            //   std::cout << "Cloud B: " << std::endl;
            //  std::cout << cloudB << std::endl;

            //  std::cout << "Columns cloud A: " << cloudA.n_cols << std::endl;
            //  std::cout << "Rows cloud A: " << cloudA.n_rows << std::endl;

           //! auto firstP = std::make_pair(cloudA, graphA);
           //! auto secondP = std::make_pair(cloudB, graphB);


            std::cout << "Starting to compute isometry" << std::endl;

            double dd = nestedNaive.Compute_Isometry(cloudA, cloudB, matching, totalTimeForSMI, totalTimeForBD);

            double dzd = nestedClassic.Compute_Isometry(twoDcloudA, twoDcloudB, matching, totalTimeForSMI, totalTimeForBD );

         //   double ddudd = nestedNaive.Compute_Isometry(cloudA, cloudC, matching, totalTimeForSMI, totalTimeForBD);


            std::cout << "Distance computed" << std::endl;

            //!double purebd = CrystData::computeEuclideanBottleneckDistanceBC(cloudA, cloudB);
            std::vector<int> matchingOne;
            double purebd = eucBd.compute(cloudA, cloudB, matchingOne);
            double rati = dd / purebd;
        //    if (rati > 2)
         //   {
         //   std::string tmplol = std::to_string(i) + "_" + std::to_string(j);
         //   failures.push_back(tmplol);

         //   }
            // double dd = pp.ComputeDistance(cloudA, cloudB, graphA, graphB, matching, res, totalTimeForSMI, totalTimeForBD);
            std::cout << i << " , " << j << " : origNaive : " << dd << " | time: " << totalTimeForBD  << " |  ratio " << rati << std::endl;
            std::cout << i << " , " << j << " : nestedClassic : " << dzd << std::endl;

           // std::cout << i << " , " << j << " : permute : " << ddudd << std::endl;

             //  long totalTimeForSMI;
           //  long totalTimeForBD;

//            double dd2 = origBoost.Compute_Isometry(firstP, secondP, matching, totalTimeForSMI, totalTimeForBD);
//            double rati2 = dd2 / purebd;
//
//
//            std::cout << i << " , " << j << " :origBoost : " << dd2 << " | time: " << totalTimeForBD  << " |  ratio " << rati2 << std::endl;
//
//            double dd3 = nestedNaive.Compute_Isometry(firstP, secondP, matching, totalTimeForSMI, totalTimeForBD);
//            double rati3 = dd3 / purebd;
//
//
//            std::cout << i << " , " << j << " : nestedNaive : " << dd3 << " | time: " << totalTimeForBD  << " |  ratio " << rati3 << std::endl;
//
//
//             double dd4 = nestedBoost.Compute_Isometry(firstP, secondP, matching, totalTimeForSMI, totalTimeForBD);
//            double rati4 = dd4 / purebd;
//
//
//            std::cout << i << " , " << j << " : nestedBoost : " << dd4 << " | time: " << totalTimeForBD  << " |  ratio " << rati4 << std::endl;





            //  if (dd2 != dd)
            //    {
            //      std::cout << " FAILURE " << i << j << " FAILURE " << std::endl;

            //  }



            //if (dd > 0.01)
            //  {
            //    failure = failure + 1;
//
            //   }



        }
        //   std::cout << i << " | Failed: " << failure << "/" << jmax << std::endl;

    }

    std::cout << "Listing failures : " << std::endl;
    for (int i = 0; i < failures.size(); i++)
    {
        std::cout << failures[i] << std::endl;


    }




}


static void AdequareTestingSCDVer3()
{
    std::cout << "Start" << std::endl;
    int jmax = 1000;
    int maxxx = 30;
    int minn = 10;
    int failureNaive = 0;
    int failureBoost = 0;
    int failureRtree = 0;
    double err = 0.5;
    bool exit = false;
    for(int i = minn; i < maxxx; i++)
    {


        std::vector<std::vector<bool>> graphA(i+1, std::vector<bool>(i+1, true));
        std::vector<std::vector<bool>> graphB = graphA;
        for (int j = 0; j < jmax ; j++)
        {
            std::vector<int> matching;
            std::vector<std::vector<int>> res;
            long totalTimeForSMI;
            long totalTimeForBD;
            std::vector<std::vector<double>> twoDcloudA;
            for (int a = 0; a < i; a++)
            {
                double rx = ((double) rand() / (RAND_MAX));
                double ry = ((double) rand() / (RAND_MAX));
                double rz = ((double) rand() / (RAND_MAX));
                double rw = ((double) rand() / (RAND_MAX));
                twoDcloudA.push_back(std::vector<double> {rx,ry,rz,rw});

            }
            std::vector<std::vector<double>> twoDcloudB = twoDcloudA;
            auto rng = std::default_random_engine {};
            std::shuffle(std::begin(twoDcloudB), std::end(twoDcloudB), rng);
            arma::mat cloudA;
            arma::mat cloudB;
            Vec2Mat( twoDcloudA, cloudA);
            Vec2Mat( twoDcloudB, cloudB);
            PermuteAbit(cloudB, err);
            BoostBottleneck boostBd;
            BottleneckCalculator<Bipartie_Graph,GraphTreeController> naiveBd;
            BottleneckCalculator<rtree_graph, rtreeController> rtreeBd;
            std::vector<int> matchingOne;
            std::vector<int> matchingTwo;
            std::vector<int> matchingThree;


            double dd = naiveBd.bottleneck_between_clouds(cloudA, cloudB, matchingOne);
            double dd2 = boostBd.bottleneck_between_clouds(cloudA, cloudB, matchingTwo);
            double dd3 =  rtreeBd.bottleneck_between_clouds(cloudA, cloudB, matchingThree);
            std::cout << " THRESH: " << err <<  " | RESULTS " << i << ", " << j << " d1: " << dd << " d2: " << dd2 <<  " d3: " << dd3 << std::endl;
            if (dd > err)
            {
                failureNaive = failureNaive + 1;
                exit = true;
                break;
            }
            if (dd2 > err)
            {

                failureBoost = failureBoost + 1;
                exit = true;
                break;
            }
            if (dd3 > err)
            {
                exit = true;
                failureRtree = failureRtree + 1;
                break;
            }

        }
        if (exit)
        {
            break;

        }

    }
    std::cout << "Finnishing everything" << std::endl;
    std::cout << "Failurenaive : " << failureNaive << std::endl;
    std::cout << "Failureboost: " <<  failureBoost << std::endl;
    std::cout << "Failurertree : " << failureRtree << std::endl;





}


static void crztDataExprTwo()
{


//!"/home/yury/LocalProjects/XYZProject/Graphs/"
    ExperimentRunner p("/home/yury/LocalProjects/XYZProject/MolCoords/" );


    std::string inputCsv = "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/PDDOutputs/pdd_upd.csv";


//p.runFirst(inputCsv, 1);

    // p.RunThisSCD(inputCsv, "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/SCDOutputs/LayeredSCD_strengthfix_ordered_9000.csv",
    //  "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/SCDOutputs/LayeredSCD_strength_fix_unorderd_9000.csv", 30, 9000);

   // p.SingleRun(28626, 28622);

   // p.RunThisParallelSCD(inputCsv, "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/SCDOutputs/NestedSCD_9000_FullPrecision.csv",
        //                "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/SCDOutputs/NestedSCD_9000_FullPrecision_unorder.csv", 30, 9000);

        p.RunThisParallelSCD(inputCsv, "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/PBIResults/Full_10000_March27_recompute_PBI_Nested.csv",
                        "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/PBIResults/Unordered_March27_10000_recompute_PBI_Nested.csv", 30, 10000);


}

//! "/home/yury/LocalProjects/XYZProject/Graphless/SDVUpd/PDDTests.csv"

void ComputeSingleExp()
{

//  ExperimentRunner p("/home/yury/LocalProjects/XYZProject/MolCoords/" );
ExperimentRunner p;
p.RunSpecial("/home/yury/LocalProjects/XYZProject/Graphless/SDVUpd/CollectionOfXYZ/0/dsgdb9nsd_023787.csv", "/home/yury/LocalProjects/XYZProject/Graphless/SDVUpd/CollectionOfXYZ/0/dsgdb9nsd_023788.csv");
    std::string inputCsv = "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/PDDOutputs/pdd_upd.csv";

   // p.SingleRun(28626, 28622);

}

static void crztDataExprThree()
{

//ExperimentRunner p("/home/yury/LocalProjects/January2024Project/mols/" );
  //  std::string inputCsv = "/home/yury/LocalProjects/January2024Project/emdpdd_limited_extraDetail.csv";

 //   p.SingleRun(0, inputCsv);

//!p.RunThisParallelSCD(inputCsv, "/home/yury/LocalProjects/January2024Project/full_ordered_500.csv","/home/yury/LocalProjects/January2024Project/full_unordered_500.csv", 30, 500);

}

static void launcherTwoDTests()
{

ExtraExperiments::experimentsConducting(500,"/home/yury/LocalProjects/AploadMay07/TwoDExperiments/Voltwo/");
//!ExtraExperiments::PerturbationExperiments(30,3,25,0.1,10,10,"/home/yury/LocalProjects/JanuaryLate2024Project/Graphs/3d_30_10_10_0.1_out.csv");
//std::string folder = "/home/yury/LocalProjects/JanuaryLate2024Project/GraphsFeb2/";
//std::cout << "Starting pertubation experiments" << std::endl;
//ExtraExperiments::PerturbationExperiments(30,2,25,0.5,10,10,folder + "distance_2d.csv");
//std::cout << "Pertubation experiments finished, starting time experiments" << std::endl;
//ExtraExperiments::TimeExperiments(30, 2, 10, 4, 3, 0.1, folder + "time_2d.csv");
//std::cout << "All jobs done" << std::endl;

//!ExtraExperiments::PerturbationExperiments(30,2,100,0.1,10,10,"/home/yury/LocalProjects/JanuaryLate2024Project/Graphs/2d_100_10_10_0.1_recompute.csv");

    //!ExtraExperiments::TimeExperiments(30, 2, 10, 20, 5 , 0.1, "/home/yury/LocalProjects/JanuaryLate2024Project/Graphs/twodtime.csv");
    //!ExtraExperiments::TimeExperiments(30, 3, 10, 4, 3, 0.1, "/home/yury/LocalProjects/JanuaryLate2024Project/Graphs/threedtime.csv");

}

static void tupleTest()
{

std::tuple<int, int, int> primitiveTuple;
std::get<0>(primitiveTuple) = 1;
std::get<1>(primitiveTuple) = 137;
std::get<2>(primitiveTuple) = 90;
std::cout << std::get<0>(primitiveTuple) << std::endl;
std::cout << std::get<1>(primitiveTuple) << std::endl;
std::cout << std::get<2>(primitiveTuple) << std::endl;


}

static void CompNewTests()
{
std::string inputCsv = "/home/yury/Downloads/NestedSCD_9000_FullPrecision_Graph.csv";
ExperimentRunner p("/home/yury/LocalProjects/XYZProject/MolCoords/", "/home/yury/LocalProjects/XYZProject/Graphs/");
//p.RunThisParallelSCD(inputCsv, "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/SCDOutputs/NestedSCD_9000_FullPrecision.csv",
        //                "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/SCDOutputs/NestedSCD_9000_FullPrecision_unorder.csv", 30, 9000);
p.RunThisParallelSCDUpdate(inputCsv, "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/SCDOutputs/NestedSCD_30_November1_FullPrecision.csv", 30, 30);

}

static void debugPerturb()
{

//ExtraExperiments::PerturbationExperiments(30,2,25,0.5,10,10,folder + "distance_2d.csv");
double  res;
double  ratiod;
ExtraExperiments::PerturbExperimentsSingle(2, 25,  0.1, 1, 1, res , ratiod);

}

static void JustTesting()
{

double  res = 0;
double  ratiod = 0;
bool correctMatching = true;
double realdist = 0;

int numcloud = 0;
while(correctMatching)
{

ExtraExperiments::PerturbExperimentsSingleNew(2, 4,  0.05, numcloud, numcloud, res , ratiod, correctMatching,realdist);
numcloud = numcloud + 1;
}
std::cout << "Res : " << res << " Ratiod: " << ratiod << std::endl;
std::cout << " real dist: " << realdist << std::endl;
std::cout << "Error in numcloud: " << numcloud << std::endl;

std::cout << "The result is : " << correctMatching << std::endl;



}

static void AnotherTestTest()
{

ExtraExperiments::PerturbSingleNCDSmart(2, 10,  0.05, 0, 0);


}

static void SampleTest()
{


ExtraExperiments::SampleTest("/home/yury/LocalProjects/dbg/Matrices/A0.mat", "/home/yury/LocalProjects/dbg/Matrices/B0.mat");
}

static void TestLul()
{

double res = 0;
 double ratiod = 0;
 double realDist = 0;
 double ratiodistbd = 0;
 long time1;
 long time2;
ExtraExperiments::PerturbSingleNCDSmartComplicated(2, 8,  0.05, 0, 0,  res,  ratiod, realDist, ratiodistbd, time1,time2);
std::cout << "Res : " << res << " ratio: " << ratiod << " realdist: " << realDist << " raitodistbd: " << ratiodistbd << std::endl;


}

static void nextTest()
{
ExtraExperiments::MatchingExperiments(30, 2, 10, 10, 5 , 0.05, "/home/yury/LocalProjects/FebruaryLate2024Project/twodmatch_evening.csv");

}

static void rotationFinderTest()
{

 arma::mat A(2,4);
 arma::mat B(2,4);
 A(1,1) = 0;
 A(1,2) = 0;
 A(1,3) = 0;
 A(0,0) = -2;
 A(0,1) = -1;
 A(0,2) = 1;
 A(0,3) = 2;

 B(1,0) = -1;
 B(1,1) = 0;
 B(1,2) = 2;
 B(1,3) = 3;

 B(0,1) = 0;
 B(0,2) = 0;
 B(0,3) = 0;

 std::cout << "Matrix A: " << std::endl;
 std::cout << A << std::endl;

 std::cout << "Matrix B: " << std::endl;
 std::cout << B << std::endl;

 arma::mat R;
 arma::vec t;

 ExtraExperiments::findRotation(A,B,R,t);

 std::cout << "Rotation matrix found: " << std::endl;
 std::cout << R << std::endl;
 std::cout << "Transpose vector: " << std::endl;
 std::cout << t << std::endl;

 arma::mat B2 = R*A;

 ExtraExperiments::ShiftBy(B2,t);

 std::cout << "B2 is calculated as: " << std::endl;
 std::cout << B2 << std::endl;

}

static void TestOfTestz()
{

double boxS = 10;
int dim = 3;
int numPoints = 10;
double curEpsion = 0.1;
 int repeat = 0;
 int numEpsilon = 0;

  double res = 0;
   double ratiod = 0;
    double realDist = 0;
    double ratiodistbd = 0;

    double totalTimeForSMIstdA = 0;
     double totalTimeForBDA = 0;

 ExtraExperiments::PerturbSingleNCDTest(dim, numPoints,  curEpsion, repeat, numEpsilon, res, ratiod, realDist, ratiodistbd, totalTimeForSMIstdA, totalTimeForBDA, boxS);

 std::cout << "Res: " << res << std::endl;
 std::cout << "Ratiod " << ratiod << std::endl;
 std::cout << "realDist " << realDist << std::endl;
 std::cout << "ratiodistbd " << ratiodistbd << std::endl;
 std::cout << "totalTimeSMI " << totalTimeForSMIstdA << std::endl;
 std::cout << "TotalTimeBDA " << totalTimeForBDA << std::endl;


}
static void rotTestAnother()
{

ExtraExperiments::rotationTest(2, 5);
}

static void finalTest()
{
 //!ExtraExperiments::MatchingExperiments(30, 2, 10, 10, 1 , 0.05, "/home/yury/LocalProjects/FebruaryLate2024Project/twodmatch_15_Feb.csv");

//!ExtraExperiments::MatchingExperiments2DBB(30, 2, 10, 10, 10, 0.01, 3, "/home/yury/LocalProjects/FebruaryLate2024Project/MassiveHeapMaps/DataFast/TestTwo/");
//!ExtraExperiments::MatchingExperiments2DBBA(30, 3, 10, 3, 10, 0.01, 3, "/home/yury/LocalProjects/FebruaryLate2024Project/MassiveHeapMaps/DataFast/TestFour/");
std::cout << "Starting 2D " << std::endl;
ExtraExperiments::MatchingExperiments2DBBA(30, 2, 10, 10, 10, 0.01, 3, "/home/yury/LocalProjects/FebruaryLate2024Project/MassiveHeapMaps/DataMarch27/2D/");
std::cout << "Starting 3D " << std::endl;
ExtraExperiments::MatchingExperiments2DBBA(30, 3, 10, 3, 10, 0.01, 3, "/home/yury/LocalProjects/FebruaryLate2024Project/MassiveHeapMaps/DataMarch27/3D/");
std::cout << "Experimenmts finnished" << std::endl;

//!ExtraExperiments::MatchingExperiments2DBB(30, 3, 10, 4, 10, 0.01, 3, "/home/yury/LocalProjects/FebruaryLate2024Project/MassiveHeapMaps/DataFast/TestFive/");

//ExtraExperiments::MatchingExperiments2D(30, 2, 10, 10, 1, 0.005, 1, "/home/yury/LocalProjects/FebruaryLate2024Project/MassiveHeapMaps/Data/");

}

void ComputeSingleExpFiveSix()
{

//  ExperimentRunner p("/home/yury/LocalProjects/XYZProject/MolCoords/" );
ExperimentRunner p;
p.SpecialzRun("/home/yury/LocalProjects/March25Debug/0/dsgdb9nsd_129001.xyz", "/home/yury/LocalProjects/March25Debug/0/dsgdb9nsd_129022.xyz");
    std::string inputCsv = "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/PDDOutputs/pdd_upd.csv";

   // p.SingleRun(28626, 28622);

}

void quickTestzzz()
{

BottleneckCalculator< Bipartie_Graph , GraphTreeController  > p;

 std::vector<int> matching;

 arma::mat A(2,3);
 arma::mat B(2,3);

 A.col(0) = arma::vec{1,1};
 A.col(1) = arma::vec{1,5};
 A.col(2) = arma::vec{3,3};

 B.col(0) = arma::vec{1.05,5.05};
 B.col(1) = arma::vec{0.95, 0.95};
 B.col(2) = arma::vec{3.05,3.05};

double aza = p.bottleneck_between_clouds( A,  B, matching);

std::cout << " bottleneck distance between cloud: " << aza << std::endl;


 for (int i =0 ; i < matching.size() ; i++)
 {
    std::cout << matching[i] << ", ";



 }


}

void runVerdict()
{

 // ExperimentRunner p("/home/yury/LocalProjects/AploadMay07/Clouds/Numpy/" );

    ExperimentRunner p;
    p.Read2DFilesIntoMatrix("/home/yury/LocalProjects/AploadMay07/Clouds/Numpy/");
    std::string inputCsv = "/home/yury/LocalProjects/AploadMay07/ResultsPDD/may21_atomnum_pdd.csv";
    std::string outputCSV = "/home/yury/LocalProjects/AploadMay07/ResultsNCD/May21_80k_recompute.csv";
    std::string outputLimitedCSV = "/home/yury/LocalProjects/AploadMay07/ResultsNCD/May21_limited_80k_recompute.csv";


//p.runFirst(inputCsv, 1);

    // p.RunThisSCD(inputCsv, "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/SCDOutputs/LayeredSCD_strengthfix_ordered_9000.csv",
    //  "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/SCDOutputs/LayeredSCD_strength_fix_unorderd_9000.csv", 30, 9000);

   // p.SingleRun(28626, 28622);

   // p.RunThisParallelSCD(inputCsv, "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/SCDOutputs/NestedSCD_9000_FullPrecision.csv",
        //                "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/SCDOutputs/NestedSCD_9000_FullPrecision_unorder.csv", 30, 9000);

        p.RunThisParallelSCD2D(inputCsv, outputCSV, outputLimitedCSV, 30, 80000);



}

void RunMatchingExperiments()
{

//! memory: "/home/yury/LocalProjects/AploadMay07/LengthAndSpeedTests/LengthVer2/out_match.csv";

int threads = 30;
int dim = 2;
int numPoints = 100;
int stepsEps = 10;
double stepSizeEps = 0.1;
int repeat = 10;
std::string outputCore = "/home/yury/Projects/CloudsProject/Output/SpeedAndPrecision/LengthVer2/Clouds/";
std::string outputDist = "/home/yury/Projects/CloudsProject/Output/SpeedAndPrecision/LengthVer2/out_dist.csv";
std::string outputFinalMatch = "/home/yury/Projects/CloudsProject/Output/SpeedAndPrecision/LengthVer2/out_match.csv";

ExtraExperiments::MatchingExperimentsMatch(threads, dim, numPoints, stepsEps, stepSizeEps, repeat, outputCore, outputDist, outputFinalMatch);

std::cout << "***** Matching, Time ***** fininshed!" << std::endl;
//! static void MatchingExperimentsTime(int threads, int dim, int numPoints, int stepPoints, double epsilon, int repeat, std::string outputCore, std::string outputDist)

//ExtraExperiments::MatchingExperimentsTime(threads, dim, 10, 10, 0.05, repeat, "/home/yury/Projects/CloudsProject/Output/SpeedAndPrecision/Time/Clouds/","/home/yury/Projects/CloudsProject/Output/SpeedAndPrecision/Time/out_time.csv" );

}

void GraphIsomorphismTest()
{

GraphExperiments exper;
std::string graphA = "/home/yury/LocalProjects/ReaGraphStuff/TestGraphs/FullGraph3.csv";
std::string graphB = "/home/yury/LocalProjects/ReaGraphStuff/TestGraphs/FullGraph3B.csv";
exper.launchTest(graphA, graphB);






}

void GraphTestz()
{
GraphExperiments exper;
std::string graphA = "/home/yury/LocalProjects/ReaGraphStuff/TestGraphs/debug/Graph4_1.csv";
//!std::string graphA = "/home/yury/LocalProjects/ReaGraphStuff/TestGraphs/Graph6.csv";
std::string ffm = "/home/yury/LocalProjects/ReaGraphStuff/TestGraphs/debug/Graph4debug/debug4_1_3.txt";
std::vector<std::vector<int>> edgesTurned;

//! Graph::edge_descriptor e = this->my_find_edge(1, 3, in);

//!edgesTurned.push_back(std::vector<int>{1,3});
edgesTurned.push_back(std::vector<int>{0,1});
exper.testAdjancency(graphA, ffm, edgesTurned, 3);


}

void runSCDSuperspecialTests()
{


 ExperimentRunner p("/home/yury/LocalProjects/XYZProject/MolCoords/");
    std::string inputCsv = "/home/yury/ProjectSeptember10/pdd_sorted_with_srv.csv";
    //! "/home/yury/LocalProjects/AploadMay07/ResultsPDD/may21_atomnum_pdd.csv";
    std::string outputCSV = "/home/yury/ProjectSeptember10/outscd_3000.csv";
    std::string outputLimitedCSV =  "/home/yury/ProjectSeptember10/outscd_limited_3000.csv";


//p.runFirst(inputCsv, 1);

    // p.RunThisSCD(inputCsv, "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/SCDOutputs/LayeredSCD_strengthfix_ordered_9000.csv",
    //  "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/SCDOutputs/LayeredSCD_strength_fix_unorderd_9000.csv", 30, 9000);

   // p.SingleRun(28626, 28622);

   // p.RunThisParallelSCD(inputCsv, "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/SCDOutputs/NestedSCD_9000_FullPrecision.csv",
        //                "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/SCDOutputs/NestedSCD_9000_FullPrecision_unorder.csv", 30, 9000);


      //  p.SingleRunSpezziale(28626, 28622);

    p.RunThisParallelSCD3DHard(inputCsv, outputCSV, outputLimitedCSV, 30, 3000);




}


void runSCDSuperspecialTestsTwo()
{


 ExperimentRunner p("/home/yury/LocalProjects/XYZProject/MolCoords/", "/home/yury/LocalProjects/XYZProject/Graphs/");
    std::string inputCsv = "/media/yury/Linux work/QM9Graph/PDDResults/pdd_atomnum.csv";
    //! "/home/yury/LocalProjects/AploadMay07/ResultsPDD/may21_atomnum_pdd.csv";
    std::string outputCSV = "/media/yury/Linux work/QM9Graph/NCD/NCD10000_Real.csv";
    std::string outputLimitedCSV =  "/media/yury/Linux work/QM9Graph/NCD/NCDLimited10000_Real.csv";


//p.runFirst(inputCsv, 1);

    // p.RunThisSCD(inputCsv, "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/SCDOutputs/LayeredSCD_strengthfix_ordered_9000.csv",
    //  "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/SCDOutputs/LayeredSCD_strength_fix_unorderd_9000.csv", 30, 9000);

   // p.SingleRun(28626, 28622);

   // p.RunThisParallelSCD(inputCsv, "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/SCDOutputs/NestedSCD_9000_FullPrecision.csv",
        //                "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/SCDOutputs/NestedSCD_9000_FullPrecision_unorder.csv", 30, 9000);


      //  p.SingleRunSpezziale(28626, 28622);
   //! p.SingleRunSpezzialeVerTwo(28626, 28622);
    p.SingleRunSpezzialeVerTwo(28626, 28622);
   //! p.RunThisParallelSCD3DHardVerTwo(inputCsv, outputCSV, outputLimitedCSV, 30, 10000);




}

void runSCDSuperspecialTestsThree()
{


 ExperimentRunner p("/home/yury/LocalProjects/XYZProject/MolCoords/", "/home/yury/LocalProjects/XYZProject/Graphs/");
    std::string inputCsv = "/media/yury/Linux work/QM9Graph/NCD/NCD10000_Real.csv";
    //! "/home/yury/LocalProjects/AploadMay07/ResultsPDD/may21_atomnum_pdd.csv";
    std::string outputCSV = "/media/yury/Linux work/QM9Graph/NCD/NCD60_TOP.csv";
    std::string outputLimitedCSV =  "/media/yury/Linux work/QM9Graph/NCD/NCDLimited60_TOP.csv";


//p.runFirst(inputCsv, 1);

    // p.RunThisSCD(inputCsv, "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/SCDOutputs/LayeredSCD_strengthfix_ordered_9000.csv",
    //  "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/SCDOutputs/LayeredSCD_strength_fix_unorderd_9000.csv", 30, 9000);

   // p.SingleRun(28626, 28622);

   // p.RunThisParallelSCD(inputCsv, "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/SCDOutputs/NestedSCD_9000_FullPrecision.csv",
        //                "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/SCDOutputs/NestedSCD_9000_FullPrecision_unorder.csv", 30, 9000);


      //  p.SingleRunSpezziale(28626, 28622);
   //! p.SingleRunSpezzialeVerTwo(28626, 28622);
   // p.SingleRunSpezzialeVerTwo(28626, 28622);


 //!   p.RunThisParallelSCD3DHardVerThree(inputCsv, outputCSV, outputLimitedCSV, 30, 60);
     p.RunThisParallelSCD3DHardVerThree(inputCsv, outputCSV, outputLimitedCSV, 30, 60);





}

void RunAdditionalTests()
{


ExperimentRunner p("/home/yury/LocalProjects/XYZProject/MolCoords/", "/home/yury/LocalProjects/XYZProject/Graphs/");

std::vector<std::vector<int>> inputPairs;
inputPairs.push_back(std::vector<int>{ 132851,132477});
inputPairs.push_back(std::vector<int>{123532,24546});
inputPairs.push_back(std::vector<int>{123532,24520});
inputPairs.push_back(std::vector<int>{123531,24512});

p.RunThisParallelSCD3DHardVerFour( inputPairs, 4, "/media/yury/Linux work/QM9Graph/NCD/NCDSpecialOut.csv");

}

void runVerdictTwoTwo()
{

 // ExperimentRunner p("/home/yury/LocalProjects/AploadMay07/Clouds/Numpy/" );

    ExperimentRunner p("/home/yury/LocalProjects/XYZProject/MolCoords/");
    //p.Read2DFilesIntoMatrix("/home/yury/LocalProjects/AploadMay07/Clouds/Numpy/");
    std::string inputCsv = "/media/yury/Linux work/QM9Graph/PDDResults/pdd_atomnum.csv";
    std::string outputCSV = "/media/yury/Linux work/QM9Graph/PBI/PBI_10k.csv";
    std::string outputLimitedCSV = "/media/yury/Linux work/QM9Graph/PBI_10k_limited.csv";


//p.runFirst(inputCsv, 1);

    // p.RunThisSCD(inputCsv, "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/SCDOutputs/LayeredSCD_strengthfix_ordered_9000.csv",
    //  "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/SCDOutputs/LayeredSCD_strength_fix_unorderd_9000.csv", 30, 9000);

   // p.SingleRun(28626, 28622);

   // p.RunThisParallelSCD(inputCsv, "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/SCDOutputs/NestedSCD_9000_FullPrecision.csv",
        //                "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/SCDOutputs/NestedSCD_9000_FullPrecision_unorder.csv", 30, 9000);

        p.RunThisParallelSCD2D(inputCsv, outputCSV, outputLimitedCSV, 30, 10000);



}

void GraphIsomorphismTestSuper()
{

GraphExperiments exper;
std::string graphA = "/home/yury/Projects/Graph2D/GraphsCSVClose/category_0/dsgdb9nsd_005350.npy.csv";
std::string graphB = "/home/yury/Projects/Graph2D/GraphsCSVClose/category_0/dsgdb9nsd_005354.npy.csv";

std::string graphAOne = "/home/yury/Projects/Graph2D/GraphsCSVClose/category_57/dsgdb9nsd_133321.npy.csv";
std::string graphBOne = "/home/yury/Projects/Graph2D/GraphsCSVClose/category_57/dsgdb9nsd_133328.npy.csv";
//exper.launchTest(graphA, graphB);
exper.launchSimpleTest(graphAOne, graphBOne);






}
void RunMassiveExperimentsTest()
{
GraphExperiments exper;
//!std::string xyzFiles = "/home/yury/Projects/Graph2D/CSV_Graphs/";
std::string xyzFiles = "/home/yury/Projects/Graph2D/NoHStuff/CSV_Graphs/";
//!std::string csvInput = "/home/yury/Projects/Graph2D/SRV/PDD.csv";
std::string csvInput = "/home/yury/Projects/Graph2D/NoHStuff/Results/PDD.csv";
//std::string csvOutput = "/home/yury/Projects/Graph2D/ExperimentResults/QM9GraphDistanceFlat_20000.csv";
std::string csvOutput = "/home/yury/Projects/Graph2D/NoHStuff/Results/QM9GraphDistanceFlat_G1_100k.csv";
std::string csvOutputVer2 = "/home/yury/Projects/Graph2D/NoHStuff/Results/QM9GraphDistanceFlat_G2_10k.csv";
int threads = 30;
int maxxnumm = 10000;
//exper.MassiveExperimentsRun(xyzFiles, csvOutput, csvOutputVer2, threads, maxxnumm);
exper.MassiveExperimentsRun(xyzFiles, csvInput, csvOutput, threads, maxxnumm, 14500);
//exper.SingularRun( xyzFiles ,  2805, 1727 );



}

void RunMassiveExperimentsTestV2()
{
GraphExperiments exper;
//!std::string xyzFiles = "/home/yury/Projects/Graph2D/CSV_Graphs/";
//!std::string xyzFiles = "/home/yury/Projects/Graph2D/NoHStuff/CSV_Graphs/";
std::string xyzFiles = "/home/yury/Projects/Graph2D/DataOriginal/CSV/";
//!std::string csvInput = "/home/yury/Projects/Graph2D/SRV/PDD.csv";
std::string csvInput = "/home/yury/Projects/Graph2D/ExpResultsVer3/PDD.csv";
//std::string csvOutput = "/home/yury/Projects/Graph2D/ExperimentResults/QM9GraphDistanceFlat_20000.csv";
std::string csvOutput = "/home/yury/Projects/Graph2D/ExpResultsVer3/NCD1_upd.csv";
std::string csvOutputVer2 = "/home/yury/Projects/Graph2D/ExpResultsVer3/NCD2_upd.csv";

int threads = 30;
//int maxxnumm = 16453;
int maxxnumm = 5000;
exper.MassiveExperimentsRun(xyzFiles, csvOutput, csvOutputVer2, threads, maxxnumm);
//exper.MassiveExperimentsRun(xyzFiles, csvInput, csvOutput, threads, maxxnumm);
//exper.SingularRun( xyzFiles ,  2805, 1727 );



}


void RunMassiveExperimentsTestV3()
{
GraphExperiments exper;
//!std::string xyzFiles = "/home/yury/Projects/Graph2D/CSV_Graphs/";
//!std::string xyzFiles = "/home/yury/Projects/Graph2D/NoHStuff/CSV_Graphs/";
//!std::string xyzFiles = "/home/yury/Projects/Graph2D/DataOriginal/CSV/";
std::string xyzFiles = "/home/yury/Projects/Graph2D/NoHStuffV2/DataOriginal/CSV/";
//!std::string csvInput = "/home/yury/Projects/Graph2D/SRV/PDD.csv";
std::string csvInput = "/home/yury/Projects/Graph2D/NoHStuffV2/Outputs/PDD.csv";
//std::string csvOutput = "/home/yury/Projects/Graph2D/ExperimentResults/QM9GraphDistanceFlat_20000.csv";
std::string csvOutput = "/home/yury/Projects/Graph2D/NoHStuffV2/Outputs/NCD1.csv";
std::string csvOutputVer2 = "/home/yury/Projects/Graph2D/NoHStuffV2/Outputs/NCD2.csv";

int threads = 30;
//int maxxnumm = 16453;
int maxxnumm = 96001;
//exper.MassiveExperimentsRun(xyzFiles, csvInput, csvOutput, threads, maxxnumm) ;
int maxxnummtwo = 20000;
exper.MassiveExperimentsRunTwo(xyzFiles, csvOutput, csvOutputVer2, threads, maxxnummtwo);
//exper.SingularRun( xyzFiles ,  2805, 1727 );



}

void RunMassiveExperimentsTestV5()
{

//! RunThisNCDSuperV1(std::string inputCsv, std::string outputCsv, int threads, int maxxnumm, std::vector<std::vector<std::vector<double>>> & matrices, int minnumm = 0, std::string header = "", int kk = 1)

std::string xyzFiles = "/home/yury/Projects/Graph2D/NoHStuffV2/DataOriginal/Cloud_CSV/";
std::string graphCSV = "/home/yury/Projects/Graph2D/NoHStuffV2/DataOriginal/Graph_CSV/";
//std::string graphCSV = "/home/yury/Projects/Graph2D/DataOriginal/FullGraph_CSV/";


ExperimentRunner p(xyzFiles,graphCSV);//!std::string xyzFiles = "/home/yury/Projects/Graph2D/CSV_Graphs/";




std::string csvInput = "/home/yury/Projects/Graph2D/NoHStuffV2/OutputsV2/PDD.csv";
std::string csvOutput = "/home/yury/Projects/Graph2D/NoHStuffV2/OutputsV2/NCD1.csv";
std::string csvOutputVer2 = "/home/yury/Projects/Graph2D/NoHStuffV2/OutputsV2/NCD2.csv";

std::string headerOne ="indx_A,indx_B,Number_Atoms,Time_SRV,SRV,SDV,SDV_time,emd/pdd,time_emd,time_NCD1(microseconds),NBD_NCD1_time(microseconds),NBM_NCD1";
std::string headerTwo ="indx_A,indx_B,Number_Atoms,Time_SRV,SRV,SDV,SDV_time,emd/pdd,time_emd,time_NCD1(microseconds),NBD_NCD1_time(microseconds),NBM_NCD1,time_NCD2(microseconds),NBD_NCD2_time(microseconds),NBM_NCD2";


int threads = 30;
//int maxxnumm = 16453;
int maxxnumm = 96001;
int maxxnummtwo = 20000;
// oid RunThisNCDSuperV1(std::string inputCsv, std::string outputCsv, int threads, int maxxnumm, int minnumm = 0, std::string header = "", int kk = 1)

//p.SingleRunSpezzialeVerTwo(340, 336);

p.RunThisNCDSuperV1(csvInput, csvOutput, threads, maxxnumm, 0, headerOne, 1);
p.RunThisNCDSuperV1(csvOutput, csvOutputVer2, threads, maxxnummtwo, 0, headerTwo, 2);

//exper.MassiveExperimentsRun(xyzFiles, csvInput, csvOutput, threads, maxxnumm) ;
//!exper.MassiveExperimentsRunTwo(xyzFiles, csvOutput, csvOutputVer2, threads, maxxnummtwo);
//exper.SingularRun( xyzFiles ,  2805, 1727 );

}

void RunMassiveExperimentsTestV6()
{

//ExperimentRunner p;

 //z.GrabTheThings();
 //z.WhatWhat();

std::string cloudInput = "/home/yury/Projects/Graph2D/ExpResultsV5/FurtherThings/Cloud_CSV/";
std::string graphInput = "/home/yury/Projects/Graph2D/ExpResultsV5/FurtherThings/Graph_CSV/";
ExperimentRunner p(cloudInput,graphInput);
// NCDGrabber z("orca");
// z.WhatWhat();
//ExperimentRunner p;
p.SaveInfos();


}

void RunMassiveExperimentsTestV4()
{

//! RunThisNCDSuperV1(std::string inputCsv, std::string outputCsv, int threads, int maxxnumm, std::vector<std::vector<std::vector<double>>> & matrices, int minnumm = 0, std::string header = "", int kk = 1)

std::string xyzFiles = "/home/yury/Projects/Graph2D/DataOriginal/Cloud_CSV/";
std::string graphCSV = "/home/yury/Projects/Graph2D/DataOriginal/Graph_CSV/";
//std::string graphCSV = "/home/yury/Projects/Graph2D/DataOriginal/FullGraph_CSV/";


ExperimentRunner p(xyzFiles,graphCSV);//!std::string xyzFiles = "/home/yury/Projects/Graph2D/CSV_Graphs/";




std::string csvInput = "/home/yury/Projects/Graph2D/ExpResultsVer4/PDD.csv";
std::string csvOutput = "/home/yury/Projects/Graph2D/ExpResultsVer4/NCD1.csv";
std::string csvOutputVer2 = "/home/yury/Projects/Graph2D/ExpResultsVer4/NCD2.csv";

std::string headerOne ="indx_A,indx_B,Number_Atoms,Time_SRV,SRV,SDV,SDV_time,emd/pdd,time_emd,time_NCD1(microseconds),NBD_NCD1_time(microseconds),NBM_NCD1";
std::string headerTwo ="indx_A,indx_B,Number_Atoms,Time_SRV,SRV,SDV,SDV_time,emd/pdd,time_emd,time_NCD1(microseconds),NBD_NCD1_time(microseconds),NBM_NCD1,time_NCD2(microseconds),NBD_NCD2_time(microseconds),NBM_NCD2";


int threads = 30;
//int maxxnumm = 16453;
int maxxnumm = 16453;
int maxxnummtwo = 5000;
// oid RunThisNCDSuperV1(std::string inputCsv, std::string outputCsv, int threads, int maxxnumm, int minnumm = 0, std::string header = "", int kk = 1)

//p.SingleRunSpezzialeVerTwo(340, 336);

//p.RunThisNCDSuperV1(csvInput, csvOutput, threads, maxxnumm, 0, headerOne, 1);
p.RunThisNCDSuperV1(csvOutput, csvOutputVer2, threads, maxxnummtwo, 0, headerTwo, 2);

//exper.MassiveExperimentsRun(xyzFiles, csvInput, csvOutput, threads, maxxnumm) ;
//!exper.MassiveExperimentsRunTwo(xyzFiles, csvOutput, csvOutputVer2, threads, maxxnummtwo);
//exper.SingularRun( xyzFiles ,  2805, 1727 );



}

void PrepareExperiments()
{
int lululu = 3290;
ExperimentRunner p("/home/yury/Projects/Graph2D/CSV_Clouds/","/home/yury/Projects/Graph2D/GraphStructure/" );
//p.RunThisParallelSCD3DHardVerTwo("/home/yury/Projects/Graph2D/SRV/PDD.csv","/home/yury/Projects/Graph2D/ExperimentResults/QM9GraphDistanceFlat_3290_binary.csv",
//"/home/yury/Projects/Graph2D/ExperimentResults/QM9GraphDistanceFlat_3290_limited_binary.csv", 30, 3290);

p.RunThisParallelSCD3DHardVerTwo("/home/yury/Projects/Graph2D/ExperimentResults/QM9GraphDistanceFlat_16453_binary.csv","/home/yury/Projects/Graph2D/ExperimentResults/QM9GraphDistanceFlat_3290_binary.csv",
"/home/yury/Projects/Graph2D/ExperimentResults/QM9GraphDistanceFlat_3290_limited_binary.csv", 30, 3290);

}


static void crztDataExprFourz()
{
    std::string path = "/media/yury/36b7eb78-3a80-44e9-9269-a30359bf49c3";


//!"/home/yury/LocalProjects/XYZProject/Graphs/"
    ExperimentRunner p(path + "/home/yury/LocalProjects/XYZProject/MolCoords/" );


    std::string inputCsv = path + "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/PDDOutputs/pdd_upd.csv";


//p.runFirst(inputCsv, 1);

    // p.RunThisSCD(inputCsv, "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/SCDOutputs/LayeredSCD_strengthfix_ordered_9000.csv",
    //  "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/SCDOutputs/LayeredSCD_strength_fix_unorderd_9000.csv", 30, 9000);

   // p.SingleRun(28626, 28622);

   // p.RunThisParallelSCD(inputCsv, "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/SCDOutputs/NestedSCD_9000_FullPrecision.csv",
        //                "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/SCDOutputs/NestedSCD_9000_FullPrecision_unorder.csv", 30, 9000);

        p.RunThisParallelSCD(inputCsv, "/home/yury/Projects/CloudsProject/Output/May23_TOP1000_NewInvariant.csv",
                        "/home/yury/Projects/CloudsProject/Output/Asdfds.csv", 30, 1000, 0);


}


static void crztDataExprFivez()
{
    std::string path = "/media/yury/36b7eb78-3a80-44e9-9269-a30359bf49c3";


//!"/home/yury/LocalProjects/XYZProject/Graphs/"
    ExperimentRunner p(path + "/home/yury/LocalProjects/XYZProject/MolCoords/" );


    std::string inputCsv = "/home/yury/Projects/CloudsProject/Data/HierarchicalProject/PDD/ToComputeBoys.csv";


//p.runFirst(inputCsv, 1);

    // p.RunThisSCD(inputCsv, "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/SCDOutputs/LayeredSCD_strengthfix_ordered_9000.csv",
    //  "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/SCDOutputs/LayeredSCD_strength_fix_unorderd_9000.csv", 30, 9000);

   // p.SingleRun(28626, 28622);

   // p.RunThisParallelSCD(inputCsv, "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/SCDOutputs/NestedSCD_9000_FullPrecision.csv",
        //                "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/SCDOutputs/NestedSCD_9000_FullPrecision_unorder.csv", 30, 9000);

        p.RunThisThingExperimentsMay(inputCsv, "/home/yury/Projects/CloudsProject/Data/HierarchicalProject/NDP/Then100-40000.csv",30, 100, 40000);


}

static void findMatching()
{
  //  std::string path = "/media/yury/36b7eb78-3a80-44e9-9269-a30359bf49c3";

//ExperimentRunner p(path + "/home/yury/LocalProjects/XYZProject/MolCoords/" );

ExperimentRunner p;
//p.runSingelton(171,172 );
//p.runSingelton(70953, 74129 );

p.runExtra("/home/yury/Projects/CloudsProject/Geometry/dsgdb9nsd_032787.xyz", "/home/yury/Projects/CloudsProject/Geometry/dsgdb9nsd_032730.xyz");


}



//! RunThisParallelSCD3DHardVerThree

//!      void limitedExperRunSeven(int x, int y,std::vector<std::vector<std::string>> & parsedCsv, std::vector<std::vector<std::string>> & outVector )
//!

int main(int argc, char *argv[])
{
//! /home/yury/LocalProjects/XYZProject/Graphless/SDVUpd/PDDTests.csv

 srand (time(NULL));
 //!RunMassiveExperimentsTestV6();
 //!RunMassiveExperimentsTest();
 //!RunMassiveExperimentsTestV2();
 //!RunMassiveExperimentsTestV3();
 //!RunMassiveExperimentsTestV4();
 //! RunMassiveExperimentsTestV5();
//!crztDataExprFourz();
//findMatching();
//!crztDataExprFourz();

 //findMatching();
 //crztDataExprFourz();
 //findMatching();
 //crztDataExprFivez();

 //!PrepareExperiments();

 //!RunMatchingExperiments();
 //!crztDataExprFourz();


// RunMassiveExperimentsTest();
 //!GraphIsomorphismTestSuper();
//! RunMatchingExperiments();
 //!GraphTestz();
//! runVerdictTwoTwo();
 //!RunAdditionalTests();
 //!runSCDSuperspecialTestsThree();
 //runSCDSuperspecialTests();
 //runSCDSuperspecialTestsTwo();
 //!GraphIsomorphismTest();
 //!GraphTestz();
 //!RunMatchingExperiments();
 //!ExtraExperiments::MatchingExperimentsMatch
 //!runVerdict();
 //!launcherTwoDTests();
 //!ComputeSingleExpFiveSix();

 //std::string firstA = argv[1];
 //std::string secondB = argv[2];
// std::cout << "Argumentsssss " << std::endl;
 //std::cout << firstA << std::endl;
 //std::cout << secondB << std::endl;
 //!ExperimentRunner p;

 //! quickTestzzz();

 //! p.SpecialzRun(firstA, secondB);
//! ComputeSingleExpFiveSix();

 //!finalTest();
 //!crztDataExprTwo();
 //!ComputeSingleExpFiveSix();
 //!finalTest();
// crztDataExprTwo();
//! finalTest();
// crztDataExprTwo();
 //!TestOfTestz();
     //!ComputeSingleExp();
   // SampleTest();
//    cout << "Hello world!" << endl;
 //   long p = 123345453;
  //  double q = p;
 //   std::cout << q << std::endl;
  //  std::cout << "xDouble = " << std::setprecision(15) << q << std::endl;

    //!finalTest();
  //!  TestOfTestz();
 //! crztDataExprTwo();
  //!  nextTest();
 //!   TestLul();
    //!rotTestAnother();
//!    rotationFinderTest();
//!TestLul();
//!nextTest();
   //! AnotherTestTest();
   //! JustTesting();
  //!  launcherTwoDTests();
    // debugPerturb();
   // launcherTwoDTests();
 //   launcherTwoDTests();
  //!  CompNewTests();
    //!ComputeSingleExp();
 //!   tupleTest();
   // AdequareTestingSCDVer2();

    //!launcherTwoDTests();

  //!  crztDataExprTwo();
    //!AdequareTestingSCDVer2();
    //!AdequareTestingSCDVer2();
   // crztDataExprTwo();
 //!   crztDataExprTwo();
 //! crztDataExprThree();
   //! AdequareTestingSCDVer3();
  //! AdequareTestingSCDVer2();
  // arma::vec p{3,0};
  // arma::vec q{0,4};
  // std::cout << arma::norm(p-q,"inf") << std::endl;
    return 0;
}
