#ifndef EXPERIMENTRUNNER_H
#define EXPERIMENTRUNNER_H

#include "BS_thread_pool.hpp"
#include <string>
#include "CloudIsometry/Invariants/PBI.h"
#include "CloudIsometry/Invariants/PBI_Update.h"
#include "CloudIsometry/Invariants/PBI_VD.h"
#include "CloudIsometry/Invariants/PBI_Old.h"
#include "CloudIsometry/Invariants/SCD_Updated.h"
#include "CloudIsometry/Invariants/NCDGrabber.h"
#include "CloudIsometry/Invariants/NCDGrab.h"
#include "CloudIsometry/Nested_Abstract.h"
#include "CloudIsometry/Nested_Updated.h"
#include "CloudIsometry/Invariants/NCD_Graph_Connections.h"
#include <iostream>
#include <filesystem>
#include "ExtraExperiments.h"


class ExperimentRunner
{
public:
    ExperimentRunner() {}



    ExperimentRunner(std::string xyzFiles)
    {

        namespace fs = std::filesystem;
        using namespace std;

        std::set<fs::path> sorted_by_name;

        for (auto &entry : fs::directory_iterator(xyzFiles))
            sorted_by_name.insert(entry.path());


        for (auto &filename : sorted_by_name)
        {
            std::string opp = filename.c_str();
            std::cout << filename.c_str() << std::endl;
            // std::vector<std::string> strings;
            // customSplit(opp, '/', strings);
            // std::string finall = strings[strings.size()-1];
            // out << "Finall (standard): " << finall << std::endl;
            // out << "Read path (standard): " << opp<< std::endl;
            //  out << "========" << std::endl;
            //   std::cout << opp << std::endl;
            arma::mat k;
            ReadFileIntoMatrix(opp, k);
            matrices.push_back(k);
            //!epicMap[finall] = k;

        }

        std::cout << "Finnished up with matrices " << std::endl;

        //  sorted_by_name = std::set<fs::path>();

        //  std::set<fs::path> sorted_by_name_graph;

//        for (auto &entry : fs::directory_iterator(graphFiles))
//            sorted_by_name.insert(entry.path());
//        bool bam = true;



        std::cout << "Starting graphs " << std::endl;

        for (int i = 0; i < matrices.size(); i++)
        {
            std::vector<std::vector<bool>> graphG(matrices[i].n_cols + 1, std::vector<bool>(matrices[i].n_cols + 1, true));
            graphs.push_back(graphG);

        }

//        for (auto &filename : sorted_by_name)
//        {
//            if (bam)
//            {
//                std::string opp = filename.c_str();
//                std::cout << filename.c_str() << std::endl;
//                std::vector<std::vector<bool>> graph;
//                std::vector<std::vector<bool>> graphTmp;
//                //!arma::mat k;
//                //!ReadFileIntoMatrix(opp, k);
//                ReadGraphFileIntoMatrix(opp, graphTmp);
//                tableExpansion(graphTmp, graph);
//
//
//                graphs.push_back(graph);
//
//                //!epicMap[finall] = k;
//            }
//
//        }
////
//        std::cout << "Finnished up everything" << std::endl;

    }

    ExperimentRunner(std::string xyzFiles, std::string graphFiles)
    {

        namespace fs = std::filesystem;
        using namespace std;

        std::set<fs::path> sorted_by_name;

        for (auto &entry : fs::directory_iterator(xyzFiles))
            sorted_by_name.insert(entry.path());


        for (auto &filename : sorted_by_name)
        {
            std::string opp = filename.c_str();
            std::cout << filename.c_str() << std::endl;
            // std::vector<std::string> strings;
            // customSplit(opp, '/', strings);
            // std::string finall = strings[strings.size()-1];
            // out << "Finall (standard): " << finall << std::endl;
            // out << "Read path (standard): " << opp<< std::endl;
            //  out << "========" << std::endl;
            //   std::cout << opp << std::endl;
            arma::mat k;
            //  ReadFileIntoMatrix(opp, k);
            ReadFileIntoMatrix2D(opp, k);
            cloudPaths.push_back(opp);
            matrices.push_back(k);
            //!epicMap[finall] = k;

        }

        std::cout << "Finnished u with matrices " << std::endl;

        sorted_by_name = std::set<fs::path>();

        //  std::set<fs::path> sorted_by_name_graph;

        for (auto &entry : fs::directory_iterator(graphFiles))
            sorted_by_name.insert(entry.path());
        bool bam = true;



        std::cout << "Starting graphs " << std::endl;

        for (auto &filename : sorted_by_name)
        {
            if (bam)
            {
                std::string opp = filename.c_str();
                std::cout << filename.c_str() << std::endl;
                std::vector<std::vector<bool>> graph;
                std::vector<std::vector<bool>> graphTmp;
                //!arma::mat k;
                //!ReadFileIntoMatrix(opp, k);
                ReadGraphFileIntoMatrix(opp, graphTmp);
                graphPaths.push_back(opp);
                tableExpansion(graphTmp, graph);


                graphs.push_back(graph);

                //!epicMap[finall] = k;
            }

        }
//
//        std::cout << "Finnished up everything" << std::endl;

    }

    void readCSV(std::string fname, std::vector<std::vector<double>> & matt)
    {

        namespace fs = std::filesystem;
        using namespace std;
        //  std::cout << "Beggining to read" << std::endl;

        std::vector<std::vector<std::string>> content;
        std::vector<string> row;
        std::vector<double> rowd;
        string line, word;

        std::cout << "Preparing to open " << std::endl;

        std::fstream file (fname, ios::in);

        std::cout << "Preparing to read: " << fname << std::endl;

        if(file.is_open())
        {
            std::cout << "Reading ... " << std::endl;
            while(getline(file, line))
            {
                std::cout << "Accesed " << std::endl;
                std::cout << "Reading line: " << line << std::endl;
                row.clear();
                rowd.clear();

                stringstream str(line);

                while(getline(str, word, ','))
                {
                    row.push_back(word);
                    rowd.push_back(stod(word));

                }
                content.push_back(row);
                matt.push_back(rowd);
            }
        }
        std::cout << "content: " << content.size() << " , " << content[0].size() << std::endl;


    }

    void Read2DFilesIntoMatrix(std::string xyzFiles)
    {


        namespace fs = std::filesystem;
        using namespace std;

        std::set<fs::path> sorted_by_name;

        for (auto &entry : fs::directory_iterator(xyzFiles))
            sorted_by_name.insert(entry.path());


        for (auto &filename : sorted_by_name)
        {
            std::string opp = filename.c_str();
            std::cout << filename.c_str() << std::endl;
            // std::vector<std::string> strings;
            // customSplit(opp, '/', strings);
            // std::string finall = strings[strings.size()-1];
            // out << "Finall (standard): " << finall << std::endl;
            // out << "Read path (standard): " << opp<< std::endl;
            //  out << "========" << std::endl;
            //   std::cout << opp << std::endl;
            arma::mat k;
            ReadFileIntoMatrix2D(opp, k);
            matrices.push_back(k);
            //!epicMap[finall] = k;

        }

        std::cout << "Finnished u pwith matrices " << std::endl;


    }


    void ReadFileIntoMatrix2D(std::string fname, arma::mat & matt)
    {

        namespace fs = std::filesystem;
        using namespace std;
        //  std::cout << "Beggining to read" << std::endl;

        std::vector<std::vector<std::string>> content;
        vector<string> row;
        string line, word;

        std::fstream file (fname, ios::in);
        if(file.is_open())
        {
            while(getline(file, line))
            {
                row.clear();

                stringstream str(line);

                while(getline(str, word, ','))
                    row.push_back(word);
                content.push_back(row);
            }
        }
        //  std::cout << content.size() << std::endl;
        matt.resize(2,content.size()-1);
        //  std::cout << "Matt cols : " << matt.n_cols << std::endl;
        // std::cout << "Matt rows : " << matt.n_rows << std::endl;

        for(int i=1; i<content.size(); i++)
        {
            for(int j=0; j<content[i].size(); j++)
            {
                matt.col(i-1)(j) = std::stod(content[i][j]);
            }
            // cout<<"\n";
        }
    }


    void ReadFileIntoMatrix(std::string fname, arma::mat & matt)
    {

        namespace fs = std::filesystem;
        using namespace std;
        //  std::cout << "Beggining to read" << std::endl;

        std::vector<std::vector<std::string>> content;
        vector<string> row;
        string line, word;

        std::fstream file (fname, ios::in);
        if(file.is_open())
        {
            while(getline(file, line))
            {
                row.clear();

                stringstream str(line);

                while(getline(str, word, ','))
                    row.push_back(word);
                content.push_back(row);
            }
        }
        //  std::cout << content.size() << std::endl;
        matt.resize(3,content.size()-1);
        //  std::cout << "Matt cols : " << matt.n_cols << std::endl;
        // std::cout << "Matt rows : " << matt.n_rows << std::endl;

        for(int i=1; i<content.size(); i++)
        {
            for(int j=0; j<content[i].size(); j++)
            {
                matt.col(i-1)(j) = std::stod(content[i][j]);
            }
            // cout<<"\n";
        }
    }

    void ReadGraphFileIntoMatrix(std::string fname, std::vector<std::vector<bool>> & matt)
    {
        namespace fs = std::filesystem;
        using namespace std;


        vector<vector<string>> content;
        vector<string> row;
        string line, word;
        bool skipfirstt = true;
        fstream file (fname, ios::in);
        if(file.is_open())
        {
            while(getline(file, line))
            {
                if (skipfirstt == true)
                {
                    skipfirstt = false;
                    continue;

                }
                row.clear();

                stringstream str(line);

                while(getline(str, word, ','))
                    row.push_back(word);
                content.push_back(row);
            }
        }
        matt = std::vector<std::vector<bool>>(content.size(),std::vector<bool>(content[0].size() - 1, false));
        //  std::cout << "Content size: " << content.size() << ", " <<content[0].size() << std::endl;
        // std::cout << "matt size " << matt.size() << std::endl;

        //    std::cout << "Matt cols : " << matt.n_cols << std::endl;
//       std::cout << "Matt rows : " << matt.n_rows << std::endl;

        for(int i=0; i<content.size(); i++)
        {
            for(int j=1; j<content[i].size(); j++)
            {
                // std::cout << content[i][j] << ", ";
                // std::cout << "(i,j)" << i << ", " << j << std::endl;
                if (std::stoi(content[i][j]) == 1)
                {
                    matt[i][j-1] = true;
                }
                else
                {
                    matt[i][j-1] = false;

                }
            }
            //std::cout << std::endl;
            // std::cout << " ========================= " << std::endl;
            // cout<<"\n";
        }
        //  std::cout << "Saved up succesfully" << std::endl;
    }

    void readPairs(std::vector<std::vector<std::string>> & stringPairs, std::string fname)
    {

        namespace fs = std::filesystem;
        using namespace std;

        vector<vector<string>> content;
        vector<string> row;
        string line, word;

        fstream file (fname, ios::in);
        if(file.is_open())
        {
            while(getline(file, line))
            {
                row.clear();

                stringstream str(line);

                while(getline(str, word, ','))
                    row.push_back(word);
                content.push_back(row);
            }
        }
        stringPairs = std::vector<std::vector<std::string>>(content.size() - 1,std::vector<std::string>(2));

        // std::cout << "Matt cols : " << matt.n_cols << std::endl;
        // std::cout << "Matt rows : " << matt.n_rows << std::endl;

        for(int i=1; i<content.size(); i++)
        {
            for(int j=1; j<content[i].size(); j++)
            {

                stringPairs[i-1][j-1] = content[i][j];
            }
            // cout<<"\n";
        }



    }

    void tableExpansion(std::vector<std::vector<bool>> & graphA, std::vector<std::vector<bool>> & graphOut)
    {
        // std::cout << "Table expansion started " << std::endl;
        // std::cout << "GraphA[0,0] " << graphA[0][0] << std::endl;
        // std::cout << "Graph A size: " << graphA.size() << std::endl;
        //std::cout << "Graph A size: " << graphA[0].size() << std::endl;
        //     for (int i = 0; i < graphA.size(); i++)
        //   {
        //       for (int j = 0; j < graphA[0].size(); j++)
        //      {
        // std::cout << i  << " , " << j << std::endl;
        //          std::cout << graphA[i][j] << ", ";
        //      }
        //     std::cout << std::endl;



        //  }
        graphOut = std::vector<std::vector<bool>>(graphA.size()+1, std::vector<bool>(graphA[0].size()+1, true));
        for (int i = 0; i < graphA.size(); i++)
        {
            for (int j = 0; j < graphA[0].size(); j++)
            {
                // std::cout << i  << " , " << j << std::endl;
                graphOut[i][j] = graphA[i][j];
            }



        }
        //  std::cout << "Table expansion ended " << std::endl;

    }

    template <typename T>
    std::string to_string_with_precision(const T a_value, const int n = 6)
    {
        std::ostringstream out;
        out.precision(n);
        out << std::fixed << a_value;
        return std::move(out).str();
    }

    void RunSpecial(std::string A, std::string B)
    {
        arma::mat kA;
        ReadFileIntoMatrix(A, kA);
        arma::mat kB;
        ReadFileIntoMatrix(B, kB);


        typedef BottleneckCalculator<Bipartie_Graph,GraphTreeController> bipartiebd;
        typedef bd_general_interface<bipartiebd, bipartiebd> naiveBdCalc;
        EuclideanBDCompute<bipartiebd> eucBd;


        //  std::vector<std::vector<std::string>> parsedCsv;
        // parseCSV(parsedCsv, inputCsv);

//        EuclideanBDCompute<bipartiebd> eucBd;
        //!      Nested_Base<SCD, naiveBdCalc, std::pair<arma::mat, std::vector<std::vector<bool>>>> nestedNaive;

        Nested_Base<PBI, naiveBdCalc, arma::mat> nestedNaive;
        Nested_Abstract<SCD_Updated, std::tuple<double,arma::mat, arma::mat>, naiveBdCalc, std::pair<arma::mat, std::vector<std::vector<bool>>>> nestedSCD;



        long totalTimeForSMIstd;
        long totalTimeForBD;
        std::vector<int> matching;
        std::vector<std::vector<int>> res;
        double bd = 0;

        std::vector<std::vector<bool>> graphA(kA.n_cols+1, std::vector<bool>(kA.n_cols+1, true));
        std::vector<std::vector<bool>> graphB(kB.n_cols+1, std::vector<bool>(kB.n_cols+1, true));

        auto firstP = std::make_pair(kA, graphA);
        auto secondP = std::make_pair(kB, graphB);



        // std::cout << "Clouds are : " << std::endl;
        //  std::cout << "Dimension of A: " << matrices[u].n_cols << ", " << matrices[u].n_rows << std::endl;
        //   std::cout << "Dimension of B: " << matrices[v].n_cols << ", " << matrices[v].n_rows << std::endl;



        std::vector<int> matchingOne;
        bd = nestedNaive.Compute_Isometry(kA, kB, matching, totalTimeForSMIstd, totalTimeForBD);
        //!bd = nestedSCD.Compute_Isometry(firstP,secondP, matching, totalTimeForSMIstd, totalTimeForBD,3);
        double purebd = eucBd.compute(kA, kB, matchingOne);


        //std::cout << "Computed isometry " << std::endl;

        std::cout << "Bottleneck was computed as " << bd << std::endl;
        std::cout << "Purebd was : " << purebd << std::endl;
        std::cout << "Ratio: " << bd / purebd << std::endl;


    }


    std::string get_filename_from_path(const std::string& full_path)
    {
        return std::filesystem::path(full_path).filename().string();
    }

    void SaveInfos()
    {

        //!     std::vector<std::string> cloudPaths;
//!    std::vector<std::string> graphPaths;

        //! std::string fp, std::string sp    std::vector<std::string> cloudPaths;
        //  std::vector<std::string> graphPaths;
        //!GrabTheThings(std::pair<arma::mat, std::vector<std::vector<bool>>>& cloudAunc,  int k, std::string pathOut)
        //!     std::vector<std::string> cloudPaths;
//!   std::vector<std::string> graphPaths;

//!     std::vector<arma::mat> matrices;
        //!  std::vector<std::vector<std::vector<bool>>> graphs;
        std::cout << "Reached this " << std::endl;
//NCDGrabber pq("orca");

        for (int i = 0; i < matrices.size(); i++ )
        {

            auto firstP = std::make_pair(matrices[i], graphs[i]);
            // auto secondP = std::make_pair(matrices[i], graphB);
              std::string A = get_filename_from_path(cloudPaths[i]);
            //  std::string B = get_filename_from_path(cloudPaths[i]);
            std::cout << "Starting " << i << std::endl;
             NCDGrabber pq("orca");
        pq.WhatWhat();
            // pq.WhatWhat();
         //   pq.GrabTheThings(firstP,  1, "/home/yury/Projects/Graph2D/ExpResultsV5/FurtherThings/Outputs/NCD1/" + A);
            std::cout << "lets go ! " << std::endl;
            //   p.WhatWhat();
              pq.GrabTheThings(firstP,  2, "/home/yury/Projects/Graph2D/ExpResultsV5/FurtherThings/Outputs/NCD2/" + A);
        }



    }

    template <typename T>
    void print2DVector(const std::vector<std::vector<T>>& vec)
    {
        for (const auto& row : vec)
        {
            for (const auto& item : row)
            {
                std::cout << item << " ";
            }
            std::cout << '\n';
        }
    }

    void SingleRunSpezzialeVerTwo(int u, int v)
    {

        typedef BottleneckCalculator<Bipartie_Graph,GraphTreeController> bipartiebd;
        typedef bd_general_interface<bipartiebd, bipartiebd> naiveBdCalc;

        //  std::vector<std::vector<std::string>> parsedCsv;
        // parseCSV(parsedCsv, inputCsv);

        EuclideanBDCompute<bipartiebd> eucBd;
        //!      Nested_Base<SCD, naiveBdCalc, std::pair<arma::mat, std::vector<std::vector<bool>>>> nestedNaive;

        Nested_Updated<NCD_Graph_Connections<naiveBdCalc>, std::pair<arma::mat, std::vector<std::vector<bool>>>> nestedNaive;

        //! std::vector<std::vector<std::string>> parsedCsv;

        long totalTimeForSMIstd;
        long totalTimeForBD;
        std::vector<int> matching;
        std::vector<std::vector<int>> res;
        double bd = 0;

        std::cout << " HERE WE ARE " << std::endl;

        auto firstP = std::make_pair(matrices[u], graphs[u]);
        auto secondP = std::make_pair(matrices[v], graphs[v]);

        std::cout << "U: " << u << " : " << "V: " << v << std::endl;

        std::cout << cloudPaths[u] << std::endl;
        std::cout << cloudPaths[v] << std::endl;

        std::cout << matrices[u]<< std::endl;
        std::cout << matrices[u].n_cols << std::endl;

        std::cout << " ==== " << std::endl;

        std::cout << matrices[v] << std::endl;
        std::cout << matrices[v].n_cols << std::endl;

        std::cout << " ==== " << std::endl;

        print2DVector(graphs[u]);

        std::cout << " ==== " << std::endl;

        print2DVector(graphs[v]);

        std::cout << " ==== " << std::endl;

        std::cout << graphPaths[u] << std::endl;
        std::cout << graphPaths[v] << std::endl;

        std::cout << "Starting to compute" << std::endl;



        // std::cout << "Clouds are : " << std::endl;
        //  std::cout << "Dimension of A: " << matrices[u].n_cols << ", " << matrices[u].n_rows << std::endl;
        //   std::cout << "Dimension of B: " << matrices[v].n_cols << ", " << matrices[v].n_rows << std::endl;




        bd = nestedNaive.Compute_Isometry(firstP, secondP, matching, totalTimeForSMIstd, totalTimeForBD, 1);



        //std::cout << "Computed isometry " << std::endl;

        std::cout << "Bottleneck was computed as " << bd << std::endl;

        std::string tmpd = std::to_string(bd);

        std::cout << "Its string analog is : " << tmpd << std::endl;

        std::string tmpda = to_string_with_precision<double>(bd, 15);

        std::cout << "With full precision it is : " << tmpda << std::endl;

//
//                double bdtruncate = -1;
//
//            long totalTimeForSMIstdfake;
//            long totalTimeForBDfake;
//
//
//            std::vector<int> matchingOne;
//            double purebd = eucBd.compute(matrices[u], matrices[v], matchingOne);
//            //!CrystData::computeEuclideanBottleneckDistanceBC(matrices[u], matrices[v]);
//
//
//            double haussdorff = computeHaudsorff(matrices[u], matrices[v]);
//
//            tmp.push_back(parsedCsv[i][1]);
//            tmp.push_back(parsedCsv[i][2]);
//            tmp.push_back(parsedCsv[i][3]);
//            tmp.push_back(std::to_string(haussdorff));
//            tmp.push_back(std::to_string(purebd));
//            tmp.push_back(parsedCsv[i][4]);
//            tmp.push_back(parsedCsv[i][0]);
//            tmp.push_back(std::to_string(bd));
//            double ratiod = bd / purebd;
//            tmp.push_back(std::to_string(ratiod));
//            tmp.push_back(std::to_string(totalTimeForSMIstd));
//            tmp.push_back(std::to_string(totalTimeForBD));
//            tmp.push_back(std::to_string(bdtruncate));
//
//
//
//            outVector.push_back(tmp);



    }

    void SingleRunSpezzialeVerThree(int u, int v, std::vector<std::vector<std::string>> & strP, int i)
    {

        typedef BottleneckCalculator<Bipartie_Graph,GraphTreeController> bipartiebd;
        typedef bd_general_interface<bipartiebd, bipartiebd> naiveBdCalc;

        std::cout << "Launching " << i << std::endl;

        //  std::vector<std::vector<std::string>> parsedCsv;
        // parseCSV(parsedCsv, inputCsv);

        EuclideanBDCompute<bipartiebd> eucBd;
        //!      Nested_Base<SCD, naiveBdCalc, std::pair<arma::mat, std::vector<std::vector<bool>>>> nestedNaive;

        Nested_Updated<NCD_Graph_Connections<naiveBdCalc>, std::pair<arma::mat, std::vector<std::vector<bool>>>> nestedNaive;

        //! std::vector<std::vector<std::string>> parsedCsv;

        long totalTimeForSMIstd;
        long totalTimeForBD;
        std::vector<int> matching;
        std::vector<std::vector<int>> res;
        double bd = 0;

        auto firstP = std::make_pair(matrices[u], graphs[u]);
        auto secondP = std::make_pair(matrices[v], graphs[v]);



        // std::cout << "Clouds are : " << std::endl;
        //  std::cout << "Dimension of A: " << matrices[u].n_cols << ", " << matrices[u].n_rows << std::endl;
        //   std::cout << "Dimension of B: " << matrices[v].n_cols << ", " << matrices[v].n_rows << std::endl;



        bd = nestedNaive.Compute_Isometry(firstP, secondP, matching, totalTimeForSMIstd, totalTimeForBD, 3);

        //std::cout << "Computed isometry " << std::endl;

        // std::cout << "Bottleneck was computed as " << bd << std::endl;

        std::string tmpd = std::to_string(bd);
//
        //   std::cout << "Its string analog is : " << tmpd << std::endl;

        std::string tmpda = to_string_with_precision<double>(bd, 15);

        std::cout << "Combo: (" << u << ", " << v << " ): With full precision it is : " << tmpda << std::endl;
        //! strP[i] = tmpda;

        std::vector<std::string> tmpVec;
        tmpVec.push_back(std::to_string(u));
        tmpVec.push_back(std::to_string(v));
        tmpVec.push_back(tmpda);
        strP[i] = tmpVec;






    }



    void SingleRunSpezziale(int u, int v)
    {

        typedef BottleneckCalculator<Bipartie_Graph,GraphTreeController> bipartiebd;
        typedef bd_general_interface<bipartiebd, bipartiebd> naiveBdCalc;

        //  std::vector<std::vector<std::string>> parsedCsv;
        // parseCSV(parsedCsv, inputCsv);

        EuclideanBDCompute<bipartiebd> eucBd;
        //!      Nested_Base<SCD, naiveBdCalc, std::pair<arma::mat, std::vector<std::vector<bool>>>> nestedNaive;

        Nested_Updated<NCD_Smart_Update<naiveBdCalc>, arma::mat> nestedNaive;

        //! std::vector<std::vector<std::string>> parsedCsv;

        long totalTimeForSMIstd;
        long totalTimeForBD;
        std::vector<int> matching;
        std::vector<std::vector<int>> res;
        double bd = 0;

        //! auto firstP = std::make_pair(matrices[u], graphs[u]);
        //!auto secondP = std::make_pair( matrices[v], graphs[v]);



        // std::cout << "Clouds are : " << std::endl;
        //  std::cout << "Dimension of A: " << matrices[u].n_cols << ", " << matrices[u].n_rows << std::endl;
        //   std::cout << "Dimension of B: " << matrices[v].n_cols << ", " << matrices[v].n_rows << std::endl;



        bd = nestedNaive.Compute_Isometry(matrices[u],matrices[v], matching, totalTimeForSMIstd, totalTimeForBD, 2);

        //std::cout << "Computed isometry " << std::endl;

        std::cout << "Bottleneck was computed as " << bd << std::endl;

        std::string tmpd = std::to_string(bd);

        std::cout << "Its string analog is : " << tmpd << std::endl;

        std::string tmpda = to_string_with_precision<double>(bd, 15);

        std::cout << "With full precision it is : " << tmpda << std::endl;

//
//                double bdtruncate = -1;
//
//            long totalTimeForSMIstdfake;
//            long totalTimeForBDfake;
//
//
//            std::vector<int> matchingOne;
//            double purebd = eucBd.compute(matrices[u], matrices[v], matchingOne);
//            //!CrystData::computeEuclideanBottleneckDistanceBC(matrices[u], matrices[v]);
//
//
//            double haussdorff = computeHaudsorff(matrices[u], matrices[v]);
//
//            tmp.push_back(parsedCsv[i][1]);
//            tmp.push_back(parsedCsv[i][2]);
//            tmp.push_back(parsedCsv[i][3]);
//            tmp.push_back(std::to_string(haussdorff));
//            tmp.push_back(std::to_string(purebd));
//            tmp.push_back(parsedCsv[i][4]);
//            tmp.push_back(parsedCsv[i][0]);
//            tmp.push_back(std::to_string(bd));
//            double ratiod = bd / purebd;
//            tmp.push_back(std::to_string(ratiod));
//            tmp.push_back(std::to_string(totalTimeForSMIstd));
//            tmp.push_back(std::to_string(totalTimeForBD));
//            tmp.push_back(std::to_string(bdtruncate));
//
//
//
//            outVector.push_back(tmp);



    }


    void SingleRun(int u, int v)
    {

        typedef BottleneckCalculator<Bipartie_Graph,GraphTreeController> bipartiebd;
        typedef bd_general_interface<bipartiebd, bipartiebd> naiveBdCalc;

        //  std::vector<std::vector<std::string>> parsedCsv;
        // parseCSV(parsedCsv, inputCsv);

        EuclideanBDCompute<bipartiebd> eucBd;
        //!      Nested_Base<SCD, naiveBdCalc, std::pair<arma::mat, std::vector<std::vector<bool>>>> nestedNaive;

        Nested_Base<PBI, naiveBdCalc, arma::mat> nestedNaive;

        //! std::vector<std::vector<std::string>> parsedCsv;

        long totalTimeForSMIstd;
        long totalTimeForBD;
        std::vector<int> matching;
        std::vector<std::vector<int>> res;
        double bd = 0;

        auto firstP = std::make_pair(matrices[u], graphs[u]);
        auto secondP = std::make_pair( matrices[v], graphs[v]);



        // std::cout << "Clouds are : " << std::endl;
        //  std::cout << "Dimension of A: " << matrices[u].n_cols << ", " << matrices[u].n_rows << std::endl;
        //   std::cout << "Dimension of B: " << matrices[v].n_cols << ", " << matrices[v].n_rows << std::endl;



        bd = nestedNaive.Compute_Isometry(matrices[u],matrices[v], matching, totalTimeForSMIstd, totalTimeForBD);

        //std::cout << "Computed isometry " << std::endl;

        std::cout << "Bottleneck was computed as " << bd << std::endl;

        std::string tmpd = std::to_string(bd);

        std::cout << "Its string analog is : " << tmpd << std::endl;

        std::string tmpda = to_string_with_precision<double>(bd, 15);

        std::cout << "With full precision it is : " << tmpda << std::endl;

//
//                double bdtruncate = -1;
//
//            long totalTimeForSMIstdfake;
//            long totalTimeForBDfake;
//
//
//            std::vector<int> matchingOne;
//            double purebd = eucBd.compute(matrices[u], matrices[v], matchingOne);
//            //!CrystData::computeEuclideanBottleneckDistanceBC(matrices[u], matrices[v]);
//
//
//            double haussdorff = computeHaudsorff(matrices[u], matrices[v]);
//
//            tmp.push_back(parsedCsv[i][1]);
//            tmp.push_back(parsedCsv[i][2]);
//            tmp.push_back(parsedCsv[i][3]);
//            tmp.push_back(std::to_string(haussdorff));
//            tmp.push_back(std::to_string(purebd));
//            tmp.push_back(parsedCsv[i][4]);
//            tmp.push_back(parsedCsv[i][0]);
//            tmp.push_back(std::to_string(bd));
//            double ratiod = bd / purebd;
//            tmp.push_back(std::to_string(ratiod));
//            tmp.push_back(std::to_string(totalTimeForSMIstd));
//            tmp.push_back(std::to_string(totalTimeForBD));
//            tmp.push_back(std::to_string(bdtruncate));
//
//
//
//            outVector.push_back(tmp);



    }


    void limitedExperRunTwo(int x, int y,std::vector<std::vector<std::string>> & parsedCsv, std::vector<std::vector<std::string>> & outVector )
    {
        //LayeredSCD lap(0.43,2,"");
        int tt = matrices.size();
        y = std::min(y, tt);
        //  std::string debugStream = "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/GraphTypeInvariants/EMDPDD/debug.txt";
        typedef BottleneckCalculator<Bipartie_Graph,GraphTreeController> bipartiebd;
        typedef bd_general_interface<bipartiebd, bipartiebd> naiveBdCalc;
        Nested_Abstract<SCD_Updated, std::tuple<double,arma::mat, arma::mat>, naiveBdCalc, std::pair<arma::mat, std::vector<std::vector<bool>>>> nestedSCD;

        EuclideanBDCompute<bipartiebd> eucBd;
        //!  Nested_Base<SCD, naiveBdCalc, std::pair<arma::mat, std::vector<std::vector<bool>>>> nestedNaive;
        Nested_Base<PBI, naiveBdCalc, arma::mat> nestedNaive;


        for (int i = x; i < y; i++)
        {
            std::cout << "Running operation (Full Graph): [" << x << " --- " << i << " --- " << y << "]" << std::endl;
            std::vector<std::string> tmp;
            int u = stoi(parsedCsv[i][0]);
            int v = stoi(parsedCsv[i][1]);
            long totalTimeForSMIstd;
            long totalTimeForBD;
            std::vector<int> matching;
            std::vector<std::vector<int>> res;
            double bd = 0;

            auto firstP = std::make_pair(matrices[u], graphs[u]);
            auto secondP = std::make_pair( matrices[v], graphs[v]);

            // std::cout << "Clouds are : " << std::endl;
            //  std::cout << "Dimension of A: " << matrices[u].n_cols << ", " << matrices[u].n_rows << std::endl;
            //   std::cout << "Dimension of B: " << matrices[v].n_cols << ", " << matrices[v].n_rows << std::endl;



            bd = nestedSCD.Compute_Isometry(firstP, secondP, matching, totalTimeForSMIstd, totalTimeForBD,3);

            //  bd = nestedNaive.Compute_Isometry(matrices[u], matrices[v], matching, totalTimeForSMIstd, totalTimeForBD);

            //std::cout << "Computed isometry " << std::endl;

            //   double bdtruncate = -1;

            // long totalTimeForSMIstdfake;
            // long totalTimeForBDfake;


            std::vector<int> matchingOne;
            double purebd = eucBd.compute(matrices[u], matrices[v], matchingOne);
            double ratioxD = bd / purebd;
            tmp.push_back(std::to_string(bd));
            tmp.push_back(std::to_string(ratioxD));
            tmp.push_back(std::to_string(totalTimeForSMIstd));
            tmp.push_back(std::to_string(totalTimeForBD));

            //!CrystData::computeEuclideanBottleneckDistanceBC(matrices[u], matrices[v]);


            //!      double haussdorff = computeHaudsorff(matrices[u], matrices[v]);

            //! to_string_with_precision<double>(haussdorff , 15)
//            tmp.push_back(parsedCsv[i][1]);
//            tmp.push_back(parsedCsv[i][2]);
//            tmp.push_back(parsedCsv[i][3]);
//            tmp.push_back(to_string_with_precision<double>(haussdorff, 15));
//            tmp.push_back(to_string_with_precision<double>(purebd, 15));
//            tmp.push_back(parsedCsv[i][4]);
//            tmp.push_back(parsedCsv[i][0]);
//            tmp.push_back(to_string_with_precision<double>(bd, 15));
//            double ratiod = bd / purebd;
//            tmp.push_back(to_string_with_precision<double>(ratiod, 15));
//            tmp.push_back(std::to_string(totalTimeForSMIstd));
//            tmp.push_back(std::to_string(totalTimeForBD));
//            tmp.push_back(std::to_string(bdtruncate));



            outVector.push_back(tmp);
        }

        // std::cout << "Outvector size: " << outVector.size() << std::endl;


    }

    void SpecialzRun(std::string Astr, std::string Bstr)
    {
        std::cout << Astr << std::endl;
        arma::mat kA;
        ReadFileIntoMatrix(Astr, kA);
        std::cout << Bstr << std::endl;
        arma::mat kB;
        ReadFileIntoMatrix(Bstr, kB);


        std::cout << "Matrix A: " << kA.n_cols << std::endl;

        std::cout << kA << std::endl;

        std::cout << "Matrix B: " << kB.n_cols << std::endl;

        std::cout << kB << std::endl;


        //LayeredSCD lap(0.43,2,"");
        //  int tt = matrices.size();
        //  y = std::min(y, tt);
        //  std::string debugStream = "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/GraphTypeInvariants/EMDPDD/debug.txt";
        typedef BottleneckCalculator<Bipartie_Graph,GraphTreeController> bipartiebd;
        typedef bd_general_interface<bipartiebd, bipartiebd> naiveBdCalc;

        EuclideanBDCompute<bipartiebd> eucBd;
        Nested_Updated<PBI_Update<naiveBdCalc>, arma::mat> nestedNaive;
        //!  Nested_Base<SCD, naiveBdCalc, std::pair<arma::mat, std::vector<std::vector<bool>>>> nestedNaive;
        //!    Nested_Base<PBI, naiveBdCalc, arma::mat> nestedNaive;



        long totalTimeForSMIstd;
        long totalTimeForBD;
        std::vector<int> matching;
        std::vector<std::vector<int>> res;
        double bd = 0;

        //!     auto firstP = std::make_pair(matrices[u], graphs[u]);
        //!    auto secondP = std::make_pair( matrices[v], graphs[v]);

        // std::cout << "Clouds are : " << std::endl;
        //  std::cout << "Dimension of A: " << matrices[u].n_cols << ", " << matrices[u].n_rows << std::endl;
        //   std::cout << "Dimension of B: " << matrices[v].n_cols << ", " << matrices[v].n_rows << std::endl;



        //!    bd = nestedNaive.Compute_Isometry(firstP, secondP, matching, totalTimeForSMIstd, totalTimeForBD);

        bd = nestedNaive.Compute_Isometry(kA, kB, matching, totalTimeForSMIstd, totalTimeForBD);

        //std::cout << "Computed isometry " << std::endl;

        double bdtruncate = -1;

        long totalTimeForSMIstdfake;
        long totalTimeForBDfake;


        std::vector<int> matchingOne;
        double purebd = eucBd.compute(kA, kB, matchingOne);
        //!CrystData::computeEuclideanBottleneckDistanceBC(matrices[u], matrices[v]);


        double haussdorff = computeHaudsorff(kA, kB);


        //! okay boomer:

        arma::mat AA;

        ExtraExperiments::RearrangeCloud(kA, AA, matching);

        arma::mat Rv;
        arma::vec tttv;

        ExtraExperiments::findRotation(AA, kB,  Rv, tttv);

        arma::mat AAA = Rv*AA;

        ExtraExperiments::ShiftBy(AAA,tttv);

        double realDist = 0;

        for (int i = 0; i < matching.size(); i++)
        {

            double distP = arma::norm(AAA.col(i)-kB.col(i), 2);

            realDist = std::max(realDist, distP);


        }




        double ratiodistbd = realDist / bd;
        double ratiod = bd / purebd;
        bdtruncate = ratiodistbd;




        std::cout << "Rotation matrix: " << std::endl;
        std::cout << Rv << std::endl;
        std::cout << "Translation vector: " << std::endl;
        std::cout << tttv << std::endl;
        std::cout << "Euclidean bottleneck distance: " << purebd << std::endl;
        std::cout << "Nested bottleneck distance: " << bd << std::endl;
        std::cout << "Ratio (NBM / Euclidean BD ): " << ratiod << std::endl;
        std::cout << "Real dist: " << realDist << std::endl;
        std::cout << "Ratio (Real dist / NBM): " << ratiodistbd << std::endl;
        std::cout << "Obtainead matching: " << std::endl;
        printVector(matching);




//        //! to_string_with_precision<double>(haussdorff , 15)
//        tmp.push_back(parsedCsv[i][1]);
//        tmp.push_back(parsedCsv[i][2]);
//        tmp.push_back(parsedCsv[i][3]);
//        tmp.push_back(to_string_with_precision<double>(haussdorff, 15));
//        tmp.push_back(to_string_with_precision<double>(purebd, 15));
//        tmp.push_back(parsedCsv[i][4]);
//        tmp.push_back(parsedCsv[i][0]);
//        tmp.push_back(to_string_with_precision<double>(bd, 15));
//        double ratiod = bd / purebd;
//        tmp.push_back(to_string_with_precision<double>(ratiod, 15));
//        tmp.push_back(std::to_string(totalTimeForSMIstd));
//        tmp.push_back(std::to_string(totalTimeForBD));
//        tmp.push_back(std::to_string(bdtruncate));
//
//
//
//        outVector.push_back(tmp);


        // std::cout << "Outvector size: " << outVector.size() << std::endl;


    }
    template<typename t>
    void printVector(std::vector<t> & in)
    {
        //  std::cout << "Size of vector: " << in.size() << std::endl;
        for (int i = 0; i < in.size(); i++)
        {

            std::cout << in[i] << ", ";
        }
        std::cout << std::endl;
    }

    //! LimitedExperRunSeven

    void limitedExperRunEight(int x, int y,std::vector<std::vector<std::string>> & parsedCsv, std::vector<std::vector<std::string>> & outVector )
    {
        //LayeredSCD lap(0.43,2,"");
        int tt = matrices.size();
        y = std::min(y, tt);
        //  std::string debugStream = "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/GraphTypeInvariants/EMDPDD/debug.txt";
        typedef BottleneckCalculator<Bipartie_Graph,GraphTreeController> bipartiebd;
        typedef bd_general_interface<bipartiebd, bipartiebd> naiveBdCalc;

        EuclideanBDCompute<bipartiebd> eucBd;
        Nested_Updated<NCD_Graph_Connections<naiveBdCalc>, std::pair<arma::mat, std::vector<std::vector<bool>>>> nestedNaive;
        Nested_Updated<NCD_Graph_Connections<naiveBdCalc>, std::pair<arma::mat, std::vector<std::vector<bool>>>> nestedNaiveTwo;
        //!  Nested_Base<SCD, naiveBdCalc, std::pair<arma::mat, std::vector<std::vector<bool>>>> nestedNaive;
        //!    Nested_Base<PBI, naiveBdCalc, arma::mat> nestedNaive;


        for (int i = x; i < y; i++)
        {
            std::cout << "Running operation (Full Graph): [" << x << " --- " << i << " --- " << y << "]" << std::endl;
            std::vector<std::string> tmp;
            int u = stoi(parsedCsv[i][1]);
            int v = stoi(parsedCsv[i][2]);
            long totalTimeForSMIstd;
            long totalTimeForBD;

            long totalTimeForSMIstdXX;
            long totalTimeForBDXX;
            std::vector<int> matching;
            std::vector<int> matchingXX;
            std::vector<std::vector<int>> res;
            double bd = 0;

            auto firstP = std::make_pair(matrices[u], graphs[u]);
            auto secondP = std::make_pair( matrices[v], graphs[v]);

            // std::cout << "Clouds are : " << std::endl;
            //  std::cout << "Dimension of A: " << matrices[u].n_cols << ", " << matrices[u].n_rows << std::endl;
            //   std::cout << "Dimension of B: " << matrices[v].n_cols << ", " << matrices[v].n_rows << std::endl;



            //!    bd = nestedNaive.Compute_Isometry(firstP, secondP, matching, totalTimeForSMIstd, totalTimeForBD);


            //   std::cout << "*Step one" << std::endl;

            double bud = nestedNaiveTwo.Compute_Isometry(firstP, secondP, matchingXX, totalTimeForSMIstdXX, totalTimeForBDXX, 1);

            //   std::cout << "*Step two " << std::endl;

            bd = nestedNaive.Compute_Isometry(firstP, secondP, matching, totalTimeForSMIstd, totalTimeForBD, 2);

            //   std::cout << "*Step three" << std::endl;

            //!std::cout << "bd computed as : " << "[" << x << " --- " << i << " --- " << y << "]: " << bd << std::endl;



            //std::cout << "Computed isometry " << std::endl;

            double bdtruncate = -1;

            long totalTimeForSMIstdfake;
            long totalTimeForBDfake;


            std::vector<int> matchingOne;



            double purebd = eucBd.compute(matrices[u], matrices[v], matchingOne);

            //  std::cout << "*Step four " << std::endl;
            //!CrystData::computeEuclideanBottleneckDistanceBC(matrices[u], matrices[v]);




            double haussdorff = computeHaudsorff(matrices[u], matrices[v]);

            //  std::cout << "*Step five" << std::endl;




            //! okay boomer:

//            arma::mat AA;
//
//            ExtraExperiments::RearrangeCloud(matrices[u], AA, matching);
//
//            arma::mat Rv;
//            arma::vec tttv;
//
//            ExtraExperiments::findRotation(AA, matrices[v],  Rv, tttv);
//
//            arma::mat AAA = Rv*AA;
//
//            ExtraExperiments::ShiftBy(AAA,tttv);

            double realDist = 0;

//            for (int i = 0; i < matching.size(); i++)
//            {
//
//                double distP = arma::norm(AAA.col(i)-matrices[v].col(i), 2);
//
//                realDist = std::max(realDist, distP);
//            }
//            double ratiodistbd = realDist / bd;

            double ratiodistbd = 0;
            bdtruncate = ratiodistbd;


            //! to_string_with_precision<double>(haussdorff , 15)
            tmp.push_back(parsedCsv[i][1]);
            tmp.push_back(parsedCsv[i][2]);
            tmp.push_back(parsedCsv[i][5]);
            tmp.push_back(to_string_with_precision<double>(haussdorff, 15));
            tmp.push_back(to_string_with_precision<double>(purebd, 15));
            tmp.push_back(parsedCsv[i][3]);
            tmp.push_back(parsedCsv[i][4]);
            tmp.push_back(parsedCsv[i][0]);
            tmp.push_back(to_string_with_precision<double>(bd, 15));
            tmp.push_back(to_string_with_precision<double>(bud, 15));
            double ratiod = bd / purebd;
            tmp.push_back(to_string_with_precision<double>(ratiod, 15));
            tmp.push_back(std::to_string(totalTimeForSMIstd));
            tmp.push_back(std::to_string(totalTimeForBD));
            tmp.push_back(std::to_string(totalTimeForSMIstdXX));
            tmp.push_back(std::to_string(totalTimeForBDXX));
            //!tmp.push_back(std::to_string(bdtruncate));
            //! totalTimeForSMIstdXX, totalTimeForBDXX

            outVector.push_back(tmp);

            //  std::cout << "Step six" << std::endl;

        }


        // std::cout << "Outvector size: " << outVector.size() << std::endl;


    }


    void limitedExperRunTen(int x, int y,std::vector<std::vector<std::string>> & parsedCsv, std::vector<std::vector<std::string>> & outVector )
    {
        //LayeredSCD lap(0.43,2,"");
        //    int tt = matrices.size();

        //y = std::min(y, tt);

        //  std::string debugStream = "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/GraphTypeInvariants/EMDPDD/debug.txt";
        typedef BottleneckCalculator<Bipartie_Graph,GraphTreeController> bipartiebd;
        typedef bd_general_interface<bipartiebd, bipartiebd> naiveBdCalc;

        EuclideanBDCompute<bipartiebd> eucBd;
        Nested_Updated<NCD_Graph_Connections<naiveBdCalc>, std::pair<arma::mat, std::vector<std::vector<bool>>>> nestedNaive;
        Nested_Updated<NCD_Graph_Connections<naiveBdCalc>, std::pair<arma::mat, std::vector<std::vector<bool>>>> nestedNaiveTwo;
        //!  Nested_Base<SCD, naiveBdCalc, std::pair<arma::mat, std::vector<std::vector<bool>>>> nestedNaive;
        //!    Nested_Base<PBI, naiveBdCalc, arma::mat> nestedNaive;


        for (int i = x; i < y; i++)
        {
            std::cout << "Running operation (Full Graph): [" << x << " --- " << i << " --- " << y << "]" << std::endl;
            std::vector<std::string> tmp;
            int u = stoi(parsedCsv[i][0]);
            int v = stoi(parsedCsv[i][1]);
            long totalTimeForSMIstd;
            long totalTimeForBD;

            long totalTimeForSMIstdXX;
            long totalTimeForBDXX;
            std::vector<int> matching;
            std::vector<int> matchingXX;
            std::vector<std::vector<int>> res;
            double bd = 0;

            auto firstP = std::make_pair(matrices[u], graphs[u]);
            auto secondP = std::make_pair( matrices[v], graphs[v]);



            double bud = nestedNaiveTwo.Compute_Isometry(firstP, secondP, matchingXX, totalTimeForSMIstdXX, totalTimeForBDXX, 1);


            bd = nestedNaive.Compute_Isometry(firstP, secondP, matching, totalTimeForSMIstd, totalTimeForBD, 2);



            double bdtruncate = -1;

            long totalTimeForSMIstdfake;
            long totalTimeForBDfake;


            std::vector<int> matchingOne;



            double purebd = eucBd.compute(matrices[u], matrices[v], matchingOne);




            double haussdorff = computeHaudsorff(matrices[u], matrices[v]);


            double realDist = 0;


            double ratiodistbd = 0;
            bdtruncate = ratiodistbd;


            //! to_string_with_precision<double>(haussdorff , 15)
            tmp.push_back(parsedCsv[i][0]);
            tmp.push_back(parsedCsv[i][1]);
            tmp.push_back(parsedCsv[i][2]);
            tmp.push_back(to_string_with_precision<double>(haussdorff, 15));
            tmp.push_back(to_string_with_precision<double>(purebd, 15));
            tmp.push_back(parsedCsv[i][3]);
            tmp.push_back(parsedCsv[i][4]);
            tmp.push_back(parsedCsv[i][5]);
            tmp.push_back(to_string_with_precision<double>(bd, 15));
            tmp.push_back(to_string_with_precision<double>(bud, 15));
            double ratiod = bd / purebd;
            tmp.push_back(to_string_with_precision<double>(ratiod, 15));
            tmp.push_back(std::to_string(totalTimeForSMIstd));
            tmp.push_back(std::to_string(totalTimeForBD));
            tmp.push_back(std::to_string(totalTimeForSMIstdXX));
            tmp.push_back(std::to_string(totalTimeForBDXX));
            //!tmp.push_back(std::to_string(bdtruncate));
            //! totalTimeForSMIstdXX, totalTimeForBDXX

            outVector.push_back(tmp);

            //  std::cout << "Step six" << std::endl;

        }


        // std::cout << "Outvector size: " << outVector.size() << std::endl;


    }
    void limitedExperRunEleven(int x, int y,std::vector<std::vector<std::string>> & parsedCsv, std::vector<std::vector<std::string>> & outVector )
    {
        //LayeredSCD lap(0.43,2,"");
        //    int tt = matrices.size();

        //y = std::min(y, tt);

        //  std::string debugStream = "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/GraphTypeInvariants/EMDPDD/debug.txt";
        typedef BottleneckCalculator<Bipartie_Graph,GraphTreeController> bipartiebd;
        typedef bd_general_interface<bipartiebd, bipartiebd> naiveBdCalc;

        EuclideanBDCompute<bipartiebd> eucBd;
        Nested_Updated<NCD_Graph_Connections<naiveBdCalc>, std::pair<arma::mat, std::vector<std::vector<bool>>>> nestedNaive;
        Nested_Updated<NCD_Graph_Connections<naiveBdCalc>, std::pair<arma::mat, std::vector<std::vector<bool>>>> nestedNaiveTwo;
        //!  Nested_Base<SCD, naiveBdCalc, std::pair<arma::mat, std::vector<std::vector<bool>>>> nestedNaive;
        //!    Nested_Base<PBI, naiveBdCalc, arma::mat> nestedNaive;


        for (int i = x; i < y; i++)
        {
            std::cout << "Running operation (Full Graph): [" << x << " --- " << i << " --- " << y << "]" << std::endl;
            std::vector<std::string> tmp;
            int u = stoi(parsedCsv[i][0]);
            int v = stoi(parsedCsv[i][1]);
            long totalTimeForSMIstd = -1;
            long totalTimeForBD = -1;

            long totalTimeForSMIstdXX;
            long totalTimeForBDXX;
            std::vector<int> matching;
            std::vector<int> matchingXX;
            std::vector<std::vector<int>> res;
            double bd = 0;

            auto firstP = std::make_pair(matrices[u], graphs[u]);
            auto secondP = std::make_pair( matrices[v], graphs[v]);



            double bud = nestedNaiveTwo.Compute_Isometry(firstP, secondP, matchingXX, totalTimeForSMIstdXX, totalTimeForBDXX, 1);


            //   bd = nestedNaive.Compute_Isometry(firstP, secondP, matching, totalTimeForSMIstd, totalTimeForBD, 2);
            bd = -1;



            double bdtruncate = -1;

            long totalTimeForSMIstdfake;
            long totalTimeForBDfake;


            std::vector<int> matchingOne;



            double purebd = eucBd.compute(matrices[u], matrices[v], matchingOne);




            double haussdorff = computeHaudsorff(matrices[u], matrices[v]);


            double realDist = 0;


            double ratiodistbd = 0;
            bdtruncate = ratiodistbd;


            //! to_string_with_precision<double>(haussdorff , 15)
            tmp.push_back(parsedCsv[i][0]);
            tmp.push_back(parsedCsv[i][1]);
            tmp.push_back(parsedCsv[i][2]);
            tmp.push_back(to_string_with_precision<double>(haussdorff, 15));
            tmp.push_back(to_string_with_precision<double>(purebd, 15));
            tmp.push_back(parsedCsv[i][3]);
            tmp.push_back(parsedCsv[i][4]);
            tmp.push_back(parsedCsv[i][5]);
            tmp.push_back(to_string_with_precision<double>(bd, 15));
            tmp.push_back(to_string_with_precision<double>(bud, 15));
            double ratiod = bud / purebd;
            tmp.push_back(to_string_with_precision<double>(ratiod, 15));
            tmp.push_back(std::to_string(totalTimeForSMIstd));
            tmp.push_back(std::to_string(totalTimeForBD));
            tmp.push_back(std::to_string(totalTimeForSMIstdXX));
            tmp.push_back(std::to_string(totalTimeForBDXX));
            //!tmp.push_back(std::to_string(bdtruncate));
            //! totalTimeForSMIstdXX, totalTimeForBDXX

            outVector.push_back(tmp);

            //  std::cout << "Step six" << std::endl;

        }


        // std::cout << "Outvector size: " << outVector.size() << std::endl;


    }

    void SingleExperimentRunSuperVer1(int i, std::vector<std::vector<std::string>> & parsedCsv, std::vector<std::vector<std::string>> & outVector, int kok)
    {
        //LayeredSCD lap(0.43,2,"");
        int tt = matrices.size();
        // y = std::min(y, tt);
        //  std::string debugStream = "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/GraphTypeInvariants/EMDPDD/debug.txt";
        typedef BottleneckCalculator<Bipartie_Graph,GraphTreeController> bipartiebd;
        typedef bd_general_interface<bipartiebd, bipartiebd> naiveBdCalc;

        EuclideanBDCompute<bipartiebd> eucBd;
        //!  Nested_Updated<NCD_Graph<naiveBdCalc>, std::vector<std::vector<double>>> nestedNaive;
        //!  Nested_Updated<NCD_Graph<naiveBdCalc>, std::vector<std::vector<double>>> nestedNaiveTwo;

        Nested_Updated<NCD_Graph_Connections<naiveBdCalc>, std::pair<arma::mat, std::vector<std::vector<bool>>>> nestedNaive;
        Nested_Updated<NCD_Graph_Connections<naiveBdCalc>, std::pair<arma::mat, std::vector<std::vector<bool>>>> nestedNaiveTwo;

        //    Nested_Updated<NCD_Graph_Connections<naiveBdCalc>, std::pair<arma::mat, std::vector<std::vector<bool>>>> nestedNaive;
        //    Nested_Updated<NCD_Graph_Connections<naiveBdCalc>, std::pair<arma::mat, std::vector<std::vector<bool>>>> nestedNaiveTwo;


        //!  Nested_Base<SCD, naiveBdCalc, std::pair<arma::mat, std::vector<std::vector<bool>>>> nestedNaive;


        //!    Nested_Base<PBI, naiveBdCalc, arma::mat> nestedNaive;



        std::cout << "Running operation (Full Graph): [ --- " << i << " --- ]" << std::endl;
        std::vector<std::string> tmp;
        int u = stoi(parsedCsv[i][0]);
        int v = stoi(parsedCsv[i][1]);

        std::cout << "Conducted on : " << u << " , " << v << std::endl;
        long totalTimeForSMIstd;
        long totalTimeForBD;

        long totalTimeForSMIstdXX;
        long totalTimeForBDXX;
        std::vector<int> matching;
        std::vector<int> matchingXX;
        std::vector<std::vector<int>> res;
        double bd = 0;

        //  auto firstP = std::make_pair(matrices[u], graphs[u]);
        //  auto secondP = std::make_pair( matrices[v], graphs[v]);
        //!  auto firstP = matrices[u];
        //!  auto secondP = matrices[v];
        auto firstP = std::make_pair(matrices[u], graphs[u]);
        auto secondP = std::make_pair(matrices[v], graphs[v]);


        // double bud = nestedNaiveTwo.Compute_Isometry(firstP, secondP, matchingXX, totalTimeForSMIstdXX, totalTimeForBDXX, 1);

        //std::cout << "Bud is computed : " << bud << std::endl;


        double bud= nestedNaive.Compute_Isometry(firstP, secondP, matching, totalTimeForSMIstd, totalTimeForBD, kok);
        ;
        double bdtruncate = -1;

        long totalTimeForSMIstdfake;
        long totalTimeForBDfake;


        std::vector<int> matchingOne;
        //  double purebd = eucBd.compute(matrices[u], matrices[v], matchingOne);
        //  double haussdorff = computeHaudsorff(matrices[u], matrices[v]);
        double realDist = 0;
        double ratiodistbd = 0;
        bdtruncate = ratiodistbd;


        for (int jj= 0; jj < parsedCsv[i].size(); jj++)
        {
            tmp.push_back(parsedCsv[i][jj]);
        }


        tmp.push_back(std::to_string(totalTimeForSMIstd));
        tmp.push_back(std::to_string(totalTimeForBD));
        tmp.push_back(to_string_with_precision<double>(bud, 15));
        //!tmp.push_back(std::to_string(bdtruncate));
        //! totalTimeForSMIstdXX, totalTimeForBDXX

        outVector.push_back(tmp);

        //  std::cout << "Step six" << std::endl;




        // std::cout << "Outvector size: " << outVector.size() << std::endl;


    }



    void RunThisNCDSuperV1(std::string inputCsv, std::string outputCsv, int threads, int maxxnumm, int minnumm = 0, std::string header = "", int kk = 1)
    {
        //  maxxnumm = std::min(maxxnumm, (int) matrices.size());
        BS::thread_pool pool(threads);
        std::cout << "Starting " << std::endl;
        std::cout << "path of inputcsv: " << inputCsv << std::endl;
        std::vector<std::vector<std::string>> parsedCSV;
        parseCSV(parsedCSV, inputCsv);
        parsedCSV.erase(parsedCSV.begin());
        parsedCSV.resize(maxxnumm);
        std::vector<std::vector<std::vector<std::string>>> vectorCollection(maxxnumm);


        //! limitedExperRun, this, i*be, (i+1)*be, std::ref(parsedCSV), std::ref(vectorCollection[i])
        std::cout << "Minnumm is : " << minnumm << std::endl;

        std::cout << "Maxxnumm is : " << maxxnumm << std::endl;
        std::cout << "ParsedCSV size is : " << parsedCSV.size() << std::endl;
        std::cout << "Number of matrices is : " << matrices.size() << std::endl;


        for (int i = minnumm; i < maxxnumm; i++)
        {

            // SingleExperimentRunVer2

            //!  limitedExperRun(i, (i+1), parsedCSV, vectorCollection[i]);
            //! pool.push_task(&GraphExperiments::SingleExperimentRun, this, i, std::ref(parsedCSV), std::ref(vectorCollection[i]), std::ref(matrices));
            pool.push_task(&ExperimentRunner::SingleExperimentRunSuperVer1, this, i, std::ref(parsedCSV), std::ref(vectorCollection[i]), kk);
            //!  GraphExperiments::SingleExperimentRun( i, parsedCSV, vectorCollection[i], matrices);
            //  pool.push_task([i, &parsedCSV, &vectorCollection[i]]
            //   {
            //      this->limitedExperRun(i, i+1,parsedCSV)
            //   }
            //  );

        }
        pool.wait_for_tasks();

        //!std::string header ="Molecule_A,Molecule_B,Num_Atoms,Hausdorff_distance,BD,L_inf on SRV, L_inf on SDV,EMD on PDD,NBM_NCD2, NBM_NCD1 ,NBM_NDP/BD,time_NDP(microseconds),NBD_time(microseconds),Match_dist/NBM_NDP";
        //! std::string header ="indx_A,indx_B,Number_Atoms,Time_SRV,SRV,SDV,SDV_time,emd/pdd,time_emd,time_NCD1(microseconds),NBD_NCD1_time(microseconds),NBM_NCD1";
        //std::string header ="indx_A,indx_B,Number_Atoms,Time_SRV,SRV,SDV,SDV_time,emd/pdd,time_emd,time_NCD1(microseconds),NBD_NCD1_time(microseconds),NBM_NCD1,time_NCD2(microseconds),NBD_NCD2_time(microseconds),NBM_NCD2";

        //! std::string header ="Bottleneck dist, Time SMI, Time Bottleneck, PDD dist, Molecule_A, Molecule_B, NumAtoms, SDV dist";
        std::vector<std::vector<std::string>> finalString;

        for (int i = 0; i < vectorCollection.size(); i++)
        {
            finalString.insert(finalString.end(), vectorCollection[i].begin(), vectorCollection[i].end());

        }

        std::cout << "Total final string size " << finalString.size() << std::endl;

        std::vector<std::vector<std::string>> finalStringLimited = finalString;

        std::sort(finalString.begin(), finalString.end(),
                  [](const std::vector<std::string> & a, const std::vector<std::string> & b) -> bool
        {
            double aa = stod(a[a.size()-1]);
            double bb = stod(b[b.size()-1]);
            return  aa < bb;
        });

        long sumsd = 0;
        long sumbd = 0;
        long avgsd =0;
        long avgbd = 0;



        for (int i = 0; i < finalString.size(); i++)
        {
            sumsd = sumsd + stol(finalString[i][1]);
            sumbd = sumbd + stol(finalString[i][2]);
        }

        avgsd = sumsd / (double)finalString.size();
        avgbd = sumbd / (double)finalString.size();

        printExperiments(finalString, header, outputCsv);

        std::cout << "Total time SMI: " << sumsd << "milliseconds" << std::endl;
        std::cout << "Total time BD: " << sumbd << "milliseconds" << std::endl;
        std::cout << "Avg time SMI: " << avgsd<< "milliseconds" << std::endl;
        std::cout << "Avg time BD: " << avgbd << "milliseconds" << std::endl;


        // finalString;
        //finalStringLimited.resize(500);

        //  printExperiments(finalStringLimited, header, outputCsvLimited);

//   double ukk = parsedCSV.size() / (double) threads;
        // int be = floor(ukk);
        // int iblan = 0;




    }




    void limitedExperRunTwelve(int x, int y,std::vector<std::vector<std::string>> & parsedCsv, std::vector<std::vector<std::string>> & outVector )
    {
        //LayeredSCD lap(0.43,2,"");
        //    int tt = matrices.size();

        //y = std::min(y, tt);

        //  std::string debugStream = "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/GraphTypeInvariants/EMDPDD/debug.txt";
        typedef BottleneckCalculator<Bipartie_Graph,GraphTreeController> bipartiebd;
        typedef bd_general_interface<bipartiebd, bipartiebd> naiveBdCalc;

        EuclideanBDCompute<bipartiebd> eucBd;
        Nested_Updated<NCD_Graph_Connections<naiveBdCalc>, std::pair<arma::mat, std::vector<std::vector<bool>>>> nestedNaive;
        Nested_Updated<NCD_Graph_Connections<naiveBdCalc>, std::pair<arma::mat, std::vector<std::vector<bool>>>> nestedNaiveTwo;
        //!  Nested_Base<SCD, naiveBdCalc, std::pair<arma::mat, std::vector<std::vector<bool>>>> nestedNaive;
        //!    Nested_Base<PBI, naiveBdCalc, arma::mat> nestedNaive;


        for (int i = x; i < y; i++)
        {
            std::cout << "Running operation (Full Graph): [" << x << " --- " << i << " --- " << y << "]" << std::endl;
            std::vector<std::string> tmp = parsedCsv[i];
            int u = stoi(parsedCsv[i][0]);
            int v = stoi(parsedCsv[i][1]);
            long totalTimeForSMIstd;
            long totalTimeForBD;

            long totalTimeForSMIstdXX;
            long totalTimeForBDXX;
            std::vector<int> matching;
            std::vector<int> matchingXX;
            std::vector<std::vector<int>> res;
            double bd = 0;

            auto firstP = std::make_pair(matrices[u], graphs[u]);
            auto secondP = std::make_pair( matrices[v], graphs[v]);



            //   double bud = nestedNaiveTwo.Compute_Isometry(firstP, secondP, matchingXX, totalTimeForSMIstdXX, totalTimeForBDXX, 1);


            bd = nestedNaive.Compute_Isometry(firstP, secondP, matching, totalTimeForSMIstd, totalTimeForBD, 2);
            //    bd = -1;



            double bdtruncate = -1;

            long totalTimeForSMIstdfake;
            long totalTimeForBDfake;


            std::vector<int> matchingOne;



            double purebd = eucBd.compute(matrices[u], matrices[v], matchingOne);




            //    double haussdorff = computeHaudsorff(matrices[u], matrices[v]);


            double realDist = 0;


            double ratiodistbd = 0;
            bdtruncate = ratiodistbd;


            tmp[8] = to_string_with_precision<double>(bd, 15);
            double ratiod = bd / purebd;
            tmp[10] = to_string_with_precision<double>(ratiod, 15);
            tmp[11] = std::to_string(totalTimeForSMIstd);
            tmp[12] = std::to_string(totalTimeForBD);



            //! to_string_with_precision<double>(haussdorff , 15)

            //!tmp.push_back(std::to_string(bdtruncate));
            //! totalTimeForSMIstdXX, totalTimeForBDXX

            outVector.push_back(tmp);

            //  std::cout << "Step six" << std::endl;

        }


        // std::cout << "Outvector size: " << outVector.size() << std::endl;


    }

    void limitedExperRunNine(int x, int y,std::vector<std::vector<std::string>> & parsedCsv, std::vector<std::vector<std::string>> & outVector )
    {
        //LayeredSCD lap(0.43,2,"");

        // std::cout << "STarting something " << x << std::endl;
        printVector(parsedCsv[x]);

        int tt = matrices.size();
        y = std::min(y, tt);
        //  std::string debugStream = "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/GraphTypeInvariants/EMDPDD/debug.txt";
        typedef BottleneckCalculator<Bipartie_Graph,GraphTreeController> bipartiebd;
        typedef bd_general_interface<bipartiebd, bipartiebd> naiveBdCalc;

        EuclideanBDCompute<bipartiebd> eucBd;
        Nested_Updated<NCD_Graph_Connections<naiveBdCalc>, std::pair<arma::mat, std::vector<std::vector<bool>>>> nestedNaive;
        Nested_Updated<NCD_Graph_Connections<naiveBdCalc>, std::pair<arma::mat, std::vector<std::vector<bool>>>> nestedNaiveTwo;
        //!  Nested_Base<SCD, naiveBdCalc, std::pair<arma::mat, std::vector<std::vector<bool>>>> nestedNaive;
        //!    Nested_Base<PBI, naiveBdCalc, arma::mat> nestedNaive;


        for (int i = x; i < y; i++)
        {
            std::cout << "Running operation (Full Graph): [" << x << " --- " << i << " --- " << y << "]" << std::endl;
            std::vector<std::string> tmp;
            int u = stoi(parsedCsv[i][0]);
            int v = stoi(parsedCsv[i][1]);
            long totalTimeForSMIstd;
            long totalTimeForBD;

            long totalTimeForSMIstdXX;
            long totalTimeForBDXX;
            std::vector<int> matching;
            std::vector<int> matchingXX;
            std::vector<std::vector<int>> res;
            double bd = 0;

            auto firstP = std::make_pair(matrices[u], graphs[u]);
            auto secondP = std::make_pair( matrices[v], graphs[v]);

            // std::cout << "Clouds are : " << std::endl;
            //  std::cout << "Dimension of A: " << matrices[u].n_cols << ", " << matrices[u].n_rows << std::endl;
            //   std::cout << "Dimension of B: " << matrices[v].n_cols << ", " << matrices[v].n_rows << std::endl;



            //!    bd = nestedNaive.Compute_Isometry(firstP, secondP, matching, totalTimeForSMIstd, totalTimeForBD);


            //   std::cout << "*Step one" << std::endl;

            //   double bud = nestedNaiveTwo.Compute_Isometry(firstP, secondP, matchingXX, totalTimeForSMIstdXX, totalTimeForBDXX, 1);

            //   std::cout << "*Step two " << std::endl;

            bd = nestedNaive.Compute_Isometry(firstP, secondP, matching, totalTimeForSMIstd, totalTimeForBD, 3);

            //   std::cout << "*Step three" << std::endl;

            //!std::cout << "bd computed as : " << "[" << x << " --- " << i << " --- " << y << "]: " << bd << std::endl;



            //std::cout << "Computed isometry " << std::endl;

            double bdtruncate = -1;

            long totalTimeForSMIstdfake;
            long totalTimeForBDfake;


            std::vector<int> matchingOne;



            double purebd = eucBd.compute(matrices[u], matrices[v], matchingOne);

            //  std::cout << "*Step four " << std::endl;
            //!CrystData::computeEuclideanBottleneckDistanceBC(matrices[u], matrices[v]);




            double haussdorff = computeHaudsorff(matrices[u], matrices[v]);

            //  std::cout << "*Step five" << std::endl;




            //! okay boomer:

//            arma::mat AA;
//
//            ExtraExperiments::RearrangeCloud(matrices[u], AA, matching);
//
//            arma::mat Rv;
//            arma::vec tttv;
//
//            ExtraExperiments::findRotation(AA, matrices[v],  Rv, tttv);
//
//            arma::mat AAA = Rv*AA;
//
//            ExtraExperiments::ShiftBy(AAA,tttv);

            double realDist = 0;

//            for (int i = 0; i < matching.size(); i++)
//            {
//
//                double distP = arma::norm(AAA.col(i)-matrices[v].col(i), 2);
//
//                realDist = std::max(realDist, distP);
//            }
//            double ratiodistbd = realDist / bd;

            double ratiodistbd = 0;
            bdtruncate = ratiodistbd;


            //! to_string_with_precision<double>(haussdorff , 15)

            for(int j = 0 ; j < 10 ; j++)
            {
                tmp.push_back(parsedCsv[i][j]);
            }
            tmp.push_back(to_string_with_precision<double>(bd, 15));
            double ratiod = bd / purebd;
            tmp.push_back(to_string_with_precision<double>(ratiod,15));
            for (int j = 10; j < parsedCsv[i].size(); j++)
            {
                tmp.push_back(parsedCsv[i][j]);
            }
//
//
//
//            tmp.push_back(parsedCsv[i][0]);
//            tmp.push_back(parsedCsv[i][2]);
//            tmp.push_back(parsedCsv[i][5]);
//            tmp.push_back(to_string_with_precision<double>(haussdorff, 15));
//            tmp.push_back(to_string_with_precision<double>(purebd, 15));
//            tmp.push_back(parsedCsv[i][3]);
//            tmp.push_back(parsedCsv[i][4]);
//            tmp.push_back(parsedCsv[i][0]);
//            tmp.push_back(to_string_with_precision<double>(bd, 15));
//            tmp.push_back(to_string_with_precision<double>(bud, 15));
//            double ratiod = bd / purebd;
//            tmp.push_back(to_string_with_precision<double>(ratiod, 15));
            tmp.push_back(std::to_string(totalTimeForSMIstd));
            tmp.push_back(std::to_string(totalTimeForBD));
            //  tmp.push_back(std::to_string(totalTimeForSMIstdXX));
            //  tmp.push_back(std::to_string(totalTimeForBDXX));
            //!tmp.push_back(std::to_string(bdtruncate));
            //! totalTimeForSMIstdXX, totalTimeForBDXX

            outVector.push_back(tmp);

            //  std::cout << "Step six" << std::endl;

        }


        // std::cout << "Outvector size: " << outVector.size() << std::endl;


    }





    void limitedExperRunSeven(int x, int y,std::vector<std::vector<std::string>> & parsedCsv, std::vector<std::vector<std::string>> & outVector )
    {
        //LayeredSCD lap(0.43,2,"");
        int tt = matrices.size();
        y = std::min(y, tt);
        //  std::string debugStream = "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/GraphTypeInvariants/EMDPDD/debug.txt";
        typedef BottleneckCalculator<Bipartie_Graph,GraphTreeController> bipartiebd;
        typedef bd_general_interface<bipartiebd, bipartiebd> naiveBdCalc;

        EuclideanBDCompute<bipartiebd> eucBd;
        Nested_Updated<NCD_Smart_Update<naiveBdCalc>, arma::mat> nestedNaive;
        Nested_Updated<NCD_Smart_Update<naiveBdCalc>, arma::mat> nestedNaiveTwo;
        //!  Nested_Base<SCD, naiveBdCalc, std::pair<arma::mat, std::vector<std::vector<bool>>>> nestedNaive;
        //!    Nested_Base<PBI, naiveBdCalc, arma::mat> nestedNaive;


        for (int i = x; i < y; i++)
        {
            std::cout << "Running operation (Full Graph): [" << x << " --- " << i << " --- " << y << "]" << std::endl;
            std::vector<std::string> tmp;
            int u = stoi(parsedCsv[i][1]);
            int v = stoi(parsedCsv[i][2]);
            long totalTimeForSMIstd;
            long totalTimeForBD;

            long totalTimeForSMIstdXX;
            long totalTimeForBDXX;
            std::vector<int> matching;
            std::vector<int> matchingXX;
            std::vector<std::vector<int>> res;
            double bd = 0;

            //!     auto firstP = std::make_pair(matrices[u], graphs[u]);
            //!    auto secondP = std::make_pair( matrices[v], graphs[v]);

            // std::cout << "Clouds are : " << std::endl;
            //  std::cout << "Dimension of A: " << matrices[u].n_cols << ", " << matrices[u].n_rows << std::endl;
            //   std::cout << "Dimension of B: " << matrices[v].n_cols << ", " << matrices[v].n_rows << std::endl;



            //!    bd = nestedNaive.Compute_Isometry(firstP, secondP, matching, totalTimeForSMIstd, totalTimeForBD);

            double bud = nestedNaiveTwo.Compute_Isometry(matrices[u], matrices[v], matchingXX, totalTimeForSMIstdXX, totalTimeForBDXX, 1);



            bd = nestedNaive.Compute_Isometry(matrices[u], matrices[v], matching, totalTimeForSMIstd, totalTimeForBD, 2);

            //!std::cout << "bd computed as : " << "[" << x << " --- " << i << " --- " << y << "]: " << bd << std::endl;



            //std::cout << "Computed isometry " << std::endl;

            double bdtruncate = -1;

            long totalTimeForSMIstdfake;
            long totalTimeForBDfake;


            std::vector<int> matchingOne;



            double purebd = eucBd.compute(matrices[u], matrices[v], matchingOne);
            //!CrystData::computeEuclideanBottleneckDistanceBC(matrices[u], matrices[v]);




            double haussdorff = computeHaudsorff(matrices[u], matrices[v]);




            //! okay boomer:

//            arma::mat AA;
//
//            ExtraExperiments::RearrangeCloud(matrices[u], AA, matching);
//
//            arma::mat Rv;
//            arma::vec tttv;
//
//            ExtraExperiments::findRotation(AA, matrices[v],  Rv, tttv);
//
//            arma::mat AAA = Rv*AA;
//
//            ExtraExperiments::ShiftBy(AAA,tttv);

            double realDist = 0;

//            for (int i = 0; i < matching.size(); i++)
//            {
//
//                double distP = arma::norm(AAA.col(i)-matrices[v].col(i), 2);
//
//                realDist = std::max(realDist, distP);
//            }
//            double ratiodistbd = realDist / bd;

            double ratiodistbd = 0;
            bdtruncate = ratiodistbd;


            //! to_string_with_precision<double>(haussdorff , 15)
            tmp.push_back(parsedCsv[i][1]);
            tmp.push_back(parsedCsv[i][2]);
            tmp.push_back(to_string_with_precision<double>(haussdorff, 15));
            tmp.push_back(to_string_with_precision<double>(purebd, 15));
            tmp.push_back(parsedCsv[i][5]);
            tmp.push_back(parsedCsv[i][3]);
            tmp.push_back(parsedCsv[i][4]);
            tmp.push_back(parsedCsv[i][0]);
            tmp.push_back(to_string_with_precision<double>(bd, 15));
            tmp.push_back(to_string_with_precision<double>(bud, 15));
            double ratiod = bd / purebd;
            tmp.push_back(to_string_with_precision<double>(ratiod, 15));
            tmp.push_back(std::to_string(totalTimeForSMIstd));
            tmp.push_back(std::to_string(totalTimeForBD));
            tmp.push_back(std::to_string(bdtruncate));


            outVector.push_back(tmp);
        }


        // std::cout << "Outvector size: " << outVector.size() << std::endl;


    }


    void limitedExperRunSix(int x, int y,std::vector<std::vector<std::string>> & parsedCsv, std::vector<std::vector<std::string>> & outVector )
    {
        //LayeredSCD lap(0.43,2,"");
        int tt = matrices.size();
        y = std::min(y, tt);
        //  std::string debugStream = "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/GraphTypeInvariants/EMDPDD/debug.txt";
        typedef BottleneckCalculator<Bipartie_Graph,GraphTreeController> bipartiebd;
        typedef bd_general_interface<bipartiebd, bipartiebd> naiveBdCalc;

        EuclideanBDCompute<bipartiebd> eucBd;
        Nested_Updated<PBI_Update<naiveBdCalc>, arma::mat> nestedNaive;
        //!  Nested_Base<SCD, naiveBdCalc, std::pair<arma::mat, std::vector<std::vector<bool>>>> nestedNaive;
        //!    Nested_Base<PBI, naiveBdCalc, arma::mat> nestedNaive;


        for (int i = x; i < y; i++)
        {
            std::cout << "Running operation (Full Graph): [" << x << " --- " << i << " --- " << y << "]" << std::endl;
            std::vector<std::string> tmp;
            int u = stoi(parsedCsv[i][1]);
            int v = stoi(parsedCsv[i][2]);
            long totalTimeForSMIstd;
            long totalTimeForBD;
            std::vector<int> matching;
            std::vector<std::vector<int>> res;
            double bd = 0;

            //!     auto firstP = std::make_pair(matrices[u], graphs[u]);
            //!    auto secondP = std::make_pair( matrices[v], graphs[v]);

            // std::cout << "Clouds are : " << std::endl;
            //  std::cout << "Dimension of A: " << matrices[u].n_cols << ", " << matrices[u].n_rows << std::endl;
            //   std::cout << "Dimension of B: " << matrices[v].n_cols << ", " << matrices[v].n_rows << std::endl;



            //!    bd = nestedNaive.Compute_Isometry(firstP, secondP, matching, totalTimeForSMIstd, totalTimeForBD);

            bd = nestedNaive.Compute_Isometry(matrices[u], matrices[v], matching, totalTimeForSMIstd, totalTimeForBD);

            std::cout << " BD computed as : " << bd << std::endl;

            //std::cout << "Computed isometry " << std::endl;

            double bdtruncate = -1;

            long totalTimeForSMIstdfake;
            long totalTimeForBDfake;


            std::vector<int> matchingOne;
            double purebd = eucBd.compute(matrices[u], matrices[v], matchingOne);
            //!CrystData::computeEuclideanBottleneckDistanceBC(matrices[u], matrices[v]);


            double haussdorff = computeHaudsorff(matrices[u], matrices[v]);


            //! okay boomer:

            double realDist = 0;


//            arma::mat AA;
//
//            ExtraExperiments::RearrangeCloud(matrices[u], AA, matching);
//
//            arma::mat Rv;
//            arma::vec tttv;
//
//            ExtraExperiments::findRotation(AA, matrices[v],  Rv, tttv);
//
//            arma::mat AAA = Rv*AA;
//
//            ExtraExperiments::ShiftBy(AAA,tttv);
//
//

//            for (int i = 0; i < matching.size(); i++)
//            {
//
//                double distP = arma::norm(AAA.col(i)-matrices[v].col(i), 2);
//
//                realDist = std::max(realDist, distP);
//            }
//



            //  double ratiodistbd = realDist / bd;
            double ratiodistbd = 0;
            bdtruncate = ratiodistbd;

            //! to_string_with_precision<double>(haussdorff , 15)
            tmp.push_back(parsedCsv[i][1]);
            tmp.push_back(parsedCsv[i][2]);
            tmp.push_back(parsedCsv[i][5]);
            tmp.push_back(to_string_with_precision<double>(haussdorff, 15));
            tmp.push_back(to_string_with_precision<double>(purebd, 15));
            tmp.push_back(parsedCsv[i][3]);
            tmp.push_back(parsedCsv[i][4]);
            tmp.push_back(parsedCsv[i][0]);
            tmp.push_back(to_string_with_precision<double>(bd, 15));
            double ratiod = bd / purebd;
            tmp.push_back(to_string_with_precision<double>(ratiod, 15));
            tmp.push_back(std::to_string(totalTimeForSMIstd));
            tmp.push_back(std::to_string(totalTimeForBD));
            tmp.push_back(std::to_string(bdtruncate));



            outVector.push_back(tmp);
        }

        // std::cout << "Outvector size: " << outVector.size() << std::endl;


    }

    void limitedExperThirteen(int x, std::vector<std::vector<std::string>> & parsedCsv, std::vector<std::vector<std::string>> & outVector )
    {
        typedef BottleneckCalculator<Bipartie_Graph,GraphTreeController> bipartiebd;
        typedef bd_general_interface<bipartiebd, bipartiebd> naiveBdCalc;
        Nested_Updated<PBI_Update<naiveBdCalc>, arma::mat> nestedNaive;

        std::cout << "Running operation (Full Graph): [ " << x <<  " ]" << std::endl;

        int u = stoi(parsedCsv[x][0]);
        int v = stoi(parsedCsv[x][1]);
        long totalTimeForSMIstd;
        long totalTimeForBD;
        std::vector<int> matching;
        std::vector<std::vector<int>> res;
        double bd = nestedNaive.Compute_Isometry(matrices[u], matrices[v], matching, totalTimeForSMIstd, totalTimeForBD);


        std::vector<std::string> tmp;


        tmp.push_back(parsedCsv[x][0]);
        tmp.push_back(parsedCsv[x][1]);
        tmp.push_back(to_string_with_precision<double>(bd, 15));


        outVector.push_back(tmp);



    }


    void runExtra(std::string p1, std::string p2)
    {


        std::string opp = p1.c_str();
        std::cout << p1.c_str() << std::endl;
        // std::vector<std::string> strings;
        // customSplit(opp, '/', strings);
        // std::string finall = strings[strings.size()-1];
        // out << "Finall (standard): " << finall << std::endl;
        // out << "Read path (standard): " << opp<< std::endl;
        //  out << "========" << std::endl;
        //   std::cout << opp << std::endl;
        arma::mat k;
        ReadFileIntoMatrix(opp, k);

        std::string opp2 = p2.c_str();
        std::cout << p2.c_str() << std::endl;
        // std::vector<std::string> strings;
        // customSplit(opp, '/', strings);
        // std::string finall = strings[strings.size()-1];
        // out << "Finall (standard): " << finall << std::endl;
        // out << "Read path (standard): " << opp<< std::endl;
        //  out << "========" << std::endl;
        //   std::cout << opp << std::endl;
        arma::mat k2;
        ReadFileIntoMatrix(opp2, k2);


        typedef BottleneckCalculator<Bipartie_Graph,GraphTreeController> bipartiebd;
        typedef bd_general_interface<bipartiebd, bipartiebd> naiveBdCalc;

        EuclideanBDCompute<bipartiebd> eucBd;
        Nested_Updated<PBI_VD<naiveBdCalc>, arma::mat> nestedNaive;
        //!  Nested_Base<SCD, naiveBdCalc, std::pair<arma::mat, std::vector<std::vector<bool>>>> nestedNaive;
        //!    Nested_Base<PBI, naiveBdCalc, arma::mat> nestedNaive;

        long totalTimeForSMIstd;
        long totalTimeForBD;
        std::vector<int> matching(k.n_cols);
        std::vector<std::vector<int>> res;

        std::cout << "Starting the boy" << std::endl;


        double bd = nestedNaive.Compute_Isometry(k, k2, matching, totalTimeForSMIstd, totalTimeForBD);

        std::cout << "method exit succesful" << std::endl;

        int a = 3-5;
        std::cout << a << std::endl;

        std::cout << "Bottleneck d computed as : " << bd << "1234323" << std::endl;
        std::cout << "Resulting into: " << std::endl;
        std::cout << "Size of matching: " << matching.size() << std::endl;
        for (int i = 0; i < matching.size() ; i++)
        {
            std::cout << matching[i] << ", ";
        }
        std::cout << "Final " << std::endl;


    }

    void runSingelton(int u, int v )
    {
        //LayeredSCD lap(0.43,2,"");

        //  std::string debugStream = "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/GraphTypeInvariants/EMDPDD/debug.txt";
        typedef BottleneckCalculator<Bipartie_Graph,GraphTreeController> bipartiebd;
        typedef bd_general_interface<bipartiebd, bipartiebd> naiveBdCalc;

        EuclideanBDCompute<bipartiebd> eucBd;
        Nested_Updated<PBI_Update<naiveBdCalc>, arma::mat> nestedNaive;
        //!  Nested_Base<SCD, naiveBdCalc, std::pair<arma::mat, std::vector<std::vector<bool>>>> nestedNaive;
        //!    Nested_Base<PBI, naiveBdCalc, arma::mat> nestedNaive;

        long totalTimeForSMIstd;
        long totalTimeForBD;
        std::cout << "Size of matrices[u] :" << matrices[u].n_cols << std::endl;
        std::cout << "Size of matrices[v] :" << matrices[v].n_cols << std::endl;
        std::vector<int> matching(matrices[u].n_cols);
        std::vector<std::vector<int>> res;

        double bd = nestedNaive.Compute_Isometry(matrices[u], matrices[v], matching, totalTimeForSMIstd, totalTimeForBD);
        std::cout << "Bottleneck d computed as : " << bd << std::endl;
        std::cout << "Resulting into: " << std::endl;

        std::cout << "We also have " << matching.size() << std::endl;

        for (int i = 0; i < matching.size() ; i++)
        {

            std::cout << matching[i] << ", ";

        }




        // std::cout << "Outvector size: " << outVector.size() << std::endl;


    }



    void limitedExperRunFive(int x, int y,std::vector<std::vector<std::string>> & parsedCsv, std::vector<std::vector<std::string>> & outVector )
    {
        //LayeredSCD lap(0.43,2,"");
        int tt = matrices.size();
        y = std::min(y, tt);
        //  std::string debugStream = "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/GraphTypeInvariants/EMDPDD/debug.txt";
        typedef BottleneckCalculator<Bipartie_Graph,GraphTreeController> bipartiebd;
        typedef bd_general_interface<bipartiebd, bipartiebd> naiveBdCalc;

        EuclideanBDCompute<bipartiebd> eucBd;
        Nested_Updated<PBI_VD<naiveBdCalc>, arma::mat> nestedNaive;
        //!  Nested_Base<SCD, naiveBdCalc, std::pair<arma::mat, std::vector<std::vector<bool>>>> nestedNaive;
        //!    Nested_Base<PBI, naiveBdCalc, arma::mat> nestedNaive;


        for (int i = x; i < y; i++)
        {
            std::cout << "Running operation (Full Graph): [" << x << " --- " << i << " --- " << y << "]" << std::endl;
            std::vector<std::string> tmp;
            int u = stoi(parsedCsv[i][1]);
            int v = stoi(parsedCsv[i][2]);
            long totalTimeForSMIstd;
            long totalTimeForBD;
            std::vector<int> matching;
            std::vector<std::vector<int>> res;
            double bd = 0;

            //!     auto firstP = std::make_pair(matrices[u], graphs[u]);
            //!    auto secondP = std::make_pair( matrices[v], graphs[v]);

            // std::cout << "Clouds are : " << std::endl;
            //  std::cout << "Dimension of A: " << matrices[u].n_cols << ", " << matrices[u].n_rows << std::endl;
            //   std::cout << "Dimension of B: " << matrices[v].n_cols << ", " << matrices[v].n_rows << std::endl;



            //!    bd = nestedNaive.Compute_Isometry(firstP, secondP, matching, totalTimeForSMIstd, totalTimeForBD);

            bd = nestedNaive.Compute_Isometry(matrices[u], matrices[v], matching, totalTimeForSMIstd, totalTimeForBD);

            //std::cout << "Computed isometry " << std::endl;

            double bdtruncate = -1;

            long totalTimeForSMIstdfake;
            long totalTimeForBDfake;


            std::vector<int> matchingOne;
            double purebd = eucBd.compute(matrices[u], matrices[v], matchingOne);
            //!CrystData::computeEuclideanBottleneckDistanceBC(matrices[u], matrices[v]);


            double haussdorff = computeHaudsorff(matrices[u], matrices[v]);


            //! okay boomer:

            arma::mat AA;

            ExtraExperiments::RearrangeCloud(matrices[u], AA, matching);

            arma::mat Rv;
            arma::vec tttv;

            ExtraExperiments::findRotation(AA, matrices[v],  Rv, tttv);

            arma::mat AAA = Rv*AA;

            ExtraExperiments::ShiftBy(AAA,tttv);

            double realDist = 0;

            for (int i = 0; i < matching.size(); i++)
            {

                double distP = arma::norm(AAA.col(i)-matrices[v].col(i), 2);

                realDist = std::max(realDist, distP);


            }




            double ratiodistbd = realDist / bd;
            bdtruncate = ratiodistbd;



            //! to_string_with_precision<double>(haussdorff , 15)
            tmp.push_back(parsedCsv[i][1]);
            tmp.push_back(parsedCsv[i][2]);
            tmp.push_back(parsedCsv[i][3]);
            tmp.push_back(to_string_with_precision<double>(haussdorff, 15));
            tmp.push_back(to_string_with_precision<double>(purebd, 15));
            tmp.push_back(parsedCsv[i][4]);
            tmp.push_back(parsedCsv[i][0]);
            tmp.push_back(to_string_with_precision<double>(bd, 15));
            double ratiod = bd / purebd;
            tmp.push_back(to_string_with_precision<double>(ratiod, 15));
            tmp.push_back(std::to_string(totalTimeForSMIstd));
            tmp.push_back(std::to_string(totalTimeForBD));
            tmp.push_back(std::to_string(bdtruncate));



            outVector.push_back(tmp);
        }

        // std::cout << "Outvector size: " << outVector.size() << std::endl;


    }

    void limitedExperRun(int x, int y,std::vector<std::vector<std::string>> & parsedCsv, std::vector<std::vector<std::string>> & outVector )
    {
        //LayeredSCD lap(0.43,2,"");
        int tt = matrices.size();
        y = std::min(y, tt);
        //  std::string debugStream = "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/GraphTypeInvariants/EMDPDD/debug.txt";
        typedef BottleneckCalculator<Bipartie_Graph,GraphTreeController> bipartiebd;
        typedef bd_general_interface<bipartiebd, bipartiebd> naiveBdCalc;

        EuclideanBDCompute<bipartiebd> eucBd;
        //!  Nested_Base<SCD, naiveBdCalc, std::pair<arma::mat, std::vector<std::vector<bool>>>> nestedNaive;
        Nested_Base<PBI, naiveBdCalc, arma::mat> nestedNaive;


        for (int i = x; i < y; i++)
        {
            std::cout << "Running operation (Full Graph): [" << x << " --- " << i << " --- " << y << "]" << std::endl;
            std::vector<std::string> tmp;
            int u = stoi(parsedCsv[i][1]);
            int v = stoi(parsedCsv[i][2]);
            long totalTimeForSMIstd;
            long totalTimeForBD;
            std::vector<int> matching;
            std::vector<std::vector<int>> res;
            double bd = 0;

            //!     auto firstP = std::make_pair(matrices[u], graphs[u]);
            //!    auto secondP = std::make_pair( matrices[v], graphs[v]);

            // std::cout << "Clouds are : " << std::endl;
            //  std::cout << "Dimension of A: " << matrices[u].n_cols << ", " << matrices[u].n_rows << std::endl;
            //   std::cout << "Dimension of B: " << matrices[v].n_cols << ", " << matrices[v].n_rows << std::endl;



            //!    bd = nestedNaive.Compute_Isometry(firstP, secondP, matching, totalTimeForSMIstd, totalTimeForBD);

            bd = nestedNaive.Compute_Isometry(matrices[u], matrices[v], matching, totalTimeForSMIstd, totalTimeForBD);

            //std::cout << "Computed isometry " << std::endl;

            double bdtruncate = -1;

            long totalTimeForSMIstdfake;
            long totalTimeForBDfake;


            std::vector<int> matchingOne;
            double purebd = eucBd.compute(matrices[u], matrices[v], matchingOne);
            //!CrystData::computeEuclideanBottleneckDistanceBC(matrices[u], matrices[v]);


            double haussdorff = computeHaudsorff(matrices[u], matrices[v]);

            //! to_string_with_precision<double>(haussdorff , 15)
            tmp.push_back(parsedCsv[i][1]);
            tmp.push_back(parsedCsv[i][2]);
            tmp.push_back(parsedCsv[i][3]);
            tmp.push_back(to_string_with_precision<double>(haussdorff, 15));
            tmp.push_back(to_string_with_precision<double>(purebd, 15));
            tmp.push_back(parsedCsv[i][4]);
            tmp.push_back(parsedCsv[i][0]);
            tmp.push_back(to_string_with_precision<double>(bd, 15));
            double ratiod = bd / purebd;
            tmp.push_back(to_string_with_precision<double>(ratiod, 15));
            tmp.push_back(std::to_string(totalTimeForSMIstd));
            tmp.push_back(std::to_string(totalTimeForBD));
            tmp.push_back(std::to_string(bdtruncate));



            outVector.push_back(tmp);
        }

        // std::cout << "Outvector size: " << outVector.size() << std::endl;


    }

    void parseCSV(std::vector<std::vector<std::string>> & parsedCsv, std::string input)
    {
        std::ifstream  data(input);
        std::string line;
        while(std::getline(data,line))
        {
            std::stringstream lineStream(line);
            std::string cell;
            std::vector<std::string> parsedRow;
            while(std::getline(lineStream,cell,','))
            {
                parsedRow.push_back(cell);
            }

            parsedCsv.push_back(parsedRow);
        }

    }

    double computeHaudsorff(arma::mat & A,arma::mat & B)
    {
        double maxxOneWayReal = 0;
        for (int i = 0; i < A.n_cols; i++)
        {
            double maxxOneWay = std::numeric_limits<double>::max();
            for (int j = 0; j < B.n_cols; j++)
            {
                double dtmp = arma::norm(A.col(i) - B.col(j), 2);
                maxxOneWay = std::min(maxxOneWay, dtmp);

            }
            maxxOneWayReal  = std::max(maxxOneWayReal,maxxOneWay );

        }
        double maxxTwoWayReal = 0;
        for (int j = 0; j < B.n_cols; j++)
        {
            double maxxTwoWay = std::numeric_limits<double>::max();
            for (int i = 0; i < A.n_cols; i++)
            {
                double dtmp = arma::norm(A.col(i) - B.col(j), 2);
                maxxTwoWay = std::min(maxxTwoWay, dtmp);

            }
            maxxTwoWayReal  = std::max(maxxTwoWayReal,maxxTwoWay);



        }


        return std::max(maxxTwoWayReal,maxxOneWayReal);
    }

    void printExperiments(std::vector<std::vector<std::string>> & values, std::string header, std::string output)
    {

        std::ofstream out(output);
        out << header << std::endl;
        for (int i = 0; i < values.size(); i++)
        {
            for (int j = 0; j < values[i].size(); j++)
            {
                out << values[i][j];
                if (j  < values[i].size() - 1)
                {
                    out << ",";
                }

            }
            out << std::endl;

        }
        out.close();

    }


    void RunThisParallelSCDUpdate(std::string inputCsv, std::string outputCsv, int threads, int maxxnumm)
    {
        BS::thread_pool pool(threads);
        std::vector<std::vector<std::string>> parsedCSV;
        parseCSV(parsedCSV, inputCsv);
        parsedCSV.erase(parsedCSV.begin());
        parsedCSV.resize(maxxnumm);
        std::vector<std::vector<std::vector<std::string>>> vectorCollection(maxxnumm);


        //! limitedExperRunTwo

        //! limitedExperRun, this, i*be, (i+1)*be, std::ref(parsedCSV), std::ref(vectorCollection[i])

        for (int i = 0; i < maxxnumm; i++)
        {

            //!  limitedExperRun(i, (i+1), parsedCSV, vectorCollection[i]);
            pool.push_task(&ExperimentRunner::limitedExperRunTwo, this, i, (i+1), std::ref(parsedCSV), std::ref(vectorCollection[i]));
            //  pool.push_task([i, &parsedCSV, &vectorCollection[i]]
            //   {
            //      this->limitedExperRun(i, i+1,parsedCSV)
            //   }
            //  );



        }

        pool.wait_for_tasks();

        std::string header ="Molecule_A, Molecule_B, NumAtoms, Hausdorff distance, BD, L_inf on SDV, EMD on PDD, NBD on NCD,  NBD_NCD / BD, Invariant build time (ms), NBD time (ms), ExtraColumn";
        header = header + ", SCDUpd, SCDUpd / BD, Invariant build time (ms), SCDUpd time (ms)";
        //! std::string header ="Bottleneck dist, Time SMI, Time Bottleneck, PDD dist, Molecule_A, Molecule_B, NumAtoms, SDV dist";
        std::vector<std::vector<std::string>> finalString;

        for (int i = 0; i < vectorCollection.size(); i++)
        {
            finalString.insert(finalString.end(), vectorCollection[i].begin(), vectorCollection[i].end());

        }

        std::cout << "Total final string size " << finalString.size() << std::endl;

        std::vector<std::vector<std::string>> finalStringLimited = finalString;

        std::sort(finalString.begin(), finalString.end(),
                  [](const std::vector<std::string> & a, const std::vector<std::string> & b) -> bool
        {
            double aa = stod(a[7]);
            double bb = stod(b[7]);
            return  aa < bb;
        });

        long sumsd = 0;
        long sumbd = 0;
        long avgsd =0;
        long avgbd = 0;



        for (int i = 0; i < finalString.size(); i++)
        {
            sumsd = sumsd + stol(finalString[i][1]);
            sumbd = sumbd + stol(finalString[i][2]);
        }

        avgsd = sumsd / (double)finalString.size();
        avgbd = sumbd / (double)finalString.size();

        printExperiments(finalString, header, outputCsv);

        std::cout << "Total time SMI: " << sumsd << "milliseconds" << std::endl;
        std::cout << "Total time BD: " << sumbd << "milliseconds" << std::endl;
        std::cout << "Avg time SMI: " << avgsd<< "milliseconds" << std::endl;
        std::cout << "Avg time BD: " << avgbd << "milliseconds" << std::endl;


        // finalString;
        //finalStringLimited.resize(500);

        //!  printExperiments(finalStringLimited, header, outputCsvLimited);

//   double ukk = parsedCSV.size() / (double) threads;
        // int be = floor(ukk);
        // int iblan = 0;




    }

    void RunThisThingExperimentsMay(std::string inputCsv, std::string outputCsv,  int threads, int start, int endd)
    {
        BS::thread_pool pool(threads);
        std::vector<std::vector<std::string>> parsedCSV;
        parseCSV(parsedCSV, inputCsv);
        parsedCSV.erase(parsedCSV.begin());
        std::vector<std::vector<std::vector<std::string>>> vectorCollection(endd - start);


        //    parsedCSV.resize(maxxnumm);

        for (int i = start; i < endd; i++)
        {

            //!  limitedExperRun(i, (i+1), parsedCSV, vectorCollection[i]);
            pool.push_task(&ExperimentRunner::limitedExperThirteen, this, i, std::ref(parsedCSV), std::ref(vectorCollection[i-start]));
            //  pool.push_task([i, &parsedCSV, &vectorCollection[i]]
            //   {
            //      this->limitedExperRun(i, i+1,parsedCSV)
            //   }
            //  );



        }


        pool.wait_for_tasks();

        std::string header ="Molecule_A, Molecule_B,NBD";
        //! std::string header ="Bottleneck dist, Time SMI, Time Bottleneck, PDD dist, Molecule_A, Molecule_B, NumAtoms, SDV dist";
        std::vector<std::vector<std::string>> finalString;

        for (int i = 0; i < vectorCollection.size(); i++)
        {
            finalString.insert(finalString.end(), vectorCollection[i].begin(), vectorCollection[i].end());

        }




        printExperiments(finalString, header, outputCsv);

//        std::cout << "Total time SMI: " << sumsd << "milliseconds" << std::endl;
//       std::cout << "Total time BD: " << sumbd << "milliseconds" << std::endl;
        //      std::cout << "Avg time SMI: " << avgsd<< "milliseconds" << std::endl;
        //    std::cout << "Avg time BD: " << avgbd << "milliseconds" << std::endl;










    }

    void RunThisParallelSCD(std::string inputCsv, std::string outputCsv, std::string outputCsvLimited, int threads, int maxxnumm, int minnumm = 0)
    {
        BS::thread_pool pool(threads);
        std::vector<std::vector<std::string>> parsedCSV;
        parseCSV(parsedCSV, inputCsv);
        parsedCSV.erase(parsedCSV.begin());
        parsedCSV.resize(maxxnumm);
        std::vector<std::vector<std::vector<std::string>>> vectorCollection(maxxnumm-minnumm);


        //! limitedExperRun, this, i*be, (i+1)*be, std::ref(parsedCSV), std::ref(vectorCollection[i])

        for (int i = minnumm; i < maxxnumm; i++)
        {

            //!  limitedExperRun(i, (i+1), parsedCSV, vectorCollection[i]);
            pool.push_task(&ExperimentRunner::limitedExperRunFive, this, i, (i+1), std::ref(parsedCSV), std::ref(vectorCollection[i-minnumm]));
            //  pool.push_task([i, &parsedCSV, &vectorCollection[i]]
            //   {
            //      this->limitedExperRun(i, i+1,parsedCSV)
            //   }
            //  );



        }

        pool.wait_for_tasks();

        std::string header ="Molecule_A, Molecule_B, NumAtoms, Hausdorff distance, BD, L_inf on SDV, EMD on PDD, NBD on NCD,  NBD_NCD / BD, Invariant build time (ms), NBD time (ms), ExtraColumn";
        //! std::string header ="Bottleneck dist, Time SMI, Time Bottleneck, PDD dist, Molecule_A, Molecule_B, NumAtoms, SDV dist";
        std::vector<std::vector<std::string>> finalString;

        for (int i = 0; i < vectorCollection.size(); i++)
        {
            finalString.insert(finalString.end(), vectorCollection[i].begin(), vectorCollection[i].end());

        }

        std::cout << "Total final string size " << finalString.size() << std::endl;

        std::vector<std::vector<std::string>> finalStringLimited = finalString;

        std::sort(finalString.begin(), finalString.end(),
                  [](const std::vector<std::string> & a, const std::vector<std::string> & b) -> bool
        {
            double aa = stod(a[7]);
            double bb = stod(b[7]);
            return  aa < bb;
        });

        long sumsd = 0;
        long sumbd = 0;
        long avgsd =0;
        long avgbd = 0;



        for (int i = 0; i < finalString.size(); i++)
        {
            sumsd = sumsd + stol(finalString[i][1]);
            sumbd = sumbd + stol(finalString[i][2]);
        }

        avgsd = sumsd / (double)finalString.size();
        avgbd = sumbd / (double)finalString.size();

        printExperiments(finalString, header, outputCsv);

        std::cout << "Total time SMI: " << sumsd << "milliseconds" << std::endl;
        std::cout << "Total time BD: " << sumbd << "milliseconds" << std::endl;
        std::cout << "Avg time SMI: " << avgsd<< "milliseconds" << std::endl;
        std::cout << "Avg time BD: " << avgbd << "milliseconds" << std::endl;


        // finalString;
        //finalStringLimited.resize(500);

        printExperiments(finalStringLimited, header, outputCsvLimited);

//   double ukk = parsedCSV.size() / (double) threads;
        // int be = floor(ukk);
        // int iblan = 0;




    }

    void RunThisParallelSCD2D(std::string inputCsv, std::string outputCsv, std::string outputCsvLimited, int threads, int maxxnumm)
    {
        BS::thread_pool pool(threads);
        std::vector<std::vector<std::string>> parsedCSV;
        parseCSV(parsedCSV, inputCsv);
        parsedCSV.erase(parsedCSV.begin());
        parsedCSV.resize(maxxnumm);
        std::vector<std::vector<std::vector<std::string>>> vectorCollection(maxxnumm);


        //! limitedExperRun, this, i*be, (i+1)*be, std::ref(parsedCSV), std::ref(vectorCollection[i])

        for (int i = 0; i < maxxnumm; i++)
        {

            //!  limitedExperRun(i, (i+1), parsedCSV, vectorCollection[i]);
            pool.push_task(&ExperimentRunner::limitedExperRunSix, this, i, (i+1), std::ref(parsedCSV), std::ref(vectorCollection[i]));
            //  pool.push_task([i, &parsedCSV, &vectorCollection[i]]
            //   {
            //      this->limitedExperRun(i, i+1,parsedCSV)
            //   }
            //  );



        }

        pool.wait_for_tasks();

        std::string header ="Molecule_A,Molecule_B,Num_Atoms,Hausdorff_distance,BD,L_inf on SRV,L_inf on SDV,EMD on PDD,NBM_NDP,NBM_NDP/BD,time_NDP(microseconds),NBD_time(microseconds),Match_dist/NBM_NDP";
        //! std::string header ="Bottleneck dist, Time SMI, Time Bottleneck, PDD dist, Molecule_A, Molecule_B, NumAtoms, SDV dist";
        std::vector<std::vector<std::string>> finalString;

        for (int i = 0; i < vectorCollection.size(); i++)
        {
            finalString.insert(finalString.end(), vectorCollection[i].begin(), vectorCollection[i].end());

        }

        std::cout << "Total final string size " << finalString.size() << std::endl;

        std::vector<std::vector<std::string>> finalStringLimited = finalString;

        std::sort(finalString.begin(), finalString.end(),
                  [](const std::vector<std::string> & a, const std::vector<std::string> & b) -> bool
        {
            double aa = stod(a[8]);
            double bb = stod(b[8]);
            return  aa < bb;
        });

        long sumsd = 0;
        long sumbd = 0;
        long avgsd =0;
        long avgbd = 0;



        for (int i = 0; i < finalString.size(); i++)
        {
            sumsd = sumsd + stol(finalString[i][1]);
            sumbd = sumbd + stol(finalString[i][2]);
        }

        avgsd = sumsd / (double)finalString.size();
        avgbd = sumbd / (double)finalString.size();

        printExperiments(finalString, header, outputCsv);

        std::cout << "Total time SMI: " << sumsd << "milliseconds" << std::endl;
        std::cout << "Total time BD: " << sumbd << "milliseconds" << std::endl;
        std::cout << "Avg time SMI: " << avgsd<< "milliseconds" << std::endl;
        std::cout << "Avg time BD: " << avgbd << "milliseconds" << std::endl;


        // finalString;
        //finalStringLimited.resize(500);

        printExperiments(finalStringLimited, header, outputCsvLimited);

//   double ukk = parsedCSV.size() / (double) threads;
        // int be = floor(ukk);
        // int iblan = 0;




    }

    void RunThisParallelSCD3DHard(std::string inputCsv, std::string outputCsv, std::string outputCsvLimited, int threads, int maxxnumm)
    {
        BS::thread_pool pool(threads);
        std::vector<std::vector<std::string>> parsedCSV;
        parseCSV(parsedCSV, inputCsv);
        parsedCSV.erase(parsedCSV.begin());
        parsedCSV.resize(maxxnumm);
        std::vector<std::vector<std::vector<std::string>>> vectorCollection(maxxnumm);


        //! limitedExperRun, this, i*be, (i+1)*be, std::ref(parsedCSV), std::ref(vectorCollection[i])

        for (int i = 0; i < maxxnumm; i++)
        {

            //!  limitedExperRun(i, (i+1), parsedCSV, vectorCollection[i]);
            pool.push_task(&ExperimentRunner::limitedExperRunSeven, this, i, (i+1), std::ref(parsedCSV), std::ref(vectorCollection[i]));
            //  pool.push_task([i, &parsedCSV, &vectorCollection[i]]
            //   {
            //      this->limitedExperRun(i, i+1,parsedCSV)
            //   }
            //  );



        }

        pool.wait_for_tasks();

        std::string header ="Molecule_A,Molecule_B,Num_Atoms,Hausdorff_distance,BD,L_inf on SRV, L_inf on SDV,EMD on PDD,NBM_NCD2, NBM_NCD1 ,NBM_NDP/BD,time_NDP(microseconds),NBD_time(microseconds),Match_dist/NBM_NDP";
        //! std::string header ="Bottleneck dist, Time SMI, Time Bottleneck, PDD dist, Molecule_A, Molecule_B, NumAtoms, SDV dist";
        std::vector<std::vector<std::string>> finalString;

        for (int i = 0; i < vectorCollection.size(); i++)
        {
            finalString.insert(finalString.end(), vectorCollection[i].begin(), vectorCollection[i].end());

        }

        std::cout << "Total final string size " << finalString.size() << std::endl;

        std::vector<std::vector<std::string>> finalStringLimited = finalString;

        std::sort(finalString.begin(), finalString.end(),
                  [](const std::vector<std::string> & a, const std::vector<std::string> & b) -> bool
        {
            double aa = stod(a[8]);
            double bb = stod(b[8]);
            return  aa < bb;
        });

        long sumsd = 0;
        long sumbd = 0;
        long avgsd =0;
        long avgbd = 0;



        for (int i = 0; i < finalString.size(); i++)
        {
            sumsd = sumsd + stol(finalString[i][1]);
            sumbd = sumbd + stol(finalString[i][2]);
        }

        avgsd = sumsd / (double)finalString.size();
        avgbd = sumbd / (double)finalString.size();

        printExperiments(finalString, header, outputCsv);

        std::cout << "Total time SMI: " << sumsd << "milliseconds" << std::endl;
        std::cout << "Total time BD: " << sumbd << "milliseconds" << std::endl;
        std::cout << "Avg time SMI: " << avgsd<< "milliseconds" << std::endl;
        std::cout << "Avg time BD: " << avgbd << "milliseconds" << std::endl;


        // finalString;
        //finalStringLimited.resize(500);

        printExperiments(finalStringLimited, header, outputCsvLimited);

//   double ukk = parsedCSV.size() / (double) threads;
        // int be = floor(ukk);
        // int iblan = 0;




    }

    void RunThisParallelSCD3DHardVerTwo(std::string inputCsv, std::string outputCsv, std::string outputCsvLimited, int threads, int maxxnumm)
    {
        BS::thread_pool pool(threads);
        std::vector<std::vector<std::string>> parsedCSV;
        parseCSV(parsedCSV, inputCsv);
        parsedCSV.erase(parsedCSV.begin());
        parsedCSV.resize(maxxnumm);
        std::vector<std::vector<std::vector<std::string>>> vectorCollection(maxxnumm);

        std::cout << "MAXNUMM IS : " << maxxnumm << std::endl;


        //! limitedExperRun, this, i*be, (i+1)*be, std::ref(parsedCSV), std::ref(vectorCollection[i])

        for (int i = 0; i < maxxnumm; i++)
        {
            //! limitedExperRunEight

            //!  limitedExperRun(i, (i+1), parsedCSV, vectorCollection[i]);
            // ExperimentRunner::limitedExperRunTwelve( i, (i+1), std::ref(parsedCSV), std::ref(vectorCollection[i]));
            pool.push_task(&ExperimentRunner::limitedExperRunTwelve, this, i, (i+1), std::ref(parsedCSV), std::ref(vectorCollection[i]));
            //  pool.push_task([i, &parsedCSV, &vectorCollection[i]]
            //   {
            //      this->limitedExperRun(i, i+1,parsedCSV)
            //   }
            //  );



        }

        pool.wait_for_tasks();

        //!std::string header ="Molecule_A,Molecule_B,Num_Atoms,Hausdorff_distance,BD,L_inf on SRV, L_inf on SDV,EMD on PDD,NBM_NCD2, NBM_NCD1 ,NBM_NDP/BD,time_NDP(microseconds),NBD_time(microseconds),Match_dist/NBM_NDP";
        std::string header ="Molecule_A,Molecule_B,Num_Atoms,Hausdorff_distance,BD,L_inf on SRV, L_inf on SDV,EMD on PDD,NBM_NCD2, NBM_NCD1 ,NBM_NDP/BD,time_NCD2(microseconds),NBD_NCD2_time(microseconds),time_NCD1(microseconds),NBD_NCD1_time(microseconds)";

        //! std::string header ="Bottleneck dist, Time SMI, Time Bottleneck, PDD dist, Molecule_A, Molecule_B, NumAtoms, SDV dist";
        std::vector<std::vector<std::string>> finalString;

        for (int i = 0; i < vectorCollection.size(); i++)
        {
            finalString.insert(finalString.end(), vectorCollection[i].begin(), vectorCollection[i].end());

        }

        std::cout << "Total final string size " << finalString.size() << std::endl;

        std::vector<std::vector<std::string>> finalStringLimited = finalString;

        std::sort(finalString.begin(), finalString.end(),
                  [](const std::vector<std::string> & a, const std::vector<std::string> & b) -> bool
        {
            double aa = stod(a[8]);
            double bb = stod(b[8]);
            return  aa < bb;
        });

        long sumsd = 0;
        long sumbd = 0;
        long avgsd =0;
        long avgbd = 0;



        for (int i = 0; i < finalString.size(); i++)
        {
            sumsd = sumsd + stol(finalString[i][1]);
            sumbd = sumbd + stol(finalString[i][2]);
        }

        avgsd = sumsd / (double)finalString.size();
        avgbd = sumbd / (double)finalString.size();

        printExperiments(finalString, header, outputCsv);

        std::cout << "Total time SMI: " << sumsd << "milliseconds" << std::endl;
        std::cout << "Total time BD: " << sumbd << "milliseconds" << std::endl;
        std::cout << "Avg time SMI: " << avgsd<< "milliseconds" << std::endl;
        std::cout << "Avg time BD: " << avgbd << "milliseconds" << std::endl;


        // finalString;
        //finalStringLimited.resize(500);

        printExperiments(finalStringLimited, header, outputCsvLimited);

//   double ukk = parsedCSV.size() / (double) threads;
        // int be = floor(ukk);
        // int iblan = 0;




    }

    void RunThisParallelSCD3DHardVerFour(std::vector<std::vector<int>> & inputPairs, int threads, std::string outputfile)
    {
        BS::thread_pool pool(threads);
        std::vector<std::vector<std::string>> answers(inputPairs.size());
        for (int i = 0; i < inputPairs.size(); i++)
        {
            pool.push_task(&ExperimentRunner::SingleRunSpezzialeVerThree, this, inputPairs[i][0], inputPairs[i][1], std::ref(answers), i);

        }
        pool.wait_for_tasks();
        std::string header = "molecule_A, molecule_B, scd3";
        printExperiments(answers, header, outputfile);




    }


    void RunThisParallelSCD3DHardVerThree(std::string inputCsv, std::string outputCsv, std::string outputCsvLimited, int threads, int maxxnumm)
    {
        BS::thread_pool pool(threads);
        std::vector<std::vector<std::string>> parsedCSV;
        //std::cout << "Preparing to read" << std::endl;
        parseCSV(parsedCSV, inputCsv);
        // std::cout << "Reading complete" << std::endl;
        parsedCSV.erase(parsedCSV.begin());
        parsedCSV.resize(maxxnumm);
        std::vector<std::vector<std::vector<std::string>>> vectorCollection(maxxnumm);


        //! limitedExperRun, this, i*be, (i+1)*be, std::ref(parsedCSV), std::ref(vectorCollection[i])

        for (int i = 0; i < maxxnumm; i++)
        {

            //!  limitedExperRun(i, (i+1), parsedCSV, vectorCollection[i]);
            pool.push_task(&ExperimentRunner::limitedExperRunNine, this, i, (i+1), std::ref(parsedCSV), std::ref(vectorCollection[i]));
            //  pool.push_task([i, &parsedCSV, &vectorCollection[i]]
            //   {
            //      this->limitedExperRun(i, i+1,parsedCSV)
            //   }
            //  );



        }

        pool.wait_for_tasks();

        //!std::string header ="Molecule_A,Molecule_B,Num_Atoms,Hausdorff_distance,BD,L_inf on SRV, L_inf on SDV,EMD on PDD,NBM_NCD2, NBM_NCD1 ,NBM_NDP/BD,time_NDP(microseconds),NBD_time(microseconds),Match_dist/NBM_NDP";
        std::string header ="Molecule_A,Molecule_B,Num_Atoms,Hausdorff_distance,BD,L_inf on SRV, L_inf on SDV,EMD on PDD,NBM_NCD2, NBM_NCD1, NBM_NCD3, NBM_NCD3 / BD,NBM_NCD3/BD,time_NCD2(microseconds),NBD_NCD2_time(microseconds),time_NCD1(microseconds),NBD_NCD1_time(microseconds),time_NCD3,NCD3_time";

        //! std::string header ="Bottleneck dist, Time SMI, Time Bottleneck, PDD dist, Molecule_A, Molecule_B, NumAtoms, SDV dist";
        std::vector<std::vector<std::string>> finalString;

        for (int i = 0; i < vectorCollection.size(); i++)
        {
            finalString.insert(finalString.end(), vectorCollection[i].begin(), vectorCollection[i].end());

        }

        std::cout << "Total final string size " << finalString.size() << std::endl;

        std::vector<std::vector<std::string>> finalStringLimited = finalString;

        std::sort(finalString.begin(), finalString.end(),
                  [](const std::vector<std::string> & a, const std::vector<std::string> & b) -> bool
        {
            double aa = stod(a[10]);
            double bb = stod(b[10]);
            return  aa < bb;
        });

        long sumsd = 0;
        long sumbd = 0;
        long avgsd =0;
        long avgbd = 0;



        for (int i = 0; i < finalString.size(); i++)
        {
            sumsd = sumsd + stol(finalString[i][1]);
            sumbd = sumbd + stol(finalString[i][2]);
        }

        avgsd = sumsd / (double)finalString.size();
        avgbd = sumbd / (double)finalString.size();

        printExperiments(finalString, header, outputCsv);

        std::cout << "Total time SMI: " << sumsd << "milliseconds" << std::endl;
        std::cout << "Total time BD: " << sumbd << "milliseconds" << std::endl;
        std::cout << "Avg time SMI: " << avgsd<< "milliseconds" << std::endl;
        std::cout << "Avg time BD: " << avgbd << "milliseconds" << std::endl;


        // finalString;
        //finalStringLimited.resize(500);

        printExperiments(finalStringLimited, header, outputCsvLimited);

//   double ukk = parsedCSV.size() / (double) threads;
        // int be = floor(ukk);
        // int iblan = 0;




    }





    virtual ~ExperimentRunner() {}

protected:

private:
    std::vector<std::string> cloudPaths;
    std::vector<std::string> graphPaths;
    std::vector<arma::mat> matrices;
    std::vector<std::vector<std::vector<bool>>> graphs;
};

#endif // EXPERIMENTRUNNER_H
