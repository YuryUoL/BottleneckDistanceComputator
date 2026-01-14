#ifndef GRAPHEXPERIMENTS_H
#define GRAPHEXPERIMENTS_H

#include "header.hpp"

#include "CloudIsometry/Invariants/NCD_Smart.h"
#include "CloudIsometry/Invariants/NCD_Graph.h"
#include "CloudIsometry/Invariants/NCD_Smart_Update.h"
#include "CloudIsometry/Invariants/PBI_Update.h"
#include "CloudIsometry/Nested_Updated.h"
#include "Bottleneckdistance/Misc/GraphTreeControllerVector.h"
#include <set>

class GraphExperiments
{
public:
    GraphExperiments() {}

    void recursiveVertexAdition(int minVal, std::vector<int> & setPotential, Graph & in, int i, std::unordered_set<int> setA, Graph & out, std::vector<std::vector<int>> & MetaSets,
                                std::vector<std::unordered_set<int>> & CorrespondingSets, std::vector<int> & vertexIndxTracker)
    {
        std::cout << i << " - " << minVal << " : [ ";

        for (int jj : setA)
        {
            std::cout  << jj << ", ";

        }
        std::cout << "]" << std::endl;

        if (setA.size() % 2 == 0)
        {
            int ni = boost::add_vertex(out);
            vertexIndxTracker.push_back(i);
            CorrespondingSets.push_back(setA);
            MetaSets[i].push_back(ni);

        }
        if (setA.size() == setPotential.size())
        {
            return;
        }
        for (int jjj = minVal ; jjj < setPotential.size() ; jjj++)
        {
            std::unordered_set<int> setB = setA;
            setB.insert(setPotential[jjj]);
            recursiveVertexAdition(jjj+1, setPotential, in, i, setB, out, MetaSets, CorrespondingSets, vertexIndxTracker);
        }

    }

    void invertGraph(Graph & in, Graph::edge_descriptor & e, Graph & out, int numVertices,  std::vector<std::vector<int>> & MetaSets,
                     std::vector<std::unordered_set<int>> & CorrespondingSets, std::vector<int> & vertexIndxTracker)
    {
        std::cout << "Source: " << boost::source(e,in) << " - Target: "  << boost::target(e,in) << std::endl;
        // std::cout << boost::target(e,in) << std::endl;
        for (int i = 0; i < numVertices; i++)
        {
            boost::add_vertex(out);
        }

        std::vector<std::string> lulz;

        for (int i = 0; i < boost::num_vertices(out); i++)
        {
            std::string wal = "";
            //   std::cout << i << " : " << vertexIndxTracker[i]+1 << " : [" ;
            wal = wal + std::to_string(vertexIndxTracker[i]+1) + " : [";
            for (int jj : CorrespondingSets[i])
            {
                //     std::cout  << jj+1 << ", ";
                wal = wal + std::to_string(jj+1) + ", ";

            }
            // std::cout << "]" << std::endl;
            wal = wal + "]";
            lulz.push_back(wal);
        }


        for (auto eit = boost::edges(in).first; eit != boost::edges(in).second; eit++)
        {

            int aa = boost::source(*eit, in);
            int bb = boost::target(*eit, in);
            for (int i = 0; i < MetaSets[aa].size(); i++)
            {
                for (int j = 0; j < MetaSets[bb].size(); j++)
                {
                    int aiter = MetaSets[aa][i];
                    int biter = MetaSets[bb][j];
                    if (*eit != e)
                    {
                        bool partOne = ((CorrespondingSets[biter].find(aa) != CorrespondingSets[biter].end()) && (CorrespondingSets[aiter].find(bb) != CorrespondingSets[aiter].end()));
                        bool partTwo = (CorrespondingSets[biter].find(aa) == CorrespondingSets[biter].end()) && (CorrespondingSets[aiter].find(bb) == CorrespondingSets[aiter].end());
                        if (partOne || partTwo)
                        {
                            boost::add_edge(aiter, biter,1.0,out);
                            std::cout << "Adding edge " << "[" << aiter << ", " << biter << "] : " << lulz[aiter] << " --- " << lulz[biter] << std::endl;
                        }
                    }
                    else
                    {
                        bool partOne = ((CorrespondingSets[biter].find(aa) == CorrespondingSets[biter].end()) && (CorrespondingSets[aiter].find(bb) != CorrespondingSets[aiter].end()));
                        bool partTwo = (CorrespondingSets[biter].find(aa) != CorrespondingSets[biter].end()) && (CorrespondingSets[aiter].find(bb) == CorrespondingSets[aiter].end());

                        if (partOne || partTwo)
                        {
                            boost::add_edge(aiter, biter,1.0, out);
                            std::cout << "Adding edge " << "[" << aiter << ", " << biter << "] : " << lulz[aiter] << " --- " << lulz[biter] << std::endl;
                        }


                    }

                }
            }

        }


    }

    void invertGraph(Graph & in, Graph & out, int numVertices,  std::vector<std::vector<int>> & MetaSets,
                     std::vector<std::unordered_set<int>> & CorrespondingSets, std::vector<int> & vertexIndxTracker, std::vector<std::vector<int>> & invertEdges)
    {
        //!    std::cout << "Source: " << boost::source(e,in) << " - Target: "  << boost::target(e,in) << std::endl;
        // std::cout << boost::target(e,in) << std::endl;
        for (int i = 0; i < numVertices; i++)
        {
            boost::add_vertex(out);
        }

        std::vector<std::string> lulz;

        for (int i = 0; i < boost::num_vertices(out); i++)
        {
            std::string wal = "";
            //   std::cout << i << " : " << vertexIndxTracker[i]+1 << " : [" ;
            wal = wal + std::to_string(vertexIndxTracker[i]+1) + " : [";
            for (int jj : CorrespondingSets[i])
            {
                //     std::cout  << jj+1 << ", ";
                wal = wal + std::to_string(jj+1) + ", ";

            }
            // std::cout << "]" << std::endl;
            wal = wal + "]";
            lulz.push_back(wal);
        }

        std::set<Graph::edge_descriptor> edges;

        for (int i = 0 ; i< invertEdges.size(); i++)
        {
            Graph::edge_descriptor e = my_find_edge(invertEdges[i][0], invertEdges[i][1], in);
            edges.insert(e);
        }

        for (auto eit = boost::edges(in).first; eit != boost::edges(in).second; eit++)
        {

            int aa = boost::source(*eit, in);
            int bb = boost::target(*eit, in);
            for (int i = 0; i < MetaSets[aa].size(); i++)
            {
                for (int j = 0; j < MetaSets[bb].size(); j++)
                {
                    int aiter = MetaSets[aa][i];
                    int biter = MetaSets[bb][j];
                    if (edges.find(*eit) == edges.end())
                    {
                        bool partOne = ((CorrespondingSets[biter].find(aa) != CorrespondingSets[biter].end()) && (CorrespondingSets[aiter].find(bb) != CorrespondingSets[aiter].end()));
                        bool partTwo = (CorrespondingSets[biter].find(aa) == CorrespondingSets[biter].end()) && (CorrespondingSets[aiter].find(bb) == CorrespondingSets[aiter].end());
                        if (partOne || partTwo)
                        {
                            boost::add_edge(aiter, biter,1.0,out);
                            std::cout << "Adding edge " << "[" << aiter << ", " << biter << "] : " << lulz[aiter] << " --- " << lulz[biter] << std::endl;
                        }
                    }
                    else
                    {
                        std::cout << "Found bad guy: " << "[" << aa << ", " << bb << "]" << std::endl;
                        bool partOne = ((CorrespondingSets[biter].find(aa) == CorrespondingSets[biter].end()) && (CorrespondingSets[aiter].find(bb) != CorrespondingSets[aiter].end()));
                        bool partTwo = (CorrespondingSets[biter].find(aa) != CorrespondingSets[biter].end()) && (CorrespondingSets[aiter].find(bb) == CorrespondingSets[aiter].end());

                        if (partOne || partTwo)
                        {
                            boost::add_edge(aiter, biter,1.0, out);
                            std::cout << "Adding edge " << "[" << aiter << ", " << biter << "] : " << lulz[aiter] << " --- " << lulz[biter] << std::endl;
                        }


                    }

                }
            }

        }


    }

    Graph::edge_descriptor  my_find_edge(Graph::vertex_descriptor v, Graph::vertex_descriptor u, Graph const& g)
    {
        for (auto e : boost::make_iterator_range(boost::out_edges(v, g)))
        {
            if (target(e, g) == u)
                return e;
        }
        throw std::domain_error("my_find_edge: not found");
    }


    void constructFurerGraph(Graph & in, Graph & out, Graph & outTwo, std::vector<std::vector<int>> & MetaSets, std::vector<std::unordered_set<int>> & CorrespondingSets,
                             std::vector<int> & vertexIndxTracker, std::vector<std::vector<int>> & changess)
    {
        std::cout << "Launching thingh" << std::endl;
        std::cout << "Numvertices : " << boost::num_vertices(in) << std::endl;
        for (int i = 0; i < boost::num_vertices(in); i++)
        {
            std::cout << "Starting : " <<  i << " : " << std::endl;
            MetaSets.push_back(std::vector<int>());

            //     AdjacencyIterator j, endd;
//            for (boost::tie(j, endd) = boost::adjacent_vertices(i, in); j != endd; j++)
//            {
//                std::cout << *j << ", ";
//            }

            std::vector<int> setPotential;
            for (auto j = boost::adjacent_vertices(i,in).first; j != boost::adjacent_vertices(i,in).second; j++)
            {
                //  std::cout << *j << ", ";
                setPotential.push_back(*j);
            }
            std::unordered_set<int> setA;
            std::cout << "Go for recursion " << std::endl;
            recursiveVertexAdition(0, setPotential, in, i, setA, out, MetaSets, CorrespondingSets, vertexIndxTracker);
            std::cout << "Recursion finnished  " << std::endl;
            //!  int k = boost::add_vertex(out);
            //! std::cout << " added vertex: " << k << std::endl;
            //   std::cout << std::endl;



        }
        std::vector<std::string> lulz;

        for (int i = 0; i < boost::num_vertices(out); i++)
        {
            std::string wal = "";
            std::cout << i << " : " << vertexIndxTracker[i]+1 << " : [" ;
            wal = wal + std::to_string(vertexIndxTracker[i]+1) + " : [";
            for (int jj : CorrespondingSets[i])
            {
                std::cout  << jj+1 << ", ";
                wal = wal + std::to_string(jj+1) + ", ";

            }
            std::cout << "]" << std::endl;
            wal = wal + "]";
            lulz.push_back(wal);
        }

        std::cout << "Metasets: " << std::endl;
        for (int i = 0; i < boost::num_vertices(in); i++)
        {

            std::cout << i+1 << " : [ ";

            for (int jj = 0; jj < MetaSets[i].size(); jj++)
            {

                std::cout << MetaSets[i][jj] << ", ";

            }


            std::cout << "]" << std::endl;


        }
        std::cout << "Edges added " << std::endl;

        // Graph::edge_descriptor e = *(boost::edges(in).first);

        //!  Graph::edge_descriptor e = this->my_find_edge(1, 3, in);

        for (auto eit = boost::edges(in).first; eit != boost::edges(in).second; eit++)
        {

            int aa = boost::source(*eit, in);
            int bb = boost::target(*eit, in);

            // if (((aa == 1) && (bb == 3)) || ((aa == 3) && (bb == 1)))
            //  {
            //      e = *eit;

            //   }

            for (int i = 0; i < MetaSets[aa].size(); i++)
            {
                for (int j = 0; j < MetaSets[bb].size(); j++)
                {
                    int aiter = MetaSets[aa][i];
                    int biter = MetaSets[bb][j];
                    bool partOne = ((CorrespondingSets[biter].find(aa) != CorrespondingSets[biter].end()) && (CorrespondingSets[aiter].find(bb) != CorrespondingSets[aiter].end()));
                    bool partTwo = (CorrespondingSets[biter].find(aa) == CorrespondingSets[biter].end()) && (CorrespondingSets[aiter].find(bb) == CorrespondingSets[aiter].end());
                    if (partOne || partTwo)
                    {
                        boost::add_edge(aiter, biter, 1.0, out);
                        std::cout << "Adding edge " << "[" << aiter << ", " << biter << "] : " << lulz[aiter] << " --- " << lulz[biter] << std::endl;


                    }

                }
            }

        }



        invertGraph(in,outTwo, boost::num_vertices(out),  MetaSets, CorrespondingSets, vertexIndxTracker, changess);




        std::cout << "Printing all edges" << std::endl;
//         for (auto eit = boost::edges(in).first; eit != boost::edges(in).second; eit++)
//        {
//            Graph::edge_descriptor e = *eit;
////!            std::cout << "Source: " << boost::source(e,in) << " - "  << boost::target(e,in) << std::endl;
//
//            invertGraph(in, e);
//
//
//        }

        // boost::adjacent_vertices()
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

    void SingularRun(std::string xyzFiles, int a, int b)
    {

        namespace fs = std::filesystem;
        using namespace std;

        std::set<fs::path> sorted_by_name;

        for (auto &entry : fs::directory_iterator(xyzFiles))
            sorted_by_name.insert(entry.path());

        std::vector<std::vector<std::vector<double>>> matrices;
         std::vector<std::string> fileNamez;


        for (auto &filename : sorted_by_name)
        {
            std::string opp = filename.c_str();
            fileNamez.push_back(opp);
         //!   std::cout << filename.c_str() << std::endl;
            // std::vector<std::string> strings;
            // customSplit(opp, '/', strings);
            // std::string finall = strings[strings.size()-1];
            // out << "Finall (standard): " << finall << std::endl;
            // out << "Read path (standard): " << opp<< std::endl;
            //  out << "========" << std::endl;
            //   std::cout << opp << std::endl;
            //arma::mat k;
            std::vector<std::vector<double>> k;
            // ReadFileIntoMatrix(opp, k);
            readCSV(opp, k );
            matrices.push_back(k);
            //!epicMap[finall] = k;

        }
        std::cout << "A: " << fileNamez[a] << std::endl;
        std::cout << "B: " << fileNamez[b] << std::endl;

        long totalTimeForSMIstdXX;
        long totalTimeForBDXX;
        std::vector<int> matchingXX;
        std::vector<std::vector<int>> res;
        double bd = 0;

      //  auto firstP = std::make_pair(matrices[u], graphs[u]);
      //  auto secondP = std::make_pair( matrices[v], graphs[v]);
        auto firstP = matrices[a];
        auto secondP = matrices[b];

             typedef BottleneckCalculator<Bipartie_Graph,GraphTreeController> bipartiebd;
        typedef bd_general_interface<bipartiebd, bipartiebd> naiveBdCalc;

        EuclideanBDCompute<bipartiebd> eucBd;
        Nested_Updated<NCD_Graph<naiveBdCalc>, std::vector<std::vector<double>>> nestedNaive;
        Nested_Updated<NCD_Graph<naiveBdCalc>, std::vector<std::vector<double>>> nestedNaiveTwo;


        double bud = nestedNaiveTwo.Compute_Isometry(firstP, secondP, matchingXX, totalTimeForSMIstdXX, totalTimeForBDXX, 1);

        std::cout << "Final answer is : " << bud << std::endl;




    }

    void MassiveExperimentsRun(std::string xyzFiles, std::string csvInput, std::string csvOutput, int threads, int maxxnumm, int minnumm = 0)
    {

        namespace fs = std::filesystem;
        using namespace std;

        std::set<fs::path> sorted_by_name;

        for (auto &entry : fs::directory_iterator(xyzFiles))
            sorted_by_name.insert(entry.path());

        std::vector<std::vector<std::vector<double>>> matrices;




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
            //arma::mat k;
            std::vector<std::vector<double>> k;
            // ReadFileIntoMatrix(opp, k);
            readCSV(opp, k );
            matrices.push_back(k);
            //!epicMap[finall] = k;

        }

        std::cout << "Finnished up with matrices now conclude experiments" << std::endl;

        std::cout << "One  two three" << std::endl;
      //!  RunThisNCD(csvInput, csvOutput, threads, maxxnumm, matrices, minnumm );
            std::string header ="Molecule_A,Molecule_B,Num_Atoms,L_inf on SRV, L_inf on SDV,EMD on PDD, NBM_NCD1,time_NCD1(microseconds),NBD_NCD1_time(microseconds)";
       RunThisNCDSuperV1(csvInput, csvOutput, threads, maxxnumm, matrices, minnumm, header, 1);
    }

      void MassiveExperimentsRunTwo(std::string xyzFiles, std::string csvInput, std::string csvOutput, int threads, int maxxnumm, int minnumm = 0)
    {

        namespace fs = std::filesystem;
        using namespace std;

        std::set<fs::path> sorted_by_name;

        for (auto &entry : fs::directory_iterator(xyzFiles))
            sorted_by_name.insert(entry.path());

        std::vector<std::vector<std::vector<double>>> matrices;




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
            //arma::mat k;
            std::vector<std::vector<double>> k;
            // ReadFileIntoMatrix(opp, k);
            readCSV(opp, k );
            matrices.push_back(k);
            //!epicMap[finall] = k;

        }

        std::cout << "Finnished up with matrices now conclude experiments" << std::endl;

        std::cout << "One  two three" << std::endl;
      //!  RunThisNCD(csvInput, csvOutput, threads, maxxnumm, matrices, minnumm );
      std::string header ="indx_A,indx_B,Number_Atoms,Time_SRV,SRV,SDV,SDV_time,emd/pdd,time_emd,time_NCD1(microseconds),NBD_NCD1_time(microseconds),NBM_NCD1,time_NCD2(microseconds),NBD_NCD2_time(microseconds),NBM_NCD2";
       RunThisNCDSuperV1(csvInput, csvOutput, threads, maxxnumm, matrices, minnumm,header, 2);
    }

     void SingleExperimentRunVer2(int i, std::vector<std::vector<std::string>> & parsedCsv, std::vector<std::vector<std::string>> & outVector,std::vector<std::vector<std::vector<double>>> & matrices, int k)
    {
        //LayeredSCD lap(0.43,2,"");
        int tt = matrices.size();
     // y = std::min(y, tt);
        //  std::string debugStream = "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/GraphTypeInvariants/EMDPDD/debug.txt";
         typedef BottleneckCalculator<Bipartie_Graph,GraphTreeController> bipartiebd;
        typedef bd_general_interface<bipartiebd, bipartiebd> naiveBdCalc;

        EuclideanBDCompute<bipartiebd> eucBd;
        Nested_Updated<NCD_Graph<naiveBdCalc>, std::vector<std::vector<double>>> nestedNaive;
        Nested_Updated<NCD_Graph<naiveBdCalc>, std::vector<std::vector<double>>> nestedNaiveTwo;
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
        auto firstP = matrices[u];
        auto secondP = matrices[v];


        //double bud = nestedNaiveTwo.Compute_Isometry(firstP, secondP, matchingXX, totalTimeForSMIstdXX, totalTimeForBDXX, 1);


        bd = nestedNaive.Compute_Isometry(firstP, secondP, matching, totalTimeForSMIstd, totalTimeForBD, k);
        double bdtruncate = -1;

        long totalTimeForSMIstdfake;
        long totalTimeForBDfake;


        std::vector<int> matchingOne;
      //  double purebd = eucBd.compute(matrices[u], matrices[v], matchingOne);
      //  double haussdorff = computeHaudsorff(matrices[u], matrices[v]);
        double realDist = 0;
        double ratiodistbd = 0;
        bdtruncate = ratiodistbd;


            //! to_string_with_precision<double>(haussdorff , 15)
        tmp.push_back(parsedCsv[i][0]);
        tmp.push_back(parsedCsv[i][1]);
        tmp.push_back(parsedCsv[i][2]);
      //  tmp.push_back(to_string_with_precision<double>(haussdorff, 15));
      //  tmp.push_back(to_string_with_precision<double>(purebd, 15));
        tmp.push_back(parsedCsv[i][3]);
        tmp.push_back(parsedCsv[i][4]);
        tmp.push_back(parsedCsv[i][5]);
        tmp.push_back(to_string_with_precision<double>(bd, 15));
        tmp.push_back(parsedCsv[i][6]);
        tmp.push_back(parsedCsv[i][7]);
        tmp.push_back(parsedCsv[i][8]);
     //   tmp.push_back(to_string_with_precision<double>(bud, 15));
    //    double ratiod = bd / purebd;
      //  tmp.push_back(to_string_with_precision<double>(ratiod, 15));
     //   tmp.push_back(std::to_string(totalTimeForSMIstd));
     //   tmp.push_back(std::to_string(totalTimeForBD));
        tmp.push_back(std::to_string(totalTimeForSMIstdXX));
        tmp.push_back(std::to_string(totalTimeForBDXX));
        //!tmp.push_back(std::to_string(bdtruncate));
            //! totalTimeForSMIstdXX, totalTimeForBDXX

        outVector.push_back(tmp);

            //  std::cout << "Step six" << std::endl;




        // std::cout << "Outvector size: " << outVector.size() << std::endl;


    }

    void SingleExperimentRunSuperVer1(int i, std::vector<std::vector<std::string>> & parsedCsv, std::vector<std::vector<std::string>> & outVector,std::vector<std::vector<std::vector<double>>> & matrices, int kok)
    {
        //LayeredSCD lap(0.43,2,"");
        int tt = matrices.size();
     // y = std::min(y, tt);
        //  std::string debugStream = "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/GraphTypeInvariants/EMDPDD/debug.txt";
         typedef BottleneckCalculator<Bipartie_Graph,GraphTreeController> bipartiebd;
        typedef bd_general_interface<bipartiebd, bipartiebd> naiveBdCalc;

        EuclideanBDCompute<bipartiebd> eucBd;
        Nested_Updated<NCD_Graph<naiveBdCalc>, std::vector<std::vector<double>>> nestedNaive;
        Nested_Updated<NCD_Graph<naiveBdCalc>, std::vector<std::vector<double>>> nestedNaiveTwo;
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
        auto firstP = matrices[u];
        auto secondP = matrices[v];


       // double bud = nestedNaiveTwo.Compute_Isometry(firstP, secondP, matchingXX, totalTimeForSMIstdXX, totalTimeForBDXX, 1);

        //std::cout << "Bud is computed : " << bud << std::endl;


        double bud= nestedNaive.Compute_Isometry(firstP, secondP, matching, totalTimeForSMIstd, totalTimeForBD, kok);
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


      void SingleExperimentRun(int i, std::vector<std::vector<std::string>> & parsedCsv, std::vector<std::vector<std::string>> & outVector,std::vector<std::vector<std::vector<double>>> & matrices )
    {
        //LayeredSCD lap(0.43,2,"");
        int tt = matrices.size();
     // y = std::min(y, tt);
        //  std::string debugStream = "/home/yury/LocalProjects/XYZProject/DistanceBasedInavriance/GraphTypeInvariants/EMDPDD/debug.txt";
         typedef BottleneckCalculator<Bipartie_Graph,GraphTreeController> bipartiebd;
        typedef bd_general_interface<bipartiebd, bipartiebd> naiveBdCalc;

        EuclideanBDCompute<bipartiebd> eucBd;
        Nested_Updated<NCD_Graph<naiveBdCalc>, std::vector<std::vector<double>>> nestedNaive;
        Nested_Updated<NCD_Graph<naiveBdCalc>, std::vector<std::vector<double>>> nestedNaiveTwo;
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
        auto firstP = matrices[u];
        auto secondP = matrices[v];


        double bud = nestedNaiveTwo.Compute_Isometry(firstP, secondP, matchingXX, totalTimeForSMIstd, totalTimeForBD, 1);


        //!bd = nestedNaive.Compute_Isometry(firstP, secondP, matching, totalTimeForSMIstd, totalTimeForBD, 2);
        double bdtruncate = -1;

        long totalTimeForSMIstdfake;
        long totalTimeForBDfake;


        std::vector<int> matchingOne;
      //  double purebd = eucBd.compute(matrices[u], matrices[v], matchingOne);
      //  double haussdorff = computeHaudsorff(matrices[u], matrices[v]);
        double realDist = 0;
        double ratiodistbd = 0;
        bdtruncate = ratiodistbd;


            //! to_string_with_precision<double>(haussdorff , 15)
        tmp.push_back(parsedCsv[i][0]);
        tmp.push_back(parsedCsv[i][1]);
        tmp.push_back(parsedCsv[i][2]);
      //  tmp.push_back(to_string_with_precision<double>(haussdorff, 15));
      //  tmp.push_back(to_string_with_precision<double>(purebd, 15));
        tmp.push_back(parsedCsv[i][3]);
        tmp.push_back(parsedCsv[i][4]);
        tmp.push_back(parsedCsv[i][5]);
      //!  tmp.push_back(to_string_with_precision<double>(bd, 15));
        tmp.push_back(to_string_with_precision<double>(bud, 15));
    //    double ratiod = bd / purebd;
      //  tmp.push_back(to_string_with_precision<double>(ratiod, 15));
     //   tmp.push_back(std::to_string(totalTimeForSMIstd));
     //   tmp.push_back(std::to_string(totalTimeForBD));
        tmp.push_back(std::to_string(totalTimeForSMIstdXX));
        tmp.push_back(std::to_string(totalTimeForBDXX));
        //!tmp.push_back(std::to_string(bdtruncate));
            //! totalTimeForSMIstdXX, totalTimeForBDXX

        outVector.push_back(tmp);

            //  std::cout << "Step six" << std::endl;




        // std::cout << "Outvector size: " << outVector.size() << std::endl;


    }

    template <typename T>
    std::string to_string_with_precision(const T a_value, const int n = 6)
    {
        std::ostringstream out;
        out.precision(n);
        out << std::fixed << a_value;
        return std::move(out).str();
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

    void RunThisNCDSuperV1(std::string inputCsv, std::string outputCsv, int threads, int maxxnumm, std::vector<std::vector<std::vector<double>>> & matrices , int minnumm = 0, std::string header = "", int kk = 1)
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
           pool.push_task(&GraphExperiments::SingleExperimentRunSuperVer1, this, i, std::ref(parsedCSV), std::ref(vectorCollection[i]), std::ref(matrices) , kk);
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



    void RunThisNCD(std::string inputCsv, std::string outputCsv, int threads, int maxxnumm, std::vector<std::vector<std::vector<double>>> & matrices , int minnumm = 0)
    {
      //  maxxnumm = std::min(maxxnumm, (int) matrices.size());
        BS::thread_pool pool(threads);
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
           pool.push_task(&GraphExperiments::SingleExperimentRunVer2, this, i, std::ref(parsedCSV), std::ref(vectorCollection[i]), std::ref(matrices), 1);
       //!  GraphExperiments::SingleExperimentRun( i, parsedCSV, vectorCollection[i], matrices);
            //  pool.push_task([i, &parsedCSV, &vectorCollection[i]]
            //   {
            //      this->limitedExperRun(i, i+1,parsedCSV)
            //   }
            //  );



        }

       pool.wait_for_tasks();

        //!std::string header ="Molecule_A,Molecule_B,Num_Atoms,Hausdorff_distance,BD,L_inf on SRV, L_inf on SDV,EMD on PDD,NBM_NCD2, NBM_NCD1 ,NBM_NDP/BD,time_NDP(microseconds),NBD_time(microseconds),Match_dist/NBM_NDP";
    std::string header ="Molecule_A,Molecule_B,Num_Atoms,L_inf on SRV,L_inf on SDV,EMD on PDD,NBM_NCD2, NBM_NCD1 ,time_NCD2(microseconds),NBD_NCD2_time(microseconds),time_NCD1(microseconds),NBD_NCD1_time(microseconds)";

     //!  std::string header ="Molecule_A,Molecule_B,Num_Atoms,L_inf on SRV, L_inf on SDV,EMD on PDD, NBM_NCD1,time_NCD1(microseconds),NBD_NCD1_time(microseconds)";

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
            double aa = stod(a[6]);
            double bb = stod(b[6]);
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


    void RunThisNCDVer2(std::string inputCsv, std::string outputCsv, int threads, int maxxnumm, std::vector<std::vector<std::vector<double>>> & matrices )
    {
      //  maxxnumm = std::min(maxxnumm, (int) matrices.size());
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
            pool.push_task(&GraphExperiments::SingleExperimentRun, this, i, std::ref(parsedCSV), std::ref(vectorCollection[i]), std::ref(matrices));
            //  pool.push_task([i, &parsedCSV, &vectorCollection[i]]
            //   {
            //      this->limitedExperRun(i, i+1,parsedCSV)
            //   }
            //  );



        }

        pool.wait_for_tasks();

        //!std::string header ="Molecule_A,Molecule_B,Num_Atoms,Hausdorff_distance,BD,L_inf on SRV, L_inf on SDV,EMD on PDD,NBM_NCD2, NBM_NCD1 ,NBM_NDP/BD,time_NDP(microseconds),NBD_time(microseconds),Match_dist/NBM_NDP";
        std::string header ="Molecule_A,Molecule_B,Num_Atoms,L_inf on SRV, L_inf on SDV,EMD on PDD,NBM_NCD2, NBM_NCD1 ,time_NCD2(microseconds),NBD_NCD2_time(microseconds),time_NCD1(microseconds),NBD_NCD1_time(microseconds)";

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
            double aa = stod(a[6]);
            double bb = stod(b[6]);
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




    void readCSV(std::string fname, std::vector<std::vector<int>> & matt)
    {

        namespace fs = std::filesystem;
        using namespace std;
        //  std::cout << "Beggining to read" << std::endl;

        std::vector<std::vector<std::string>> content;
        std::vector<string> row;
        std::vector<int> rowd;
        string line, word;


        std::fstream file (fname, ios::in);


        if(file.is_open())
        {
            std::cout << "Reading ... " << std::endl;
            while(getline(file, line))
            {

                row.clear();
                rowd.clear();

                stringstream str(line);

                while(getline(str, word, ','))
                {
                    row.push_back(word);
                    rowd.push_back(stoi(word));

                }
                content.push_back(row);
                matt.push_back(rowd);
            }
        }
        std::cout << "content: " << content.size() << " , " << content[0].size() << std::endl;


    }
    template<typename dataType>
    void printVector(std::vector<dataType> & in)
    {
        //  std::cout << "Size of vector: " << in.size() << std::endl;
        for (int i = 0; i < in.size(); i++)
        {

            std::cout << in[i] << ", ";
        }
        std::cout << std::endl;
    }

    template<typename dataType>
    void printMatrix(std::vector<std::vector<dataType>> & in)
    {
        //  std::cout << "Size of vector: " << in.size() << std::endl;
        for (int i = 0; i < in.size(); i++)
        {
            printVector(in[i]);
        }
    }

    void testAdjancency(std::string inputA, std::string ffm, std::vector<std::vector<int>> & changess, int numm = 2)
    {
        std::vector<std::vector<int>> mattA;
        std::vector<std::vector<int>> mattB;

        std::cout << "Preparing to read CSV: " << inputA << " Debug file: " << ffm << std::endl;

        readCSV(inputA, mattA);

        std::cout << "Input succesfully read" << std::endl;

        printMatrix(mattA);

        Graph ff;
        convertMatToGraph(mattA, ff);



        Graph gA;
        Graph gB;
        std::vector<std::vector<int>> dumpOne;
        std::vector<std::unordered_set<int>> dumbTwo;
        std::vector<int> vertexIndxTracker;
        constructFurerGraph(ff, gA, gB, dumpOne, dumbTwo, vertexIndxTracker, changess);

        //  convertMatToGraph(mattA, gA);
        //  convertMatToGraph(mattB, gB);
        std::vector<std::vector<double>> distmattA;
        std::vector<std::vector<double>> distmattB;
        convertGraphToDistMat(gA, distmattA);
        convertGraphToDistMat(gB, distmattB);

        std::cout << "Numvertices A: " << boost::num_vertices(gA) << " num edges: " << boost::num_edges(gA) << std::endl;
        std::cout << "Numvertices B: " << boost::num_vertices(gB) << " num edges: " << boost::num_edges(gB) << std::endl;

        std::cout << "DistMatrix A: | " << distmattA.size() << " , " << distmattA[0].size() << std::endl;
        printMatrix(distmattA);
        std::cout << "DistMatrix B: | " << distmattB.size() << " , " << distmattB[0].size() << std::endl;
        printMatrix(distmattB);
        //    long  totalTimeForSMIstd;
        //    long  totalTimeForBD;
        double bd;
        std::vector<int> matching;
        long totalTimeForSMIstd = 0;
        long totalTimeForBD = 0;
        std::cout << "Starting computations" << std::endl;


        //   std::string ffm = "/home/yury/LocalProjects/ReaGraphStuff/TestGraphs/debug/outtest.txt";
        std::ofstream k;
        k.open(ffm);
        this->computeBetweenDistMats(distmattA, distmattB, bd, matching,totalTimeForSMIstd, totalTimeForBD,ffm, numm);
        k.close();
        std::cout << "Bottlneckdistance computed as : " << bd << std::endl;
        std::cout << "Invariant construction time: " << totalTimeForSMIstd << std::endl;
        std::cout << "Invariant construction time: " <<totalTimeForBD  << std::endl;
        for (int i = 0; i < matching.size() ; i++)
        {
            std::cout << matching[i] << ", ";
        }
        std::cout << std::endl;
    }

    void vectorToArmaMat(const std::vector<std::vector<double>>& vec, arma::mat & result)
    {


        size_t rows = vec.size();
        size_t cols = vec[0].size();

        result = arma::mat(rows, cols);

        for (size_t i = 0; i < rows; ++i)
        {
            if (vec[i].size() != cols)
            {
                throw std::runtime_error("Inconsistent number of columns in input vector.");
            }
            for (size_t j = 0; j < cols; ++j)
            {
                result(i, j) = vec[i][j];
            }
        }

//   return result;
    }


    void launchSimpleTest(std::string inputA, std::string inputB)
    {
        std::vector<std::vector<double>> mattA;
        std::vector<std::vector<double>> mattB;
        readCSV(inputA, mattA);
        readCSV(inputB, mattB);
        std::cout << "Matrix A: | " << mattA.size() << " , " << mattA[0].size() << std::endl;
        printMatrix(mattA);
        std::cout << "Matrix B: | " << mattB.size() << " , " << mattB[0].size() << std::endl;
        printMatrix(mattB);

        int numm = 2;
        std::string ffm = "";
        double bd = 0;
        std::vector<int> matching;
        long totalTimeForSMIstd;
        long totalTimeForBD;
        arma::mat distmattA;
        arma::mat distmattB;
        vectorToArmaMat(mattA, distmattA);
        vectorToArmaMat(mattB, distmattB);

        typedef BottleneckCalculator<Bipartie_Graph,GraphTreeController> bipartiebd;
        typedef bd_general_interface<bipartiebd, bipartiebd> naiveBdCalc;

        EuclideanBDCompute<bipartiebd> eucBd;
        Nested_Updated<NCD_Graph<naiveBdCalc>, std::vector<std::vector<double>>> nestedNaive;


        double dd = nestedNaive.Compute_Isometry(mattA, mattB, matching, totalTimeForSMIstd, totalTimeForBD, numm);


        //   k.close();
        bd = dd;

        //   this->computeBetweenDistMats(distmattA, distmattB, bd, matching,totalTimeForSMIstd, totalTimeForBD,ffm, numm);
        //  k.close();
        std::cout << "Bottlneckdistance computed as : " << bd << std::endl;
        std::cout << "Invariant construction time: " << totalTimeForSMIstd << std::endl;
        std::cout << "Invariant construction time: " <<totalTimeForBD  << std::endl;
        for (int i = 0; i < matching.size() ; i++)
        {
            std::cout << matching[i] << ", ";
        }
        std::cout << std::endl;

    }

    void launchTest(std::string inputA, std::string inputB)
    {
        std::vector<std::vector<int>> mattA;
        std::vector<std::vector<int>> mattB;
        readCSV(inputA, mattA);
        readCSV(inputB, mattB);
        std::cout << "Matrix A: | " << mattA.size() << " , " << mattA[0].size() << std::endl;
        printMatrix(mattA);
        std::cout << "Matrix B: | " << mattB.size() << " , " << mattB[0].size() << std::endl;
        printMatrix(mattB);


        Graph gA;
        Graph gB;
        convertMatToGraph(mattA, gA);
        convertMatToGraph(mattB, gB);
        std::vector<std::vector<double>> distmattA;
        std::vector<std::vector<double>> distmattB;
        convertGraphToDistMat(gA, distmattA);
        convertGraphToDistMat(gB, distmattB);

        std::cout << "DistMatrix A: | " << distmattA.size() << " , " << distmattA[0].size() << std::endl;
        printMatrix(distmattA);
        std::cout << "DistMatrix B: | " << distmattB.size() << " , " << distmattB[0].size() << std::endl;
        printMatrix(distmattB);



        //    long  totalTimeForSMIstd;
        //    long  totalTimeForBD;
        double bd;
        std::vector<int> matching;
        long totalTimeForSMIstd = 0;
        long totalTimeForBD = 0;
        std::cout << "Starting computations" << std::endl;
        this->computeBetweenDistMats(distmattA, distmattB, bd, matching,totalTimeForSMIstd, totalTimeForBD);
        std::cout << "Bottlneckdistance computed as : " << bd << std::endl;
        std::cout << "Invariant construction time: " << totalTimeForSMIstd << std::endl;
        std::cout << "Invariant construction time: " <<totalTimeForBD  << std::endl;
        for (int i = 0; i < matching.size() ; i++)
        {
            std::cout << matching[i] << ", ";
        }
        std::cout << std::endl;



    }

//   void readCSV(std::string input , std::vector<std::vector<int>> & graph)
    //  {




    //  }

    void convertMatToGraph(std::vector<std::vector<int>> & mat, Graph & G)
    {
        for (int i = 0; i < mat.size(); i++)
        {
            for (int j = 0; j < i; j++)
            {
                if (i != j)
                {
                    if (mat[i][j] == 1)
                    {
                        boost::add_edge(i,j, 1.0, G);
                    }
                }
            }

        }



    }

    void convertGraphToDistMat(Graph & g, std::vector<std::vector<double>> & matOut)
    {
        matOut = std::vector<std::vector<double>>(boost::num_vertices(g), std::vector<double>(boost::num_vertices(g)));
        WeightMap weight_pmap = boost::get(boost::edge_weight, g);
        DistanceMatrix distances(num_vertices(g));
        DistanceMatrixMap dm(distances, g);

        // find all pairs shortest paths
        bool valid = floyd_warshall_all_pairs_shortest_paths(g, dm, boost::weight_map(weight_pmap));
        for (std::size_t i = 0; i < num_vertices(g); ++i)
        {
            for (std::size_t j = 0; j < num_vertices(g); ++j)
            {
                matOut[i][j] = distances[i][j];
                //  std::cout << "From vertex " << i+1 << " to " << j+1 << " : ";
//                if(distances[i][j] == std::numeric_limits<t_weight>::max())
//                    std::cout << "inf" << std::endl;
//                else
//                    std::cout << distances[i][j] << std::endl;
            }
            // std::cout << std::endl;
        }
    }
    void computeBetweenDistMats(std::vector<std::vector<double>> & A, std::vector<std::vector<double>> & B, double & bd, std::vector<int> & matching,
                                long & totalTimeForSMIstd, long & totalTimeForBD, int numm = 2)
    {

        typedef BottleneckCalculator<Bipartie_Graph,GraphTreeController> bipartiebd;
        typedef bd_general_interface<bipartiebd, bipartiebd> naiveBdCalc;

        EuclideanBDCompute<bipartiebd> eucBd;
        Nested_Updated<NCD_Graph<naiveBdCalc>, std::vector<std::vector<double>>> nestedNaive;
        //!  Nested_Base<SCD, naiveBdCalc, std::pair<arma::mat, std::vector<std::vector<bool>>>> nestedNaive;
        //!    Nested_Base<PBI, naiveBdCalc, arma::mat> nestedNaive;




        // std::vector<std::vector<int>> res;
        // double bd = 0;

        //!     auto firstP = std::make_pair(matrices[u], graphs[u]);
        //!    auto secondP = std::make_pair( matrices[v], graphs[v]);

        // std::cout << "Clouds are : " << std::endl;
        //  std::cout << "Dimension of A: " << matrices[u].n_cols << ", " << matrices[u].n_rows << std::endl;
        //   std::cout << "Dimension of B: " << matrices[v].n_cols << ", " << matrices[v].n_rows << std::endl;



        //!    bd = nestedNaive.Compute_Isometry(firstP, secondP, matching, totalTimeForSMIstd, totalTimeForBD);

        //   std::string ffm = "/home/yury/LocalProjects/ReaGraphStuff/TestGraphs/debug/outtest.txt";
        //   std::ofstream k;
        //  k.open(ffm);

        double dd = nestedNaive.Compute_Isometry(A, B, matching, totalTimeForSMIstd, totalTimeForBD, numm);


        //   k.close();
        bd = dd;




    }

    void computeBetweenDistMats(std::vector<std::vector<double>> & A, std::vector<std::vector<double>> & B, double & bd, std::vector<int> & matching,
                                long & totalTimeForSMIstd, long & totalTimeForBD, std::string ffm, int numm = 2)
    {

        typedef BottleneckCalculator<Bipartie_Graph,GraphTreeController> bipartiebd;
        typedef bd_general_interface<bipartiebd, bipartiebd> naiveBdCalc;

        EuclideanBDCompute<bipartiebd> eucBd;
        Nested_Updated<NCD_Graph<naiveBdCalc>, std::vector<std::vector<double>>> nestedNaive;
        //!  Nested_Base<SCD, naiveBdCalc, std::pair<arma::mat, std::vector<std::vector<bool>>>> nestedNaive;
        //!    Nested_Base<PBI, naiveBdCalc, arma::mat> nestedNaive;




        // std::vector<std::vector<int>> res;
        // double bd = 0;

        //!     auto firstP = std::make_pair(matrices[u], graphs[u]);
        //!    auto secondP = std::make_pair( matrices[v], graphs[v]);

        // std::cout << "Clouds are : " << std::endl;
        //  std::cout << "Dimension of A: " << matrices[u].n_cols << ", " << matrices[u].n_rows << std::endl;
        //   std::cout << "Dimension of B: " << matrices[v].n_cols << ", " << matrices[v].n_rows << std::endl;



        //!    bd = nestedNaive.Compute_Isometry(firstP, secondP, matching, totalTimeForSMIstd, totalTimeForBD);

        //   std::string ffm = "/home/yury/LocalProjects/ReaGraphStuff/TestGraphs/debug/outtest.txt";
        std::ofstream k;
        k.open(ffm);

        std::cout << "Running with numm: " << numm << std::endl;

        double dd = nestedNaive.Compute_Isometry(A, B, matching, totalTimeForSMIstd, totalTimeForBD, numm,k);


        k.close();
        bd = dd;




    }


    virtual ~GraphExperiments() {}

protected:

private:
};

#endif // GRAPHEXPERIMENTS_H
