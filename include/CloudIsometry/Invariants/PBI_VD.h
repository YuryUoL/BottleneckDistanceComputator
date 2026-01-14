#ifndef PBI_VD_H
#define PBI_VD_H

template<typename bdcomputatorType>
class PBI_VD
{
public:
    PBI_VD()
    {
        sortedIndices.push_back(std::vector<size_t>());
        sortedIndices.push_back(std::vector<size_t>());
        vds.push_back(std::vector<std::vector<double>>());
        vds.push_back(std::vector<std::vector<double>>());
        rads.push_back(0);
        rads.push_back(0);
        ocd.push_back(std::vector<arma::mat>());
        dists.push_back(std::vector<arma::mat>());
        ocd.push_back(std::vector<arma::mat>());
        dists.push_back(std::vector<arma::mat>());
        radclouds.push_back(std::vector<double>());
        radclouds.push_back(std::vector<double>());
        //   constC = sqrt(3);
        aConst = 14;
        furthestPoint.push_back(0);
        furthestPoint.push_back(0);
        minVall = std::numeric_limits<double>::max();
    }
    template<typename typ>
    void printVector(std::vector<typ> & in, std::ofstream & debugS)
    {
        //  std::cout << "Size of vector: " << in.size() << std::endl;
        for (int i = 0; i < in.size(); i++)
        {

            debugS << in[i] << ", ";
        }
        debugS << std::endl;
    }

    template<typename typ>
    void printVectorDouble(std::vector<typ> & in, std::vector<typ> & intwo, std::ofstream & debugS)
    {
        //  std::cout << "Size of vector: " << in.size() << std::endl;
        for (int i = 0; i < in.size(); i++)
        {

            debugS << in[i] << ", ";
        }
        for (int i = 0; i < intwo.size(); i++)
        {

            debugS << intwo[i] << ", ";
        }
        debugS << std::endl;
    }


    int NumElem(arma::mat & cloudAunc)
    {
        return cloudAunc.n_cols;

    }
    int Dimension(arma::mat & cloudAunc)
    {

        return cloudAunc.n_rows;

    }

    void initilize(arma::mat & cloudAunc, arma::mat & cloudA, int index, std::ofstream & debugS)
    {


        double radius = initilize(cloudAunc,  cloudA, index);
        debugS << "Radius " << index << " set to: " << radius << std::endl;

        //     std::cout << cloudAunc << std::endl;
        //   std::cout << "index : " << index << std::endl;
        //  std::cout << "Rads size: " << rads.size() << std::endl;
        //   std::cout << "Radus is: " << radius << std::endl;
        //  std::cout << "Rads[index]: " << rads[index] << std::endl;
        // std::cout << "Radius computed as " << radius << std:endl;
        //    return radius;
    }

    double initilize(arma::mat & cloudAunc, arma::mat & cloudA, int index)
    {
        // std::cout << "Initilizing cloud " << index << std::endl;

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

        //    std::cout << "Size is " << cloudOUT.n_cols << std::endl;

        std::vector<std::pair<double, size_t>> dists;

        for (int i = 0; i < cloudOUT.n_cols; i++)
        {

            //    std::cout << "Itearting over " << i << " / " << cloudOUT.n_cols << std::endl;
            arma::vec blablaz = cloudAunc.col(i);
            arma::vec tmpd = blablaz + centroid;
            cloudOUT.col(i) = tmpd;
            double dudu = arma::norm(tmpd,2);

            dists.emplace_back(dudu, i);

            //   std::cout << "Dudu is computed as " << dudu << std::endl;

            if (dudu > radius)
            {
                //     std::cout << "Updating radius " << std::endl;
                radius = dudu;
                //  furthestPoint[index] = i;

            }
            //   radius = std::max(radius,dudu);

            //! radius = std::max(radius, sqrt(vectorTools::SquaredNorm(cloudOUT[i])));

        }

        std::sort(dists.begin(), dists.end(), [](const auto& a, const auto& b)
        {
            return a.first > b.first;
        });

// Extract sorted indices
        //  std::vector<size_t> sortedIndices;
        sortedIndices.reserve(dists.size());
        for (const auto& p : dists)
        {
            sortedIndices[index].push_back(p.second);

        }







        cloudA = cloudOUT;
        rads[index] = radius;
        return radius;


        //     std::cout << cloudAunc << std::endl;
        //   std::cout << "index : " << index << std::endl;
        //  std::cout << "Rads size: " << rads.size() << std::endl;
        //   std::cout << "Radus is: " << radius << std::endl;
        //  std::cout << "Rads[index]: " << rads[index] << std::endl;
        // std::cout << "Radius computed as " << radius << std:endl;
        //    return radius;
    }

    void MultiDCrossProduct(arma::mat & cloud, arma::vec & out)
    {
        out = arma::vec(cloud.n_rows);
        if (cloud.n_rows == 3)
        {
            double R = std::max(arma::norm(cloud.col(0)), arma::norm(cloud.col(1)) );
            arma::vec c_P = cross( cloud.col(0), cloud.col(1));
            out = c_P / R;

            //  out = cross( cloud.col(indices[0]), cloud.col(indices[1]) );
        }
        else if (cloud.n_rows == 2)
        {
            double a = -cloud.col(0)(1);
            double b = cloud.col(0)(0);

            std::vector<double> c_P_d{a,b};
            arma::vec c_P(c_P_d);
            out = c_P;
            //return c_P;

        }
        else
        {

            // return 0;
        }

    }


    arma::vec FindLastGuy(arma::mat & cloud, std::vector<int> & indices)
    {

        if (cloud.n_rows == 3)
        {
            double R = std::max(arma::norm(cloud.col(indices[0])), arma::norm(cloud.col(indices[1])) );
            arma::vec c_P = cross( cloud.col(indices[0]), cloud.col(indices[1]) );
            return c_P / R;
        }
        else if (cloud.n_rows == 2)
        {
            double a = -cloud.col(indices[0])(1);
            double b = cloud.col(indices[0])(0);

            std::vector<double> c_P_d{a,b};
            arma::vec c_P(c_P_d);
            return c_P;

        }
        else
        {

            return 0;
        }

    }

    void computeBasisMatrix(std::vector<int> & indices,arma::mat & cloud, arma::mat & basisCollection, arma::mat & orthogonals)
    {

        basisCollection = arma::mat(cloud.n_rows, indices.size());
        // arma::vec c_P(cloud.n_rows);



        //! std::vector<double> c_P(cloud[indices[0]].size());
        // std::vector<std::vector<double>> basisCollection;
        for (int i = 0; i < indices.size(); i++)
        {
            basisCollection.col(i) = cloud.col(indices[i]);

            //!  basisCollection.push_back(cloud[indices[i]]);
        }
        //! vectorTools::multiDcrossProduct(basisCollection, c_P);
        //vectorTools::crossProduct(cloud[indices[0]], cloud[indices[1]], c_P);

        orthogonals = arma::mat(cloud.n_rows, indices.size());
        arma::vec c_P = FindLastGuy(cloud, indices);
        orthogonals.col(indices.size()-1) = c_P;
        for (int i = indices.size()-2; i >= 0;  i--)
        {
            int numbering = 0;
            arma::mat FindFromThisSet(cloud.n_rows, indices.size());
            arma::vec p;
            int realI = i + 1;
            for(int j = 0; j < indices.size() + 1; j++)
            {
                int jjj = j;
                if (jjj < realI)
                {
                    FindFromThisSet.col(numbering) = basisCollection.col(j);
                    numbering = numbering + 1;



                }
                else if (j > realI)
                {

                    FindFromThisSet.col(numbering) = orthogonals.col(j-1);
                    numbering = numbering + 1;

                }


            }

            MultiDCrossProduct(cloud, p);
            orthogonals.col(i) = p;


        }






    }

    void ConstructSRCloud( arma::mat & cloud, arma::mat & basis, arma::mat & orthogonals, std::vector<int> & indices,arma::mat & SrCloud, int index)
    {
        int m = cloud.n_cols;
        int n = cloud.n_rows;
        //!arma::vec c_P(cloud.n_rows);
        SrCloud = arma::mat(n+2, m - n + 1);
        //  std::cout << "SrCloud rows: " << SrCloud.n_rows << " , SrCloud columns: " << SrCloud.n_cols << std::endl;
        std::set<int> s(indices.begin(), indices.end());

        // std::cout << "First: " << cloud[i][0] << ", "<< cloud[i][1] << ", " << cloud[i][2] << std::endl;
        // std::cout << "Second: " << cloud[j][0] << ", "<< cloud[j][1] << ", " << cloud[j][2] << std::endl;
//        vectorTools::crossProduct(cloud[i], cloud[j], c_P);
        // std::cout << "Final: " << c_P[0] << ", "<< c_P[1] << ", " << c_P[2] << std::endl;
        int summ = 0;
//        double scalar = 0;
//        for (int i = 0; i < basis.n_cols - 1; i++)
//        {
//            scalar = std::max(scalar, arma::norm(basis.col(i), 2));
//
//        }
//        scalar = 1/scalar;

        for (int k = 0; k < cloud.n_cols ; k++)
        {
            // std::cout << "Iteration, k : " << k << " summ : " << summ << std::endl;
            if (s.find(k) == s.end())
            {
                std::vector<double> newVector;
                newVector.push_back(arma::norm(cloud.col(k), 2));

                for (int i = 0; i < basis.n_cols; i++)
                {
                    // double dottR = arma::dot(cloud.col(k), basis.col(i)) / rads[index];
                    //   dottR = dottR / pow(rads[index],(cloud.n_rows-2));
                    // newVector.push_back(dottR);

                    newVector.push_back(arma::norm(cloud.col(k) - basis.col(i), 2));


                }

                for (int i = 0; i < orthogonals.n_cols; i++)
                {
                    // double dottR = arma::dot(cloud.col(k), basis.col(i)) / rads[index];
                    //   dottR = dottR / pow(rads[index],(cloud.n_rows-2));
                    // newVector.push_back(dottR);

                    newVector.push_back(arma::norm(cloud.col(k) - orthogonals.col(i), 2));


                }



                //std::cout << "New vector size " << newVector.size() << std::endl;
                arma::vec p(newVector);
                // std::cout << "Size of the arma::vec" << p.size() << std::endl;
                SrCloud.col(summ) = p;

                //  std::cout << "Survived without error!" << std::endl;

                //SrCloud.push_back(newVector);
                summ = summ + 1;
            }

        }

        //  std::cout << "Reached end" << std::endl;

        //CrystData::PrintMatrix(SrCloud);

        // std::cout << "SrCloudMat dimensions: " << SrCloud.size() << " , " << SrCloud[i].size() << std::endl;
    }

    double ConstructDistanceMatrix(arma::mat & basisCollection, arma::mat & matrix, int index)
    {

        //!! Assume that dimension = 3
        matrix = arma::mat(basisCollection.n_cols + 1, basisCollection.n_cols + 1);
        arma::mat newBase(basisCollection.n_rows, basisCollection.n_cols + 1);

       // arma::mat newBase = basisCollection;
        arma::vec origin(basisCollection.n_rows,  arma::fill::zeros);
//        for (int i = 0; i < origin.size(); i++)
//        {
//            origin(i) = 0;
//        }

        for (int i = 0; i < basisCollection.n_cols; i++)
        {
            newBase.col(i+1) = basisCollection.col(i);


        }

        newBase.col(0) = origin;
//
        double maxRad = 0;
//
       for (int i = 0; i < basisCollection.n_cols; i++)
        {
            double pp = arma::norm(basisCollection.col(i));
            maxRad = std::max(maxRad, pp);

        }



        for (int i = 0; i < basisCollection.n_cols; i++)
        {
            for (int j = 0; j < basisCollection.n_cols; j++)
            {
                double ttt = arma::norm(newBase.col(i) - newBase.col(j), 2);
                matrix(i,j) = ttt;

                //   matrix(i,j) = ttt * ttt / maxRad;

                //!rads[index];

                //! ttt / rads[index]

//                std::vector<double> vecU = basisCollection[j];
//                vectorTools::vectorScalarMultiplication(vecU,-1);
//                std::vector<double> vecV = basisCollection[i];
//                vectorTools::vectorAddition(vecV, vecU);
//              //  double norm;
//                matrix[i][j] = vectorTools::SquaredNorm(vecV);
            }
        }
        return maxRad;
    }


    void ConstructVectorOfDistance(arma::mat & basisCollection, arma::mat & orthogonals, std::vector<double> & out)
    {

        for (int i = 1; i < basisCollection.n_cols; i++)
        {
            out.push_back(arma::norm(basisCollection.col(i), 2));
            out.push_back(arma::norm(basisCollection.col(i)-orthogonals.col(i-1), 2));
        }
        out.push_back(arma::norm(orthogonals.col(basisCollection.n_cols-1)));




    }

    void ConstructCloudInvariantCombo(arma::mat & cloud, std::vector<std::vector<int>> & combos, int index, std::ofstream & debugS)
    {
        //std::cout << "debug version of cloudinvariant combo" << std::endl;
        //   std::cout << "index : " << index << std::endl;
        //  std::cout << "Rads[index]: " << rads[index] << std::endl;
        //   std::cout << "We are entering ConstructionCloudInvariantCombo with size: " << combos.size() << std::endl;
        for (int i = 0; i < combos.size(); i++)
        {
            arma::mat basisCollection;
            arma::mat orthogonals;
            std::vector<double> vde;
            //! We are here:
            //!    CrystData::computeBasisMatrix(combos[i],cloud,basisCollection);
            computeBasisMatrix(combos[i],cloud, basisCollection, orthogonals );
            arma::mat SrCloud;
            arma::mat DistMat;


            //  CrystData::ConstructSRCloud(cloud,basisCollection,combos[i],SrCloud);
            //CrystData::ConstructDistanceMatrix(basisCollection,DistMat);
            ConstructVectorOfDistance(basisCollection, orthogonals, vde);
            ConstructSRCloud( cloud,basisCollection,orthogonals, combos[i],SrCloud,index);


            double maxRad = ConstructDistanceMatrix(basisCollection,DistMat, index);
            radclouds[index].push_back(maxRad);


            //CrystData::ConstructSRCloud(i,j,cloud,SrCloud);
            //    CrystData::ConstructDistanceMatrix(indices,cloud,DistMat);

            ocd[index].push_back(SrCloud);
            dists[index].push_back(DistMat);
            vds[index].push_back(vde);
//            if (i == 1)
//            {
//                std::cout << "Special prrint" << std::endl;
//
//
//                std::cout << "Anotehrrr" << std::endl;
//
//                PrintMatrix(SrCloudCollection[0]);
//            }

            debugS << "Current index: " << index << std::endl;
            debugS << "Combination "  << i << std::endl;
            printVector(combos[i],debugS);
            debugS << SrCloud << std::endl;
            debugS << DistMat << std::endl;




        }


    }

    void ConstructCloudInvariantCombo(arma::mat & cloud, std::vector<std::vector<int>> & combos, int index)
    {
      //  std::cout << "regular version of cloudinvariant combo" << std::endl;
        //   std::cout << "index : " << index << std::endl;
        //  std::cout << "Rads[index]: " << rads[index] << std::endl;
        //   std::cout << "We are entering ConstructionCloudInvariantCombo with size: " << combos.size() << std::endl;
        for (int i = 0; i < combos.size(); i++)
        {
            arma::mat basisCollection;
            arma::mat orthogonals;
            std::vector<double> vde;
            //! We are here:
            //!    CrystData::computeBasisMatrix(combos[i],cloud,basisCollection);
            computeBasisMatrix(combos[i],cloud, basisCollection, orthogonals );

         //   std::cout << "Basis matrix okay " << std::endl;

            arma::mat SrCloud;
            arma::mat DistMat;


            //  CrystData::ConstructSRCloud(cloud,basisCollection,combos[i],SrCloud);
            //CrystData::ConstructDistanceMatrix(basisCollection,DistMat);
            ConstructVectorOfDistance(basisCollection, orthogonals, vde);

         //   std::cout << "Constructed vector of distance" << std::endl;

            ConstructSRCloud( cloud,basisCollection,orthogonals, combos[i],SrCloud,index);

         //   std::cout << "Construycted SR cloud" << std::endl;


            double maxRad = ConstructDistanceMatrix(basisCollection,DistMat, index);
            radclouds[index].push_back(maxRad);


            //CrystData::ConstructSRCloud(i,j,cloud,SrCloud);
            //    CrystData::ConstructDistanceMatrix(indices,cloud,DistMat);

            ocd[index].push_back(SrCloud);
            dists[index].push_back(DistMat);
            vds[index].push_back(vde);
//            if (i == 1)
//            {
//                std::cout << "Special prrint" << std::endl;
//
//
//                std::cout << "Anotehrrr" << std::endl;
//
//                PrintMatrix(SrCloudCollection[0]);
//            }


         //!   printVector(combos[i],debugS);

//            if (i == 1)
//            {
//                std::cout << "Special prrint" << std::endl;
//
//
//                std::cout << "Anotehrrr" << std::endl;
//
//                PrintMatrix(SrCloudCollection[0]);
//            }



        }


    }

//    double baseCase(int iii, int jjj,std::vector<int> & A, std::vector<int> & B, std::ofstream & debugS )
//    {
//        double ddd = baseCase(iii, jjj, A, B);
//        debugS << "BaseCase Combination " << iii << " , " << jjj << std::endl;
//        printVectorDouble(A,B, debugS);
//       // printVector(B, debugS);
//        debugS << "Results into: " << ddd << " radius difference is: " << abs(rads[1] - rads[0]) << std::endl;
//        return ddd;
//
//
//    }

    double baseCase(int iii, int jjj, std::vector<int> & A, std::vector<int> & B,  std::ofstream & debugS)
    {

        double distocd;
        std::vector<int> matchd;
        //! storedMatchings[iii][jjj]

        bdcomputatorType bd_computator_tmp;

        //   std::cout << "Dists[0] : " << dists[0].size() << std::endl;
        //    std::cout << "Ocd[0]: "<< ocd[0].size() << std::endl;
        //    std::cout << "Dists[1] : " << dists[1].size() << std::endl;
        //    std::cout << "Ocd[1]: "<< ocd[1].size() << std::endl;

        //   std::cout << dists[0][iii] << std::endl;
        //   std::cout << dists[1][iii] << std::endl;
        //   std::cout << ocd[0][iii] << std::endl;
        //   std::cout << ocd[1][jjj] << std::endl;


        double permute = bd_computator_tmp.computeBetweenInvariantDistancePair(dists[0][iii], dists[1][jjj], ocd[0][iii], ocd[1][jjj], distocd, matchd);
        matchings[std::make_pair(A,B)] = matchd;
       // double tmpdu = abs(rads[1] - rads[0]);
       // double extratmpdu = 0;
     //   extratmpdu = abs(radclouds[0][iii] - radclouds[1][jjj]);
        //  for (int i = 0; i < A.size();i++)
        // {
     //   tmpdu = std::max(tmpdu, extratmpdu);


        //! vds

        double maxxWall = -1;

        for (int i = 0; i < vds[0][iii].size(); i++ )
        {
            double tmpWall = abs(vds[0][iii][i]-vds[1][jjj][i]);
            maxxWall = std::max(maxxWall, tmpWall);



        }

        double tmpdu = maxxWall;


        // }

        //  std::cout << "tmpdu: " << tmpdu << std::endl;

        if (distocd < minVall)
        {
            minVall = distocd;
            selectA = A;
            selectB = B;
            selectAindx = iii;
            selectBindx = jjj;


        }


        return std::max(permute, tmpdu);

    }

    double baseCase(int iii, int jjj, std::vector<int> & A, std::vector<int> & B)
    {

        double distocd;
        std::vector<int> matchd;
        //! storedMatchings[iii][jjj]

        bdcomputatorType bd_computator_tmp;

        //   std::cout << "Dists[0] : " << dists[0].size() << std::endl;
        //    std::cout << "Ocd[0]: "<< ocd[0].size() << std::endl;
        //    std::cout << "Dists[1] : " << dists[1].size() << std::endl;
        //    std::cout << "Ocd[1]: "<< ocd[1].size() << std::endl;

        //   std::cout << dists[0][iii] << std::endl;
        //   std::cout << dists[1][iii] << std::endl;
        //   std::cout << ocd[0][iii] << std::endl;
        //   std::cout << ocd[1][jjj] << std::endl;


        double permute = bd_computator_tmp.computeBetweenInvariantDistancePair(dists[0][iii], dists[1][jjj], ocd[0][iii], ocd[1][jjj], distocd, matchd);
        matchings[std::make_pair(A,B)] = matchd;
       // double tmpdu = abs(rads[1] - rads[0]);
       // double extratmpdu = 0;
     //   extratmpdu = abs(radclouds[0][iii] - radclouds[1][jjj]);
        //  for (int i = 0; i < A.size();i++)
        // {
     //   tmpdu = std::max(tmpdu, extratmpdu);


        //! vds

        double maxxWall = -1;

        for (int i = 0; i < vds[0][iii].size(); i++ )
        {
            double tmpWall = abs(vds[0][iii][i]-vds[1][jjj][i]);
            maxxWall = std::max(maxxWall, tmpWall);



        }

        double tmpdu = maxxWall;


        // }

        //  std::cout << "tmpdu: " << tmpdu << std::endl;

        if (distocd < minVall)
        {
            minVall = distocd;
            selectA = A;
            selectB = B;
            selectAindx = iii;
            selectBindx = jjj;


        }


        return std::max(permute, tmpdu);

    }

    template<typename T>
    arma::Mat<T> convertNestedStdVectorToArmadilloMatrix(const std::vector<std::vector<T>> &V)
    {
        arma::Mat<T> A(V.size(), V[0].size());

        for (int i = 0; i < V.size(); ++i)
        {
            A.row(i) = arma::conv_to< arma::Row<T> >::from(V[i]);
        }
        return A;
    }

    void getFilteredRange(int N, const std::vector<int>& A, std::vector<int>& output)
    {
        std::unordered_set<int> skipSet(A.begin(), A.end()); // For fast lookup
        output.clear(); // Ensure the output vector is empty before filling

        // output.reserve(N - skipSet.size()); // Optional, to reduce reallocations

        for (int i = 0; i < N; ++i)
        {
            if (skipSet.find(i) == skipSet.end())
            {
                output.push_back(i);
            }
        }
    }

    double score(std::vector<std::vector<double>> & distMat, std::vector<int> & A, std::vector<int> & B)
    {
        std::vector<int> matcher;


        double dist = bd_computator.computeForGraph(distMat,matcher);

        matchings[std::make_pair(A,B)] = matcher;





        return dist;

    }



    double score(std::vector<std::vector<double>> & distMat, std::ofstream & debugS, std::vector<int> & A, std::vector<int> & B)
    {
        std::vector<int> matcher;

        arma::mat dd = convertNestedStdVectorToArmadilloMatrix(distMat);

        double dist = bd_computator.computeForGraph(distMat,matcher);

        matchings[std::make_pair(A,B)] = matcher;

        debugS << "Score launched with" << std::endl;
        debugS << dd << std::endl;
        debugS << "Recieved : " << dist << std::endl;






        return dist;

    }

//    void findMatching(std::vector<int> & matching, std::vector<std::vector<int>> & res, int n)
//    {
//        matching = finalMatch;
//    }

    void findMatching(std::vector<int> &matching, std::vector<std::vector<int>> &res, int n)
    {
       // std::cout << "Entering findMatching method, parameter: " << n << std::endl;

        // Resize matching to expected size
        matching.resize(n);

        if (sortedIndices.empty() || sortedIndices[0].size() < 2)
        {
            std::cout << "[ERROR] sortedIndices is empty or has less than 2 elements on row 0.\n";
            return;
        }

        std::vector<int> Aempty;
        if (matchings.find({Aempty, Aempty}) == matchings.end())
        {
            std::cout << "[ERROR] matchings does not contain entry for empty A and B.\n";
            return;
        }

        std::vector<int> outerMatch = matchings[ {Aempty, Aempty}];

        int saveI = sortedIndices[0][0];
        if (saveI < 0 || saveI >= outerMatch.size())
        {
            std::cout << "[ERROR] saveI index " << saveI << " out of bounds for outerMatch of size " << outerMatch.size() << "\n";
            return;
        }

        int saveJ = outerMatch[saveI];

        std::vector<int> Ama{saveI};
        std::vector<int> Bma{saveJ};

        if (matchings.find({Ama, Bma}) == matchings.end())
        {
            std::cout << "[ERROR] matchings does not contain entry for Ama, Bma.\n";
            return;
        }

        std::vector<int> nextLayer = matchings[ {Ama, Bma}];

        std::vector<int> aFiltered;
        std::vector<int> bFiltered;

        getFilteredRange(n, Ama, aFiltered);
        getFilteredRange(n, Bma, bFiltered);

        int findThisGuy = sortedIndices[0][1];
        int foundIt = -1;

        for (int i = 0; i < aFiltered.size(); i++)
        {
            if (aFiltered[i] == findThisGuy)
            {
                foundIt = i;
                break;
            }
        }

        if (foundIt == -1)
        {
            std::cout << "[ERROR] Could not find 'findThisGuy' (" << findThisGuy << ") in aFiltered.\n";
            return;
        }

        if (foundIt >= nextLayer.size())
        {
            std::cout << "[ERROR] foundIt index " << foundIt << " out of bounds for nextLayer of size " << nextLayer.size() << "\n";
            return;
        }

        int correspondance = nextLayer[foundIt];

        std::vector<int> Amb{saveI, findThisGuy};
        std::vector<int> Bmb{saveJ, correspondance};

        if (matchings.find({Amb, Bmb}) == matchings.end())
        {
            std::cout << "[ERROR] matchings does not contain entry for Amb, Bmb.\n";
            return;
        }

        std::vector<int> finalLayer = matchings[ {Amb, Bmb}];

        std::vector<int> aFilteredSuper;
        std::vector<int> bFilteredSuper;

        getFilteredRange(n, Amb, aFilteredSuper);
        getFilteredRange(n, Bmb, bFilteredSuper);

        if (saveI >= matching.size() || findThisGuy >= matching.size())
        {
            std::cout << "[ERROR] Index out of bounds when writing saveI or findThisGuy to matching.\n";
            return;
        }

        matching[saveI] = saveJ;
        matching[findThisGuy] = correspondance;

        for (int i = 0; i < finalLayer.size(); i++)
        {
            if (i >= aFilteredSuper.size() || finalLayer[i] >= bFilteredSuper.size())
            {
                std::cout << "[ERROR] Index out of bounds in finalLayer matching: i = " << i
                          << ", finalLayer[i] = " << finalLayer[i] << "\n";
                continue; // Skip invalid entry
            }
            if (aFilteredSuper[i] >= matching.size())
            {
                std::cout << "[ERROR] aFilteredSuper[" << i << "] = " << aFilteredSuper[i]
                          << " out of bounds for matching size " << matching.size() << "\n";
                continue;
            }

            matching[aFilteredSuper[i]] = bFilteredSuper[finalLayer[i]];
        }

        std::cout << "Resulted into matching: ";
        for (int i = 0; i < matching.size(); i++)
        {
            std::cout << matching[i] << ", ";
        }
        std::cout << "\n";
    }

//    void findMatching(std::vector<int> & matching, std::vector<std::vector<int>> & res, int n)
//    {
//        matching.resize(n);
//        std::cout << "Entering findmatching method , parameter : " << n << std::endl;
//        std::vector<int> Aempty;
//        std::vector<int> outerMatch= matchings[std::make_pair(Aempty,Aempty)];
//
//        int saveI = sortedIndices[0][0];
//        int saveJ = outerMatch[saveI];
//
//        std::vector<int> Ama{saveI};
//        std::vector<int> Bma{saveJ};
//        std::vector<int> nextLayer = matchings[std::make_pair(Ama,Bma)];
//
//
//        std::vector<int> aFiltered;
//        std::vector<int> bFiltered;
//
//        getFilteredRange(n,Ama,aFiltered);
//        getFilteredRange(n,Bma,bFiltered);
//
//
//        int findThisGuy = sortedIndices[0][1];
//
//        int foundIt = 0;
//        for (int i = 0 ; i < aFiltered.size() ;  i++ )
//        {
//            if (aFiltered[i] == findThisGuy)
//            {
//                foundIt = i;
//                break;
//            }
//        }
//        int correspondance = nextLayer[foundIt];
//
//        std::vector<int> Amb{saveI, findThisGuy };
//        std::vector<int> Bmb{saveJ, correspondance };
//
//        std::vector<int> finalLayer = matchings[std::make_pair(Amb,Bmb)];
//
//
//        std::vector<int> aFilteredSuper;
//        std::vector<int> bFilteredSuper;
//
//        getFilteredRange(n,Amb,aFilteredSuper);
//        getFilteredRange(n,Bmb,bFilteredSuper);
//
//
//        matching[saveI] = saveJ;
//        matching[findThisGuy] = correspondance;
//
//        for (int i =0 ; i < finalLayer.size(); i++)
//        {
//            matching[aFilteredSuper[i]] = bFilteredSuper[finalLayer[i]];
//
//
//        }
//
//        std::cout << "Resulted into: " << std::endl;
//
//        for (int i = 0; i < matching.size() ; i++)
//        {
//
//            std::cout << matching[i] << ", ";
//
//        }
//
//    }

    //! and the rest ...












    //   selectAindx = saveI;
    //   selectBindx = saveJ;






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


    //!  matching = finalMatch;





    virtual ~PBI_VD()
    {


        //     *debugS.close();
        //    delete debugS;

    }

protected:

private:
    // PBI_Update(const PBI_Update&);
    //  PBI_Update& operator=(const PBI_Update&);
    // bool debug;
    //  std::ofstream* debugS;

    boost::unordered::unordered_map<std::pair<std::vector<int>, std::vector<int>>,std::vector<int>> matchings;
    double aConst;
    double constC;
    bdcomputatorType bd_computator;
    std::vector<double> rads;
    std::vector<std::vector<size_t>> sortedIndices;
    std::vector<std::vector<std::vector<double>>> vds;
    std::vector<std::vector<arma::mat>> ocd;
    std::vector<std::vector<arma::mat>> dists;
    std::vector<std::vector<double>> radclouds;
    //  std::vector<std::vector<std::vector<int>>> storedMatchings;
    std::vector<std::vector<double>> storedDists;
    std::vector<int> furthestPoint;
    std::vector<int> finalMatch;
    double minVall;
    std::vector<int> selectA;
    std::vector<int> selectB;
    int selectAindx;
    int selectBindx;
};


#endif // PBI_VD_H
