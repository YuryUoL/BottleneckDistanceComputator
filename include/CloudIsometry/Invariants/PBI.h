#ifndef PBI_H
#define PBI_H


class PBI
{
public:
    PBI() {}
    virtual ~PBI() {}

    int NumElem(arma::mat & cloudAunc)
    {
        return cloudAunc.n_cols;

    }
    int Dimension(arma::mat & cloudAunc)
    {

        return cloudAunc.n_rows;

    }

    double initilize(arma::mat & cloudAunc, arma::mat & cloudA)
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

        }

        cloudA = cloudOUT;

        return radius;
    }

    arma::vec MultiDCrossProduct(arma::mat & cloud, std::vector<int> & indices)
    {
        if (cloud.n_rows == 3)
        {
            arma::vec c_P = cross( cloud.col(indices[0]), cloud.col(indices[1]) );
            return c_P;
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

    void computeBasisMatrix(std::vector<int> & indices,arma::mat & cloud, arma::mat & basisCollection)
    {

        basisCollection = arma::mat(cloud.n_rows, indices.size() + 1);
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

        arma::vec c_P = MultiDCrossProduct(cloud, indices);
        basisCollection.col(indices.size()) = c_P;





    }

    void ConstructSRCloud( arma::mat & cloud, arma::mat & basis,  std::vector<int> & indices,arma::mat & SrCloud)
    {
        int m = cloud.n_cols;
        int n = cloud.n_rows;
        //!arma::vec c_P(cloud.n_rows);
        SrCloud = arma::mat(n, m - n + 1);
        //  std::cout << "SrCloud rows: " << SrCloud.n_rows << " , SrCloud columns: " << SrCloud.n_cols << std::endl;
        std::set<int> s(indices.begin(), indices.end());

        // std::cout << "First: " << cloud[i][0] << ", "<< cloud[i][1] << ", " << cloud[i][2] << std::endl;
        // std::cout << "Second: " << cloud[j][0] << ", "<< cloud[j][1] << ", " << cloud[j][2] << std::endl;
//        vectorTools::crossProduct(cloud[i], cloud[j], c_P);
        // std::cout << "Final: " << c_P[0] << ", "<< c_P[1] << ", " << c_P[2] << std::endl;
        int summ = 0;
        for (int k = 0; k < cloud.n_cols ; k++)
        {
            // std::cout << "Iteration, k : " << k << " summ : " << summ << std::endl;
            if (s.find(k) == s.end())
            {
                std::vector<double> newVector;
                for (int i = 0; i < basis.n_cols; i++)
                {
                    double dottR = arma::dot(cloud.col(k), basis.col(i));
                    newVector.push_back(dottR);
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

    void ConstructDistanceMatrix(arma::mat & basisCollection, arma::mat & matrix)
    {
        //!! Assume that dimension = 3
        matrix = arma::mat(basisCollection.n_cols, basisCollection.n_cols);

        for (int i = 0; i < basisCollection.n_cols; i++)
        {
            for (int j = 0; j < basisCollection.n_cols; j++)
            {
                double ttt = arma::norm(basisCollection.col(i) - basisCollection.col(j), 2);

                matrix(i,j) = ttt * ttt;

//                std::vector<double> vecU = basisCollection[j];
//                vectorTools::vectorScalarMultiplication(vecU,-1);
//                std::vector<double> vecV = basisCollection[i];
//                vectorTools::vectorAddition(vecV, vecU);
//              //  double norm;
//                matrix[i][j] = vectorTools::SquaredNorm(vecV);
            }
        }
    }




    void ConstructCloudInvariantCombo(arma::mat & cloud, std::vector<std::vector<int>> & combos, std::vector<arma::mat> &  OCDCloudCollection, std::vector<arma::mat> & distMatrixCollection)
    {

        //  std::cout << "We are entering ConstructionCloudInvariantCombo with size: " << combos.size() << std::endl;
        for (int i = 0; i < combos.size(); i++)
        {
            arma::mat basisCollection;
            //! We are here:
            //!    CrystData::computeBasisMatrix(combos[i],cloud,basisCollection);
            computeBasisMatrix(combos[i],cloud, basisCollection);
            arma::mat SrCloud;
            arma::mat DistMat;


            //  CrystData::ConstructSRCloud(cloud,basisCollection,combos[i],SrCloud);
            //CrystData::ConstructDistanceMatrix(basisCollection,DistMat);
            ConstructSRCloud(  cloud,basisCollection, combos[i],SrCloud);


            ConstructDistanceMatrix(basisCollection,DistMat);


            //CrystData::ConstructSRCloud(i,j,cloud,SrCloud);
            //    CrystData::ConstructDistanceMatrix(indices,cloud,DistMat);

            OCDCloudCollection.push_back(SrCloud);
            distMatrixCollection.push_back(DistMat);
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



protected:

private:
};

#endif // PBI_H
