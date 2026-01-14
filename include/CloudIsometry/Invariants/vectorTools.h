#ifndef VECTORTOOLS_H
#define VECTORTOOLS_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

//#include "experiments/CrystData.h"

class vectorTools
{
public:
    vectorTools() {}

    static void mat_to_std_vec(arma::mat &A, std::vector<std::vector<double>> & V)
    {
        arma::mat B;
        B = A.t();
        V = std::vector<std::vector<double>>(B.n_rows);
        for (size_t i = 0; i < B.n_rows; ++i)
        {
            V[i] = arma::conv_to< std::vector<double> >::from(B.row(i));
        };
        //std::cout << "Number of points: " << V.size() << std::endl;
       // std::cout << "Claimed dimension: " << V[0].size() << std::endl;
       // CrystData::PrintMatrix(V);
       // return V;
    }
    static void vec_to_mat(std::vector<std::vector<double>> & V, arma::mat & A )
    {
        A = arma::mat(V[0].size(), V.size());
        for (int i = 0; i < V.size();i++)
        {
            arma::vec p(V[i]);
            A.col(i) = p;
        }



    }

    static void multiDcrossProduct(std::vector<std::vector<double>> & vectors, std::vector<double> & c_P)
    {
        int dimm = vectors[0].size();
        if (dimm == 3)
        {
            vectorTools::crossProduct(vectors[0], vectors[1], c_P);

        }
        else if (dimm == 2)
        {
            double a = -vectors[0][1];
            double b = vectors[0][0];
            c_P = std::vector<double> {a,b};
        }
    }
    static double linfinity(std::vector<double> & v_A, std::vector<double> & v_B)
    {
        double maxx = 0;
        for (int i = 0; i < v_A.size(); i++)
        {
            double tmp = abs(v_A[i] - v_B[i]);
            maxx = std::max(tmp, maxx);

        }
        return maxx;

    }

     static double euclideanDistance(std::vector<double> & v_A, std::vector<double> & v_B)
    {
        double maxx = 0;
        for (int i = 0; i < v_A.size(); i++)
        {
            double tmp = (v_A[i] - v_B[i]) * (v_A[i] - v_B[i]);
            maxx = maxx + tmp;

        }
        return sqrt(maxx);

    }

    static void crossProduct(std::vector<double> & v_A, std::vector<double> & v_B, std::vector<double> & c_P)
    {
        c_P[0] = v_A[1] * v_B[2] - v_A[2] * v_B[1];
        c_P[1] = -(v_A[0] * v_B[2] - v_A[2] * v_B[0]);
        c_P[2] = v_A[0] * v_B[1] - v_A[1] * v_B[0];
    }
    static double dotProduct(std::vector<double> & v_A, std::vector<double> & v_B)
    {
        double sum = 0;
        for (int i= 0 ; i < v_A.size(); i++)
        {
            sum = sum + v_A[i] * v_B[i];
        }
        return sum;
    }
    static void vectorAddition(std::vector<double> & vectorA, std::vector<double> & vectorB)
    {
        for (int i = 0; i < vectorA.size(); i++)
        {
            vectorA[i] = vectorA[i] + vectorB[i];
        }

    }
    static void vectorScalarMultiplication(std::vector<double> & vectorA, double scalar)
    {
        for (int i = 0; i < vectorA.size(); i++)
        {
            vectorA[i] = scalar * vectorA[i];
        }

    }
    static void MatrixScalarMultiplication(std::vector<std::vector<double>> & matrixA, double scalar)
    {
        for (int i = 0; i < matrixA.size(); i++)
        {
            for (int j = 0; j < matrixA[i].size(); j++)
            {
                matrixA[i][j] = scalar * matrixA[i][j];
            }
        }
    }



    static double permutationSign(std::vector<int> & a)
    {
        int n = a.size();
        int cnt=0;
        for(int i=0; i<n; i++)
            for(int j=i+1; j<n; j++)
                if (a[i]>a[j]) cnt++;

        if (cnt%2 == 0)
        {
            return 1;
        }
        else
        {
            return -1;

        }

    }

    static void MatrixScalarMultiplicationSpecial(std::vector<std::vector<double>> & matrixA, double scalar)
    {
        for (int i = 0; i < matrixA.size(); i++)
        {
            for (int j = 0; j < matrixA[i].size()-1; j++)
            {
                matrixA[i][j] = scalar * matrixA[i][j];
            }
            double dimensionPower = constants::Dim-1;
            double newscalar = pow(scalar,dimensionPower);
            //std::cout << "Scalar: " << scalar << " , " << "Power: " << dimensionPower << std::endl;
            matrixA[i][matrixA[i].size()-1] = newscalar * matrixA[i][matrixA[i].size()-1];
        }

    }



    static void PermuteRowsAndMultiplyLastColumn(std::vector<std::vector<double>> & In, std::vector<std::vector<double>> & Out, double scalar, std::vector<int> & permutation)
    {
        Out = In;
        for (int i = 0; i < Out.size(); i++)
        {
            for (int j = 0; j < permutation.size(); j++ )
            {
                Out[i][permutation[j]] = In[i][j];
            }
            Out[i][Out[i].size()-1] = scalar * Out[i][Out[i].size()-1];
        }
    }

    static double SquaredNorm(std::vector<double> & vectorNorm)
    {
        double sum = 0;
        for (int i = 0; i < vectorNorm.size(); i++)
        {
            sum = sum + vectorNorm[i] * vectorNorm[i];

        }
        return sum;
    }

    static double euclideanDistanceBetweenVectors(std::vector<double> & vectorA, std::vector<double> & vectorB)
    {
        std::vector<double> A = vectorA;
        std::vector<double> B = vectorA;
        vectorTools::vectorScalarMultiplication(B, -1);
        vectorTools::vectorAddition(A,B);
        return sqrt(SquaredNorm(A));
    }

    static void permuteMatrix(std::vector<std::vector<double>> & in, std::vector<std::vector<double>> & out, std::vector<int> & permutation )
    {
        out = in;
        for (int i = 0; i < in.size(); i++)
        {
            for (int j = 0; j < in[i].size(); j++)
            {
                // std::cout << "Iterating over " << i << " and " << j << std::endl;
                // std::cout << "Gives " << permutation[i] << " and " << permutation[j] << std::endl;
                out[permutation[i]][permutation[j]] = in[i][j];
            }
        }
    }


protected:

private:
};

#endif // VECTORTOOLS_H
