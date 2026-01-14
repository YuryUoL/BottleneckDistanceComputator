#ifndef EUCLIDEANBDCOMPUTE_H
#define EUCLIDEANBDCOMPUTE_H

template<typename graph_compute_engine_type>
class EuclideanBDCompute
{
public:
    EuclideanBDCompute() {}
    virtual ~EuclideanBDCompute() {}

    double compute(arma::mat & A, arma::mat & B, std::vector<int> & matching)
    {
        std::vector<std::vector<double>> graph(A.n_cols, std::vector<double>(B.n_cols));
        for (int i = 0; i < A.n_cols; i++)
        {
            for (int j = 0; j < B.n_cols; j++)
            {
                graph[i][j] = arma::norm(A.col(i) - B.col(j),2);
            }
        }

        double bda = gt.bottleneck_on_bipartie_graph(graph, matching);

        return bda;






    }

protected:

private:
    graph_compute_engine_type gt;

};

#endif // EUCLIDEANBDCOMPUTE_H
