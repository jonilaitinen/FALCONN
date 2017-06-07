
#include <vector>

#include <falconn/lsh_nn_table.h>
#include <falconn/falconn_global.h>

typedef falconn::DenseVector<float> Point;

class FalconnSearch {
public:
    FalconnSearch();
    virtual ~FalconnSearch() {}
    
    int buildIndexFromData(std::vector<Point> dataset);
    std::vector<int> search(Point query);
    
private:
    void gen_queries(std::vector<Point> *dataset, std::vector<Point> *queries);
    void gen_answers(const std::vector<Point> &dataset, const std::vector<Point> &queries,
                     std::vector<int> *answers);
    double evaluate_num_probes(falconn::LSHNearestNeighborTable<Point> *table,
                               const std::vector<Point> &queries,
                               const std::vector<int> &answers, int num_probes);
    int find_num_probes(falconn::LSHNearestNeighborTable<Point> *table,
                        const std::vector<Point> &queries, const std::vector<int> &answers,
                        int start_num_probes);
    void normalize(std::vector<Point> *dataset);
    
    std::unique_ptr<falconn::LSHNearestNeighborTable<Point, int32_t> > table;
    
};
