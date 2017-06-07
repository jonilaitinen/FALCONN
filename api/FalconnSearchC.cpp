/*
 * FalconnSearchC.cpp
 *
 *  Created on: 20 Apr 2017
 *      Author: joni
 */

#include "FalconnSearchC.h"

#include <iostream>
#include <vector>

#include <falconn/lsh_nn_table.h>
#include <falconn/falconn_global.h>

typedef falconn::DenseVector<float> Point;

const int NUM_QUERIES = 100;
const int SEED = 4057218;
const int NUM_HASH_TABLES = 50;
const int NUM_HASH_BITS = 18;
const int NUM_ROTATIONS = 1;

std::unique_ptr<falconn::LSHNearestNeighborTable<Point, int32_t> > table;

Point testQuery;

float bytesToFloat(unsigned char b0, unsigned char b1, unsigned char b2, unsigned char b3) {
    unsigned char byte_array[] = { b0, b1, b2, b3 }; // change byte order here for big/little endian
    float result;
    std::copy(reinterpret_cast<const char*>(&byte_array[0]),
              reinterpret_cast<const char*>(&byte_array[4]),
              reinterpret_cast<char*>(&result));
    
    return result;
}

std::vector<Point> createDataVector(void* dataPtr,
		int featureCount, int featureBytes) {
	unsigned char* data = static_cast<unsigned char*>(dataPtr);
	int dataLength = featureCount * featureBytes;

    std::vector<Point> dataset;

	for (int i = 0; i < dataLength; i += featureBytes) {
		Point featureVector;
        int dimensions = featureBytes / 4;
        featureVector.resize(dimensions);
        int dimension = 0;
		for (int j = 0; j < featureBytes; j += 4) {
			float value = bytesToFloat(data[i + j], data[i + j + 1],
					data[i + j + 2], data[i + j + 3]);
			featureVector[dimension] = value;
            dimension++;
		}
		dataset.push_back(featureVector);
	}
	return dataset;
}

void normalize(std::vector<Point> *dataset) {
    for (auto &p : *dataset) {
        p.normalize();
    }
}

/*
 * Chooses a random subset of the dataset to be the queries. The queries are
 * taken out of the dataset.
 */
void gen_queries(std::vector<Point> *dataset, std::vector<Point> *queries) {
    std::mt19937_64 gen(SEED);
    queries->clear();
    for (int i = 0; i < NUM_QUERIES; ++i) {
        std::uniform_int_distribution<> u(0, dataset->size() - 1);
        int ind = u(gen);
        queries->push_back((*dataset)[ind]);
        (*dataset)[ind] = dataset->back();
        dataset->pop_back();
    }
}

/*
 * Generates answers for the queries using the (optimized) linear scan.
 */
void gen_answers(const std::vector<Point> &dataset, const std::vector<Point> &queries,
                 std::vector<int> *answers) {
    answers->resize(queries.size());
    int outer_counter = 0;
    for (const auto &query : queries) {
        float best = -10.0;
        int inner_counter = 0;
        for (const auto &datapoint : dataset) {
            float score = query.dot(datapoint);
            if (score > best) {
                (*answers)[outer_counter] = inner_counter;
                best = score;
            }
            ++inner_counter;
        }
        ++outer_counter;
    }
}

/*
 * Computes the probability of success using a given number of probes.
 */
double evaluate_num_probes(falconn::LSHNearestNeighborTable<Point> *table,
                           const std::vector<Point> &queries,
                           const std::vector<int> &answers, int num_probes) {
    table->set_num_probes(num_probes);
    int outer_counter = 0;
    int num_matches = 0;
    std::vector<int32_t> candidates;
    for (const auto &query : queries) {
        table->get_candidates_with_duplicates(query, &candidates);
        for (auto x : candidates) {
            if (x == answers[outer_counter]) {
                ++num_matches;
                break;
            }
        }
        ++outer_counter;
    }
    return (num_matches + 0.0) / (queries.size() + 0.0);
}

/*
 * Finds the smallest number of probes that gives the probability of success
 * at least 0.9 using binary search.
 */
int find_num_probes(falconn::LSHNearestNeighborTable<Point> *table,
                    const std::vector<Point> &queries, const std::vector<int> &answers,
                    int start_num_probes) {
    int num_probes = start_num_probes;
    for (;;) {
        std::cout << "trying " << num_probes << " probes" << std::endl;
        double precision = evaluate_num_probes(table, queries, answers, num_probes);
        if (precision >= 0.9) {
            break;
        }
        num_probes *= 2;
    }
    
    int r = num_probes;
    int l = r / 2;
    
    while (r - l > 1) {
        int num_probes = (l + r) / 2;
        std::cout << "trying " << num_probes << " probes" << std::endl;
        double precision = evaluate_num_probes(table, queries, answers, num_probes);
        if (precision >= 0.9) {
            r = num_probes;
        } else {
            l = num_probes;
        }
    }
    
    return r;
}

int buildIndexFromData(void* dataPtr, int featureCount, int featureBytes) {
    
    std::vector<Point> dataset =
    		createDataVector(dataPtr, featureCount, featureBytes);
    
    std::vector<Point> queries;
    std::vector<int> answers;
    
    // normalize the data points
    std::cout << "normalizing points" << std::endl;
    normalize(&dataset);
    std::cout << "done" << std::endl;
    
    // find the center of mass
    Point center = dataset[0];
    for (size_t i = 1; i < dataset.size(); ++i) {
        center += dataset[i];
    }
    center /= dataset.size();
    
    // selecting NUM_QUERIES data points as queries
    std::cout << "selecting " << NUM_QUERIES << " queries" << std::endl;
    gen_queries(&dataset, &queries);
    std::cout << "done" << std::endl;
    
    // running the linear scan
    std::cout << "running linear scan (to generate nearest neighbors)" << std::endl;
    auto t1 = std::chrono::high_resolution_clock::now();
    gen_answers(dataset, queries, &answers);
    auto t2 = std::chrono::high_resolution_clock::now();
    double elapsed_time = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count();
    std::cout << "done" << std::endl;
    std::cout << elapsed_time / queries.size() << " s per query" << std::endl;
    
    // re-centering the data to make it more isotropic
    std::cout << "re-centering" << std::endl;
    for (auto &datapoint : dataset) {
        datapoint -= center;
    }
    for (auto &query : queries) {
        query -= center;
    }
    std::cout << "done" << std::endl;
    
    // setting parameters and constructing the table
    falconn::LSHConstructionParameters params;
    params.dimension = dataset[0].size();
    params.lsh_family = falconn::LSHFamily::CrossPolytope;
    params.l = NUM_HASH_TABLES;
    params.distance_function = falconn::DistanceFunction::EuclideanSquared;
    falconn::compute_number_of_hash_functions<Point>(NUM_HASH_BITS, &params);
    params.num_rotations = NUM_ROTATIONS;
    // we want to use all the available threads to set up
    params.num_setup_threads = 0;
    params.storage_hash_table = falconn::StorageHashTable::BitPackedFlatHashTable;
    /*
     For an easy way out, you could have used the following.
     LSHConstructionParameters params
     = get_default_parameters<Point>(dataset.size(),
     dataset[0].size(),
     DistanceFunction::EuclideanSquared,
     true);
     */
    std::cout << "building the index based on the cross-polytope LSH" << std::endl;
    t1 = std::chrono::high_resolution_clock::now();
    auto lshTable = falconn::construct_table<Point>(dataset, params);
    table = std::move(lshTable);
    t2 = std::chrono::high_resolution_clock::now();
    elapsed_time = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count();
    std::cout << "done" << std::endl;
    std::cout << "construction time: " << elapsed_time << std::endl;
    
    // finding the number of probes via the binary search
    std::cout << "finding the appropriate number of probes" << std::endl;
    int num_probes = find_num_probes(&*table, queries, answers, params.l);
    std::cout << "done" << std::endl;
    std::cout << num_probes << " probes" << std::endl;
    
    std::cout << "Start search.. table: " << table.get() << std::endl;
    std::vector<int> knearest;
    Point testQuery = dataset[123];
    testQuery.normalize();
    table->find_k_nearest_neighbors(testQuery, 30, &knearest);
    std::cout << "Search complete, result size: " << knearest.size() << std::endl;
    
    return 1;
}

Result search(void* dataPtr, int dataLength) {
    
    std::cout << "Start search.. " << std::endl;
    
    Result result;
    Point query;
    
    unsigned char* data = static_cast<unsigned char*>(dataPtr);
    
    int dimensions = dataLength / 4;
    query.resize(dimensions);
    int dimension = 0;
    for(int i = 0; i < dataLength; i += 4) {
        float value = bytesToFloat(data[i], data[i + 1], data[i + 2], data[i + 3]);
        query[dimension] = value;
        std::cout << value << ", ";
        dimension++;
    }
    std::cout << std::endl;
    std::cout << "search dimensions: " << dimension << ", table: " << table.get() << std::endl;
    
    auto statistics = table->get_query_statistics();
    std::cout << "average total query time: " << statistics.average_total_query_time
    << std::endl;
    std::cout << "average lsh time: " << statistics.average_lsh_time << std::endl;
    std::cout << "average hash table time: " << statistics.average_hash_table_time
    << std::endl;
    std::cout << "average distance time: " << statistics.average_distance_time
    << std::endl;
    std::cout << "average number of candidates: "
    << statistics.average_num_candidates << std::endl;
    std::cout << "average number of unique candidates: "
    << statistics.average_num_unique_candidates << std::endl;
    
    query.normalize();
    std::vector<int> knearest;
    table->find_k_nearest_neighbors(testQuery, 30, &knearest);
    
    std::cout << "Search complete, result size: " << knearest.size() << std::endl;
    
    result.Result = new int[knearest.size()];
    result.ResultCount = knearest.size();
    
    for(int i = 0; i < knearest.size(); ++i) {
        result.Result[i] = knearest[i];
    }
    
    return result;
}
