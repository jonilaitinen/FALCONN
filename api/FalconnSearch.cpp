/*
 * FalconnSearch.cpp
 *
 *  Created on: 20 Apr 2017
 *      Author: joni
 */

#include "FalconnSearch.h"

#include <iostream>
#include <map>

const int NUM_QUERIES = 100;
const int SEED = 4057218;
const int NUM_HASH_TABLES = 50;
const int NUM_HASH_BITS = 18;
const int NUM_ROTATIONS = 1;

FalconnSearch::FalconnSearch() {
    
}

/*
 * Chooses a random subset of the dataset to be the queries. The queries are
 * taken out of the dataset.
 */
/*
void FalconnSearch::gen_queries(std::vector<Point> *dataset, std::vector<Point> *queries) {
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
*/

/*
 * Generate queries by randomly choosing points in the dataset and taking an average of
 * 8 points.
 */
void FalconnSearch::gen_queries(std::vector<Point> *dataset, std::vector<Point> *queries) {
    std::mt19937_64 gen(SEED);
    queries->clear();
    for (int i = 0; i < NUM_QUERIES; ++i) {
        std::uniform_int_distribution<> u(0, dataset->size() - 8);
        int ind = u(gen);
        Point query1 = (*dataset)[ind];
        Point query2 = (*dataset)[ind + 1];
        Point query3 = (*dataset)[ind + 2];
        Point query4 = (*dataset)[ind + 3];
        Point query5 = (*dataset)[ind + 4];
        Point query6 = (*dataset)[ind + 5];
        Point query7 = (*dataset)[ind + 6];
        Point query8 = (*dataset)[ind + 7];
        
        Point query = query1 + query2 + query3 + query4 + query5 + query6 + query7 + query8;
        query /= 8;
        
        queries->push_back(query);
    }
}


/*
 * Generates answers for the queries using the (optimized) linear scan.
 */
void FalconnSearch::gen_answers(const std::vector<Point> &dataset, const std::vector<Point> &queries,
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
double FalconnSearch::evaluate_num_probes(falconn::LSHNearestNeighborTable<Point> *table,
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
int FalconnSearch::find_num_probes(falconn::LSHNearestNeighborTable<Point> *table,
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

void FalconnSearch::normalize(std::vector<Point> *dataset) {
    for (auto &p : *dataset) {
        p.normalize();
    }
}

int FalconnSearch::buildIndexFromData(const std::vector<Point>& inDataset) {
    
    dataset = inDataset;
    std::vector<Point> queries;
    std::vector<int> answers;
    
    // normalize the data points
    std::cout << "normalizing points" << std::endl;
    normalize(&dataset);
    std::cout << "done" << std::endl;
    
    // find the center of mass
    center = dataset[0];
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
    table = falconn::construct_table<Point>(dataset, params);
    //table = std::move(lshTable);
    t2 = std::chrono::high_resolution_clock::now();
    elapsed_time = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count();
    std::cout << "done" << std::endl;
    std::cout << "construction time: " << elapsed_time << std::endl;
    
    // finding the number of probes via the binary search
    std::cout << "finding the appropriate number of probes" << std::endl;
    int num_probes = find_num_probes(&*table, queries, answers, params.l);
    std::cout << "done" << std::endl;
    std::cout << num_probes << " probes" << std::endl;
    
    return 1;
}

std::vector<int> FalconnSearch::search(Point query) {
    
    std::cout << "raw input: " << std::endl;
    for(int i = 0; i < 10; i++) {
        std::cout << query[i] << ", ";
    }
    std::cout << std::endl;
    std::cout << std::endl;
    
    query.normalize();
    
    std::cout << "normalized input: " << std::endl;
    for(int i = 0; i < 10; i++) {
        std::cout << query[i] << ", ";
    }
    std::cout << std::endl;
    std::cout << std::endl;
    
    query -= center;
    
    std::cout << "centered input: " << std::endl;
    for(int i = 0; i < 10; i++) {
        std::cout << query[i] << ", ";
    }
    std::cout << std::endl;
    std::cout << "center: " << std::endl;
    for(int i = 0; i < 10; i++) {
        std::cout << center[i] << ", ";
    }
    std::cout << std::endl;
    
    std::vector<int> knearest;
    table->find_k_nearest_neighbors(query, 30, &knearest);
    
    std::cout << "Search complete, result size: " << knearest.size() << std::endl;
    
    std::map<float, int> distanceMap;
    for(int i = 0; i < dataset.size(); i++) {
        float dist = l2_distance(query, dataset[i]);
        distanceMap.insert(std::pair<float, int>(dist, i));
    }
    
    std::cout << "Nearest results: " << std::endl;
    auto it = distanceMap.begin();
    for(int i = 0; i < knearest.size(); i++) {
        float diff = l2_distance(query, dataset[knearest[i]]);
        std::cout << knearest[i]
        << ", \tdiff: " << diff
        << ", \ttrue nearest: " << it->second
        << ", dist: " << it->first << std::endl;
        it++;
    }

    
    return knearest;
}

float FalconnSearch::l2_distance(Point p1, Point p2) {
    float sum = 0;
    for(int i = 0; i < p1.size(); i++) {
        sum += pow(p1[i] - p2[i], 2);
    }
    return sqrt(sum);
}

/*
void testAccuracy() {
    
    std::cout << "Start search.. " << std::endl;
    std::vector<int> knearest;
    Point testQuery = dataset[1];
    table->find_k_nearest_neighbors(testQuery, 30, &knearest);
    std::cout << "Search complete, result size: " << knearest.size() << std::endl;
    
    std::cout << "Test query dist: " << l2_distance(testQuery, testQuery) << std::endl;
    std::cout << "Datapoint dist: " << l2_distance(dataset[1], dataset[1]) << std::endl;
    std::cout << "Test to Datapoint dist: " << l2_distance(testQuery, dataset[1]) << std::endl;
    
    std::map<float, int> distanceMap;
    for(int i = 0; i < dataset.size(); i++) {
        float dist = l2_distance(testQuery, dataset[i]);
        distanceMap.insert(std::pair<float, int>(dist, i));
    }
    
    std::cout << "Nearest results: " << std::endl;
    auto it = distanceMap.begin();
    for(int i = 0; i < knearest.size(); i++) {
        float diff = l2_distance(testQuery, dataset[knearest[i]]);
        std::cout << knearest[i]
        << ", \tdiff: " << diff
        << ", \ttrue nearest: " << it->second
        << ", dist: " << it->first << std::endl;
        it++;
    }
}
*/

