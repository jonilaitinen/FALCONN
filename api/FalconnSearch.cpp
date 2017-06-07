/*
 * FalconnSearch.cpp
 *
 *  Created on: 20 Apr 2017
 *      Author: joni
 */

#include "FalconnSearch.h"

#include <iostream>
#include <thread>

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

int FalconnSearch::buildIndexFromData(std::vector<Point> dataset) {
    
    std::thread::id this_id = std::this_thread::get_id();
    std::cout << "Build thread id: " << this_id << std::endl;
    
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
    
    std::cout << "Start search.. table: " << table.get() << std::endl;
    std::vector<int> knearest;
    Point testQuery = dataset[123];
    testQuery.normalize();
    table->find_k_nearest_neighbors(testQuery, 30, &knearest);
    std::cout << "Search complete, result size: " << knearest.size() << std::endl;
    
    return 1;
}

std::vector<int> FalconnSearch::search(Point query) {
    
    std::thread::id this_id = std::this_thread::get_id();
    std::cout << "Search thread id: " << this_id << std::endl;
    
    query.normalize();
    std::vector<int> knearest;
    table->find_k_nearest_neighbors(query, 30, &knearest);
    
    std::cout << "Search complete, result size: " << knearest.size() << std::endl;
    
    return knearest;
}
