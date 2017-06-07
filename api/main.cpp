/*
 * main.cpp
 *
 *  Created on: 20 Apr 2017
 *      Author: joni
 */

#include "FalconnSearch.h"

#include <iostream>
#include <random>

FalconnSearch fs;

int main(int argc, char** argv) {
    
    std::cout << "Test" << std::endl;
    
    std::vector<Point> dataset;
    
    // generate random data for index
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> cluster1(0.0f, 1.0f);
    for(int i = 0; i < 10000; i++) {
        Point feature;
        feature.resize(102);
        for(int j = 0; j < 102; j++) {
            feature[j] = (cluster1(gen));
        }
        dataset.push_back(feature);
    }
    
    
    int result = fs.buildIndexFromData(dataset);
    std::cout << "Index build result: " << result << std::endl;
    
    std::vector<int> searchResult = fs.search(dataset[123]);
    std::cout << "Serch result size: " << searchResult.size() << std::endl;
    
}
