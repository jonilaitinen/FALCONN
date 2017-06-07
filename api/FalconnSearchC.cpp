/*
 * FalconnSearchC.cpp
 *
 *  Created on: 20 Apr 2017
 *      Author: joni
 */

#include "FalconnSearchC.h"

#include <iostream>
#include <vector>

#include "FalconnSearch.h"

FalconnSearch falconnSearch;

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

int buildIndexFromData(void* dataPtr, int featureCount, int featureBytes) {
    
    std::vector<Point> dataset =
    		createDataVector(dataPtr, featureCount, featureBytes);
    
    return falconnSearch.buildIndexFromData(dataset);
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
        dimension++;
    }
    std::cout << std::endl;
    
    std::vector<int> knearest = falconnSearch.search(query);
    
    result.Result = new int[knearest.size()];
    result.ResultCount = knearest.size();
    
    for(int i = 0; i < knearest.size(); ++i) {
        result.Result[i] = knearest[i];
    }
    
    return result;
}
