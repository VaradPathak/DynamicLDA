/*
 * TopicChains.cpp
 *
 *  Created on: Apr 22, 2014
 *      Author: vspathak
 */

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/detail/adjacency_list.hpp>
#include <boost/graph/graph_selectors.hpp>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using namespace std;
using namespace boost;

typedef adjacency_list<vecS, vecS, undirectedS> Graph;

// Initialize number of documents, topics and words in vocabulary
unsigned int W, D, K;

double KLDivergence(double*** Pi, int t, int k, double* M) {
	double result = 0.0;
	for (unsigned int w = 0; w < W; ++w) {
		result += log(Pi[t][w][k] / M[w]) * Pi[t][w][k];
	}
	return result;
}

double JSsimilarity(double*** Pi, int t1, int k1, int t2, int k2) {
	double result = 0.0;
	double* M = new double[W];
	for (unsigned int w = 0; w < W; ++w) {
		M[w] = (Pi[t1][w][k1] + Pi[t2][w][k2]) / 2;
	}
	result = KLDivergence(Pi, t1, k1, M) + KLDivergence(Pi, t2, k2, M);
	result = result / 2;
	return result;
}

void generateTopicLinks(Graph &G, double*** Pi, int timeSlice, int topic,
		int numTopics, int windowSize, double threshold) {
	for (int w = 0; w < windowSize; w++) {
		int numLinks = 0;
		for (int k = 0; k < numTopics; k++) {
			if ((timeSlice - 1 - w >= 0) && JSsimilarity(Pi, timeSlice, topic, timeSlice - 1 - w, k) > threshold) {
				//add edge to graph structure here
				int e1 = (timeSlice * numTopics) + topic;
				int e2 = ((timeSlice - 1 - w) * numTopics) + k;

				cout << "Adding edge " << e1 << ", " << e2 << endl;
				add_edge(e1, e2, G);
				numLinks++;
			}
		}
		if (numLinks > 0) {
			break;
		}
	}
}

void generateAllLinks(Graph &G, double*** Pi, int numTimeSlices, int numTopics,
		int windowSize, double threshold) {
	for (int t = 0; t < numTimeSlices; t++) {
		for (int k = 0; k < numTopics; k++) {
			generateTopicLinks(G, Pi, t, k, numTopics, windowSize, threshold);
		}
	}
}

int main(int argc, char* argv[]) {
	if (argc < 4) {
		printf("Usage: ./fastLDA Pi_folder num_topics WindowSize SimilarityThreshold\n");
		return 1;
	}
	string piFolder = argv[1];
	cout << "Input Pi folder: " << piFolder << endl;

	double ***Pi;
	int windowSize = 0;
	double similarityThreshold = 0;

	ifstream seqfile;
	seqfile.open("Data/seqfile.txt");
	string newline = "";
	vector<int>* months = new vector<int>();
	vector<int>* numOfDocs = new vector<int>();
	vector<int>* monthFirstIdx = new vector<int>();
	vector<int>* monthLastIdx = new vector<int>();
	int curIdx = 0;

	while (seqfile >> newline) {
		const char * ptr = strchr(newline.c_str(), ':');
		int count = atoi(ptr + 1);
		ptr = "\0";
		int yearMonth = atoi(newline.c_str());
		months->push_back(yearMonth);
		numOfDocs->push_back(count);
		monthFirstIdx->push_back(curIdx);
		monthLastIdx->push_back(curIdx + count);
		curIdx += count;
	}
	seqfile.close();

	K = atoi(argv[2]);
	windowSize = atoi(argv[3]);
	similarityThreshold = atof(argv[4]);
	W = 32468;


	printf("Number of topics: %d\n", K);
	printf("Window Size: %d\n", windowSize);
	printf("Similarity Threshold: %f\n", similarityThreshold);

	// Dynamically allocate Pi
	Pi = new double**[months->size()];
	for (unsigned int m = 0; m < months->size(); ++m) {
		Pi[m] = new double*[W];
		for (unsigned int word = 0; word < W; word++) {
			Pi[m][word] = new double[K];
			for(unsigned int k = 0; k < K; k++) {
				Pi[m][word][k] = 0;
			}
		}
	}

	//Read Pi files in Memory
	for (int timeSlice = 0; timeSlice < (int)months->size(); timeSlice++) {
		string fileName = piFolder + "/topics_" + to_string(months->at(timeSlice)) + ".txt";
		cout << "Reading File: " << fileName << endl;
		ifstream pifile;
		pifile.open(fileName);
		int topic = 0;
		while (pifile >> newline) {
			std::istringstream ss(newline);
			std::string token;

			int wordId = 0;
			while (std::getline(ss, token, ',')) {
				Pi[timeSlice][wordId][topic] = stod(token);
				wordId++;
			}
			topic++;
		}
		pifile.close();

	} //All timeSlices finished
//	for (int timeSlice = 0; timeSlice < (int) months->size(); timeSlice++) {
//		for (int k = 0; k < K; k++) {
//			for (int w = 0; w < W; w++) {
//				cout << Pi[timeSlice][w][k] << ",";
//			}
//			cout << endl;
//		}
//	}

	// MAKE CHAINS
	Graph G;
	//K is unsigned -- is this a problem?
//    generateAllLinks(G, Pi, months->size(), K, windowSize, similarityThreshold);
	generateAllLinks(G, Pi, 10, K, windowSize, similarityThreshold);

	vector<int> component(num_vertices(G));
	int num = connected_components(G, &component[0]);

	vector<int>::size_type p;
	cout << "Total number of components: " << num << endl;
	for (p = 0; p != component.size(); ++p) {
		cout << "Vertex " << p << " is in component " << component[p] << endl;
	}

	return (0);

} // End of main
