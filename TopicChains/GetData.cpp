/*
 * TopicChains.cpp
 *
 *  Created on: Apr 22, 2014
 *      Author: vspathak
 */

#include <math.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <time.h>
#include <random>
#include <cstring>
#include <map>
#include <limits>
#include <boost/config.hpp>
#include <algorithm>
#include <utility>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

using namespace std;
using namespace boost;

typedef adjacency_list <vecS, vecS, undirectedS> Graph;

// Initialize number of documents, topics and words in vocabulary
unsigned int W, D, K;

int main(int argc, char* argv[]) {
    if (argc < 4) {
        printf("Usage: ./fastLDA inputfile num_iterations num_topics\n");
        return 1;
    }

    // Initlialize expected topic counts per document
    double **nTheta;
    // Dynamically
    double **nPi;
    double *N_z;
    // Initialize estimates from each minibatch
    // Initialize step sizes
    double rhoTheta = 0;
    double rhoPhi = 0;
    double ***Pi;
    double **theta;
    double **perplexities;
    // Initlalize dirichlet prior parameters
    double alpha, eta;
    double M; // Number of documents in each minibatch
    int Cj = 0;
    unsigned int i, j, k, w, MAXITER;
    int batch_idx = 0;
    int C = 0;
    int iter = 0;
    int NNZ;
	double perplexityval, innerval;
	ofstream pfile;
	pfile.open("perplexity.txt");

    M = 100; //343 works for KOS and only for KOS
    eta = 0.01; // was 0.01
    alpha = 0.1;

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
                monthLastIdx->push_back(curIdx+count);
                curIdx += count;
    }
    seqfile.close();

    //if user also specified a minibatch size
    if (argc > 4) {
        M = atof(argv[4]);
    }

    MAXITER = atoi(argv[2]);
    K = atoi(argv[3]);

    printf("Input file: %s\n", argv[1]);
    printf("Number of iterations: %d\n", MAXITER);
    printf("Number of topics: %d\n", K);
    printf("Minibatch size: %f\n", M);
    printf("alpha:  %f\n", alpha);
    printf("eta:  %f\n", eta);

    // Read the file and store it in DATA
    FILE* fptr;
    unsigned int docnum, wnum;
    unsigned char countnum;

    fptr = fopen(argv[1], "rt");

    fscanf(fptr, "%d\n", &D);
    fscanf(fptr, "%d\n", &W);
    fscanf(fptr, "%d\n", &NNZ);

    printf("Number of documents: %d\n", D);
    printf("Vocabulary size: %d\n", W);

    // Dynamically allocate phi
    Pi = new double**[months->size()];
    for (unsigned int m = 0; m < months->size(); ++m) {
        Pi[m] = new double*[W];
        for (unsigned int word = 0; word < W; word++) {
            Pi[m][word] = new double[K];
        }
    }
//#pragma omp parallel for


    printf("allocated phi\n");

    // Dynamically allocate theta

    theta = new double*[D];
//#pragma omp parallel for
    for (i = 0; i < D; i++) {
        theta[i] = new double[K];
    }

    printf("allocated theta\n");

    vector<vector<int> > corpus;
    vector<int> corpus_size(D, 0);
    corpus.resize(D);
    vector<vector<int> > corpus_expanded;
    corpus_expanded.resize(D);

    while (!feof(fptr)) {
        fscanf(fptr, "%d %d %hhu\n", &docnum, &wnum, &countnum);

        corpus[docnum - 1].push_back(wnum - 1);
        corpus[docnum - 1].push_back(countnum);

        corpus_size[docnum - 1] += countnum;

        for (i = 0; i < countnum; i++) {
            corpus_expanded[docnum - 1].push_back(wnum - 1);
        }
    }
    fclose(fptr);


    // Initialize phi_est and all other arrays
    nPi = new double*[W];

    for (i = 0; i < W; i++) {
        nPi[i] = new double[K];
    }

    // Initialize n_z and n_z_est and other arrays
    N_z = new double[K];
    for (k = 0; k < K; k++) {
        N_z[k] = 0;
    }

    nTheta = new double*[D];
    for (i = 0; i < D; i++) {
        nTheta[i] = new double[K];
    }

    for (i = 0; i < D; i++) {
        for (k = 0; k < K; k++) {
            nTheta[i][k] = rand() % 10;
        }
    }

    perplexities = new double*[months->size()];
    for (i = 0; i < months->size(); i++) {
        perplexities[i] = new double[K];
        for (unsigned int a = 0; a < K; ++a) {
            perplexities[i][a] = 0;
        }
    }

    int*** topwords;
    topwords = new int**[months->size()];

    //Generate Numbers according to Gaussian Distribution

    for (int timeSlice = 0; timeSlice < months->size(); timeSlice++) {
        cout << (*months)[timeSlice] << " " << (*numOfDocs)[timeSlice] << endl;

        //if parallelizing this, make sure to avoid race condition (most likely use reduction)
        for (k = 0; k < K; k++) {
        	N_z[k] = 0;
            for (w = 0; w < W; w++) {
                N_z[k] += nPi[w][k];
            }
        }

        // Find the total number of word in the document
        int monthFirstDoc = monthFirstIdx->at(timeSlice);
        int monthLastDoc = monthLastIdx->at(timeSlice);

        int monthD = monthLastDoc - monthFirstDoc;

        C = 0;

        for (int var = monthFirstDoc; var < monthLastDoc; var++) {
            C += corpus_size[var];
        }

        printf("Number of words in corpus: %d\n", C);

        int firstdoc = 0;
        int lastdoc = 0;
        int DM = monthD / M;

        for (iter = 0; iter < (int)MAXITER; iter++) {
            // Decide rho_phi and rho_theta
            rhoPhi = 10 / pow((1000 + iter), 0.9);
            rhoTheta = 1 / pow((10 + iter), 0.9);

#pragma omp parallel private(batch_idx,j,k,i,w,firstdoc,lastdoc)
            {
                double *gamma = new double[K];
                double *nzHat = new double[K];
                double **nPhiHat = new double *[W];
                for (k = 0; k < K; k++) {
                    gamma[k] = 0;
                    nzHat[k] = 0;
                }
                for (i = 0; i < W; i++) {
                    nPhiHat[i] = new double[K];
                    for (k = 0; k < K; k++) {
                        nPhiHat[i][k] = 0;
                    }
                }

#pragma omp for
                for (batch_idx = 0; batch_idx < DM+1; batch_idx++) {

                    // Decide the document indices which go in each minibatch
                    firstdoc = monthFirstDoc + (batch_idx * M);
                    lastdoc = monthFirstDoc + ((batch_idx + 1) * M);

                    if (batch_idx == DM) {
                        lastdoc = monthLastDoc;
                    }
                    for (j = (unsigned)firstdoc; j < (unsigned)lastdoc; j++) {

                        // First perform the burn-in passes
                        // Iteration of burn in passes

                        // Store size of corpus in Cj
                        Cj = corpus_size[j];

                        for (i = 0; i < (corpus[j].size() / 2); i++) {// indexing is very different here!

                            int w_aj = corpus[j][2 * i];
                            int m_aj = corpus[j][(2 * i) + 1];
                            // Update gamma_ij and N_theta
                            double normSum = 0;

                            for (k = 0; k < K; k++) {
                                gamma[k] = (nPi[w_aj][k] + eta) * (nTheta[j][k] + alpha) / (N_z[k] + (eta * W));
                                normSum += gamma[k];
                            }

                            for (k = 0; k < K; k++) {
                                gamma[k] = gamma[k] / normSum;
                            }

                            for (k = 0; k < K; k++) {

                                nTheta[j][k] = (pow((1 - rhoTheta), m_aj) * nTheta[j][k])
                                        + ((1 - pow((1 - rhoTheta), m_aj)) * Cj * gamma[k]);
                            }

                        }

                        // Iteration of the main loop
                        for (i = 0; i < (corpus[j].size() / 2); i++) { // indexing is very different here!

                            int w_aj = corpus[j][2 * i];
                            int m_aj = corpus[j][(2 * i) + 1];
                            double normSum = 0;
                            for (k = 0; k < K; k++) {
                                gamma[k] = (nPi[w_aj][k] + eta) * (nTheta[j][k] + alpha) / (N_z[k] + (eta * W));
                                normSum += gamma[k];
                            }

                            for (k = 0; k < K; k++) {
                                gamma[k] = gamma[k] / normSum;
                            }

                            // Update N_theta estimates
                            for (k = 0; k < K; k++) {
                                nTheta[j][k] = (pow((1 - rhoTheta), m_aj) * nTheta[j][k])
                                        + ((1 - pow((1 - rhoTheta), m_aj)) * Cj * gamma[k]);

                                nPhiHat[w_aj][k] = nPhiHat[w_aj][k] + (C * gamma[k] / M);

                                nzHat[k] = nzHat[k] + (C * gamma[k] / M);
                            }
                        }

                    } // End of j

                    // Update the estimates matrix
                    for (k = 0; k < K; k++) {
                        for (w = 0; w < W; w++) {
                            nPi[w][k] = (1 - rhoPhi) * nPi[w][k] + rhoPhi * nPhiHat[w][k];
                        }
#pragma omp atomic
                        N_z[k] *= (1 - rhoPhi);
#pragma omp atomic
                        N_z[k] += rhoPhi * nzHat[k];
                    }

                } // End of batch_idx

                // Compute phi
#pragma omp for
                for (k = 0; k < K; k++) {
                    double normSum = 0;
                    for (w = 0; w < W; w++) {
                        nPi[w][k] += eta;
                        normSum += nPi[w][k];
                    }
//                    cout << normSum << endl;
                    for (w = 0; w < W; w++) {
                        Pi[timeSlice][w][k] = (double) nPi[w][k] / normSum;
                    }
                }

                // Compute theta
#pragma omp for
                for (int var = monthFirstDoc; var < monthLastDoc; var++) {
                    double normSum = 0;
                    for (k = 0; k < K; k++) {
                        nTheta[var][k] += alpha;
                        normSum += nTheta[var][k];
                    }
                    for (k = 0; k < K; k++) {
                        theta[i][k] = (double) nTheta[i][k] / normSum;
                    }
                }

                delete[] gamma;
                delete[] nzHat;

                for (i = 0; i < W; i++) {
                    delete[] nPhiHat[i];
                }

                delete[] nPhiHat;

            }

            // Calculate the perplexity here
			perplexityval = 0;
#pragma omp parallel for private(j,i,k) reduction(+:innerval) reduction(+:perplexityval)
			for (j = monthFirstDoc; j < monthLastDoc; j++) {
				for (i = 0; i < corpus_expanded[j].size(); i++) {
					innerval = 0;
					for (k = 0; k < K; k++) {
						innerval += (theta[j][k] * Pi[timeSlice][corpus_expanded[j][i]][k]);
					}
					perplexityval += (log(innerval) / log(2));
				}
			}
			printf("%d,%f\n", iter, pow(2, -perplexityval / C));
			perplexities[timeSlice][iter] = pow(2, -perplexityval / C);

			pfile << iter + 1 << "," << perplexities[timeSlice][iter] << endl;
			pfile.flush();

        } // End of iter

        //compute the top 100 words for each topic

        double** maxval;
        topwords[timeSlice] = new int*[K];
        maxval = new double*[K];
        for (k = 0; k < K; k++) {
            topwords[timeSlice][k] = new int[100];
            maxval[k] = new double[100];
        }
        for (k = 0; k < K; k++) {
            double oldMax = std::numeric_limits<double>::max();
            for (i = 0; i < 100; i++) {
                double max = -1;
                int max_idx = -1;
                for (w = 0; w < W; w++) {
                    if (oldMax > Pi[timeSlice][w][k] && Pi[timeSlice][w][k] > max) {
                        max = Pi[timeSlice][w][k];
                        max_idx = w;
                    }
                }
                oldMax = Pi[timeSlice][max_idx][k];
                topwords[timeSlice][k][i] = max_idx;
                maxval[k][i] = max;
            }
        }

        ofstream pifile;
		pifile.open("Pi/topics_" + to_string(months->at(timeSlice)) + ".txt");
		for (k = 0; k < K; k++) {
			for (w = 0; w < W-1; w++) {
				pifile << Pi[timeSlice][w][k] << ",";
			}
			pifile << Pi[timeSlice][W-1][k];
			pifile << endl;
		}
		pifile.close();
    }//All timeSlices finished

    string *dict;
    dict = new string[W];
//    char word;
    //retrieve the words from the file
    w = 0;
    string line;
    ifstream vocabFile(argv[5]);
    if (vocabFile.is_open()) {
        while (getline(vocabFile, line)) {
            dict[w] = line;
            w++;
        }
        vocabFile.close();
    }

//    write topics file
    for (int timeSlice = 0; timeSlice < (int) months->size(); timeSlice++) {
        ofstream tfile;
        tfile.open("output/topics_" + to_string(months->at(timeSlice)) + ".txt");
        for (k = 0; k < K; k++) {
            for (w = 0; w < 100; w++) {
                tfile << topwords[timeSlice][k][w];// << ":" << maxval[k][w] << ",";

            }
            tfile << endl;

            for (w = 0; w < 100; w++) {
                tfile << dict[topwords[timeSlice][k][w]] << ",";
            }

            tfile << endl;
        }
        tfile.close();
    }

    return (0);

} // End of main
