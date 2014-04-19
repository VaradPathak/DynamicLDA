all: SCVB0/scvb0.cpp
	g++ -g -std=c++0x -fopenmp SCVB0/scvb0.cpp -o fastLDA

serial: SCVB0/scvb0.cpp
	g++ -g -std=c++0x SCVB0/scvb0.cpp -o fastLDA

clean:
	rm -f *.o fastLDA doctopics.txt perplexity.txt
