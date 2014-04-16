all: SCVB0/scvb0.cpp
	g++ -std=c++0x -fopenmp SCVB0/scvb0.cpp -o fastLDA

serial: SCVB0/scvb0.cpp
	g++ -o fastLDA SCVB0/scvb0.cpp -std=c++0x

clean:
	rm fastLDA
