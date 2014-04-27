all: SCVB0/scvb0.cpp
	g++ -g -std=c++0x -fopenmp SCVB0/scvb0.cpp -o fastLDA

GenerateChains: TopicChains/GenerateChains.cpp
	g++ -g -std=c++0x -fopenmp TopicChains/GenerateChains.cpp -o GenerateChains
	
GetData: TopicChains/GetData.cpp
	g++ -g -std=c++0x -fopenmp TopicChains/GetData.cpp -o GetData
	
serial: SCVB0/scvb0.cpp
	g++ -g -std=c++0x SCVB0/scvb0.cpp -o fastLDA

clean:
	rm -f *.o fastLDA GetData GenerateChains
