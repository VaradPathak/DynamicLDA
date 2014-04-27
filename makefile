all: SCVB0/scvb0.cpp
	g++ -g -std=c++0x -fopenmp SCVB0/scvb0.cpp -o fastLDA

TopicChain: TopicChains/TopicChains.cpp
	g++ -g -std=c++0x -fopenmp TopicChains/TopicChains.cpp -o TopicChain
	
TopicChain_GetData: TopicChains/TopicChains_GetData.cpp
	g++ -g -std=c++0x -fopenmp TopicChains/TopicChains_GetData.cpp -o TopicChain_GetData
	
serial: SCVB0/scvb0.cpp
	g++ -g -std=c++0x SCVB0/scvb0.cpp -o fastLDA

clean:
	rm -f *.o fastLDA TopicChain TopicChain_GetData
