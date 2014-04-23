all: SCVB0/scvb0.cpp
	g++ -g -std=c++0x -fopenmp SCVB0/scvb0.cpp -o fastLDA

TopicChain: TopicChains/TopicChains.cpp
	g++ -g -std=c++0x -fopenmp TopicChains/TopicChains.cpp -o TopicChain
	
serial: SCVB0/scvb0.cpp
	g++ -g -std=c++0x SCVB0/scvb0.cpp -o fastLDA

clean:
	rm -f *.o fastLDA TopicChain doctopics.txt perplexity.txt
