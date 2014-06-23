CC = g++
CFLAGS = -std=c++0x -O3 -Wall -Wconversion

all: sigma

sigma: sigma.o contig_reader.o mapping_reader.o edge_reader.o contig.o edge.o cluster.o cluster_graph.o probability_distribution.o
	$(CC) $(CFLAGS) sigma.o contig_reader.o mapping_reader.o edge_reader.o contig.o edge.o cluster.o cluster_graph.o probability_distribution.o -o sigma

sigma.o: sigma.cpp contig_reader.h mapping_reader.h edge_reader.h contig.h edge.h cluster.h cluster_graph.h probability_distribution.h
	$(CC) $(CFLAGS) -c sigma.cpp

contig_reader.o: contig_reader.cpp contig_reader.h sigma.h contig.h
	$(CC) $(CFLAGS) -c contig_reader.cpp

mapping_reader.o: mapping_reader.cpp mapping_reader.h sigma.h contig.h
	$(CC) $(CFLAGS) -c mapping_reader.cpp

edge_reader.o: edge_reader.cpp edge_reader.h sigma.h contig.h edge.h
	$(CC) $(CFLAGS) -c edge_reader.cpp

contig.o: contig.cpp contig.h sigma.h
	$(CC) $(CFLAGS) -c contig.cpp

edge.o: edge.cpp edge.h sigma.h contig.h
	$(CC) $(CFLAGS) -c edge.cpp

cluster.o: cluster.cpp cluster.h sigma.h contig.h
	$(CC) $(CFLAGS) -c cluster.cpp

cluster_graph.o: cluster_graph.cpp cluster_graph.h sigma.h contig.h edge.h cluster.h probability_distribution.h
	$(CC) $(CFLAGS) -c cluster_graph.cpp

probability_distribution.o: probability_distribution.cpp probability_distribution.h
	$(CC) $(CFLAGS) -c probability_distribution.cpp

clean:
	-rm *.o sigma