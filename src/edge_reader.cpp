#include <cstdlib>
#include <cstdio>

#include "edge_reader.h"

#include "sigma.h"

EdgeReader::~EdgeReader() {}


OperaBundleReader::OperaBundleReader() {}

void OperaBundleReader::read(const char* edges_file, const ContigMap* contigs, EdgeSet* edges, const char* skipped_edges_file) {
	char id1[256], id2[256], line[1024];

	FILE* edges_fp = fopen(edges_file, "r");

	if (edges_fp != NULL) {
		FILE* skipped_edges_fp = fopen(skipped_edges_file, "w");

		if (skipped_edges_fp == NULL) {
			fprintf(stderr, "Error opening file: %s\n", skipped_edges_file);
			exit(EXIT_FAILURE);
		}

		while (!feof(edges_fp)) {
			fscanf(edges_fp, "%[^\n]\n", line);

			// [ID1]\t[ORIENTATION1]\t[ID2]\t[ORIENTATION2]\t[DISTANCE]\t[STDEV]\t[SIZE]\n
			if (sscanf(line, "%s\t%*c\t%s\t%*c\t%*[^\n]", id1, id2) == 2) {
				auto it1 = contigs->find(id1);
				auto it2 = contigs->find(id2);

				if (it1 != contigs->end() && it2 != contigs->end()) {
					Contig* contig1 = (*it1).second;
					Contig* contig2 = (*it2).second;

					if (contig1 != contig2) {
						edges->insert(Edge(contig1, contig2));
					}
				} else {
					fprintf(skipped_edges_fp, "%s\n", line);
				}
			}
		}

		fclose(skipped_edges_fp);

		fclose(edges_fp);
	} else {
		fprintf(stderr, "Error opening file: %s\n", edges_file);
		exit(EXIT_FAILURE);
	}
}

void OperaBundleReader::filter(const char* edges_file, const ContigMap* contigs, const char* filtered_edges_file) {
	char id1[256], id2[256], line[1024];

	FILE* edges_fp = fopen(edges_file, "r");

	if (edges_fp != NULL) {
		FILE* filtered_edges_fp = fopen(filtered_edges_file, "w");

		if (filtered_edges_fp == NULL) {
			fprintf(stderr, "Error opening file: %s\n", filtered_edges_file);
			exit(EXIT_FAILURE);
		}

		while (!feof(edges_fp)) {
			fscanf(edges_fp, "%[^\n]\n", line);

			// [ID1]\t[ORIENTATION1]\t[ID2]\t[ORIENTATION2]\t[DISTANCE]\t[STDEV]\t[SIZE]\n
			if (sscanf(line, "%s\t%*c\t%s\t%*c\t%*[^\n]", id1, id2) == 2) {
				auto it1 = contigs->find(id1);
				auto it2 = contigs->find(id2);

				if (it1 != contigs->end() && it2 != contigs->end() &&
						((*it1).second->cluster() == (*it2).second->cluster())) {
					fprintf(filtered_edges_fp, "%s\n", line);
				}
			}
		}

		fclose(filtered_edges_fp);

		fclose(edges_fp);
	} else {
		fprintf(stderr, "Error opening file: %s\n", edges_file);
		exit(EXIT_FAILURE);
	}
}