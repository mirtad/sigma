#include <cstdlib>
#include <cstdio>

#include "contig_reader.h"

#include "sigma.h"

ContigReader::~ContigReader() {}


SOAPdenovoReader::SOAPdenovoReader() {}

void SOAPdenovoReader::read(const char* contigs_file, ContigMap* contigs) {
	char id[256];
	int length;

	FILE* contigs_fp = fopen(contigs_file, "r");

	if (contigs_fp != NULL) {
		while (!feof(contigs_fp)) {
			// >[ID] length [LENGTH] cvg_[COVERAGE]_tip_[TIP]\n
			if (fscanf(contigs_fp, ">%s %*s %d %*s\n", id, &length) == 2) {
				if (length >= Sigma::contig_len_thr) {
					contigs->insert(std::make_pair(id, new Contig(id, length)));
				}
			} else {
				fscanf(contigs_fp, "%*[^\n]\n");
			}
		}

		fclose(contigs_fp);
	} else {
		fprintf(stderr, "Error opening file: %s\n", contigs_file);
		exit(EXIT_FAILURE);
	}
}


VelvetReader::VelvetReader() {}

void VelvetReader::read(const char* contigs_file, ContigMap* contigs) {
	char id[256];
	int length;

	FILE* contigs_fp = fopen(contigs_file, "r");

	if (contigs_fp != NULL) {
		while (!feof(contigs_fp)) {
			// >NODE_[ID]_length_[LENGTH]_cov_[COVERAGE]\n
			if (fscanf(contigs_fp, ">%s\n", id) == 1 && sscanf(id, "%*[^_]_%*[^_]_%*[^_]_%d_%*s", &length) == 1) {
				if (length >= Sigma::contig_len_thr) {
					contigs->insert(std::make_pair(id, new Contig(id, length)));
				}
			} else {
				fscanf(contigs_fp, "%*[^\n]\n");
			}
		}

		fclose(contigs_fp);
	} else {
		fprintf(stderr, "Error opening file: %s\n", contigs_file);
		exit(EXIT_FAILURE);
	}
}