#include <cstdlib>

#include <algorithm>
#include <vector>

#include "contig.h"

#include "sigma.h"

Contig::Contig(std::string id, int length) : id_(id), length_(length) {
	if (Sigma::contig_window_len > 0) {
		num_windows_ = (length_ - 2 * Sigma::contig_edge_len) / Sigma::contig_window_len;

		const int remainder = length_ - num_windows_ * Sigma::contig_window_len;
		left_edge_ = remainder / 2;
		right_edge_ = length_ - 1 - (remainder - left_edge_);
	} else {
		num_windows_ = 1;

		left_edge_ = Sigma::contig_edge_len;
		right_edge_ = length_ - 1 - Sigma::contig_edge_len;
	}

	modified_length_ = right_edge_ - left_edge_ + 1;

	initReadCounts();

	cluster_ = NULL;
}

Contig::Contig(std::string id, int length, int left_edge, int right_edge, int num_windows) :
	id_(id), length_(length),
	left_edge_(left_edge), right_edge_(right_edge), num_windows_(num_windows) {
	modified_length_ = right_edge_ - left_edge_ + 1;

	initReadCounts();

	cluster_ = NULL;
}

void Contig::initReadCounts() {
	sum_read_counts_ = new int[Sigma::num_samples];

	for (int sample_index = 0; sample_index < Sigma::num_samples; ++sample_index) {
		sum_read_counts_[sample_index] = 0;
	}

	read_counts_ = new int*[Sigma::num_samples];

	for (int sample_index = 0; sample_index < Sigma::num_samples; ++sample_index) {
		read_counts_[sample_index] = new int[num_windows_];

		for (int window_index = 0; window_index < num_windows_; ++window_index) {
			read_counts_[sample_index][window_index] = 0;
		}
	}
}

Contig::~Contig() {
	delete[] sum_read_counts_;

	for (int sample_index = 0; sample_index < Sigma::num_samples; ++sample_index) {
		delete[] read_counts_[sample_index];
	}

	delete[] read_counts_;
}

std::string Contig::id() const { return id_; }
int Contig::length() const { return length_; }

int Contig::modified_length() const { return modified_length_; }
int Contig::left_edge() const { return left_edge_; }
int Contig::right_edge() const { return right_edge_; }
int Contig::num_windows() const { return num_windows_; }

int* Contig::sum_read_counts() const { return sum_read_counts_; }
int** Contig::read_counts() const { return read_counts_; }

Cluster* Contig::cluster() const { return cluster_; }

void Contig::set_cluster(Cluster* cluster) { cluster_ = cluster; }


void ContigIO::save_contigs(const ContigMap* contigs, const char* sigma_contigs_file_path) {
	FILE* sigma_contigs_fp = fopen(sigma_contigs_file_path, "w");

	if (sigma_contigs_fp != NULL) {
		fprintf(sigma_contigs_fp, "%d %d %d %d\n", Sigma::num_samples, Sigma::contig_len_thr, Sigma::contig_edge_len, Sigma::contig_window_len);

		for (auto it = contigs->begin(); it != contigs->end(); ++it) {
			Contig* contig = (*it).second;

			fprintf(sigma_contigs_fp, "%s\t%d\t%d\t%d\t%d\n",
					contig->id().c_str(), contig->length(), contig->left_edge(), contig->right_edge(), contig->num_windows());

			for (int sample_index = 0; sample_index < Sigma::num_samples; ++sample_index) {
				fprintf(sigma_contigs_fp, "%d\n", contig->sum_read_counts()[sample_index]);

				for (int window_index = 0; window_index < contig->num_windows(); ++window_index) {
					fprintf(sigma_contigs_fp, "%d ", contig->read_counts()[sample_index][window_index]);
				}

				fprintf(sigma_contigs_fp, "\n");
			}
		}

		fclose(sigma_contigs_fp);
	} else {
		fprintf(stderr, "Error opening file: %s\n", sigma_contigs_file_path);
		exit(EXIT_FAILURE);
	}
}

void ContigIO::load_contigs(const char* sigma_contigs_file_path, ContigMap* contigs) {
	FILE* sigma_contigs_fp = fopen(sigma_contigs_file_path, "r");

	if (sigma_contigs_fp != NULL) {
		fscanf(sigma_contigs_fp, "%d %d %d %d\n", &Sigma::num_samples, &Sigma::contig_len_thr, &Sigma::contig_edge_len, &Sigma::contig_window_len);

		while (!feof(sigma_contigs_fp)) {
			char id[256];
			int length, left_edge, right_edge, num_windows;

			fscanf(sigma_contigs_fp, "%s\t%d\t%d\t%d\t%d\n", id, &length, &left_edge, &right_edge, &num_windows);

			Contig* contig = new Contig(id, length, left_edge, right_edge, num_windows);

			for (int sample_index = 0; sample_index < Sigma::num_samples; ++sample_index) {
				fscanf(sigma_contigs_fp, "%d\n", &contig->sum_read_counts()[sample_index]);

				for (int window_index = 0; window_index < contig->num_windows(); ++window_index) {
					fscanf(sigma_contigs_fp, "%d ", &contig->read_counts()[sample_index][window_index]);
				}

				fscanf(sigma_contigs_fp, "\n");
			}

			contigs->insert(std::make_pair(id, contig));
		}

		fclose(sigma_contigs_fp);
	} else {
		fprintf(stderr, "Error opening file: %s\n", sigma_contigs_file_path);
		exit(EXIT_FAILURE);
	}
}


double compute_vmr(ContigMap* contigs) {
	std::vector<double> vmrs;

	for (auto it = contigs->begin(); it != contigs->end(); ++it) {
		Contig* contig = (*it).second;

		if (contig->length() < 10000) continue;

		for (int sample_index = 0; sample_index < Sigma::num_samples; ++sample_index) {
			int* read_counts = contig->read_counts()[sample_index];

			double mean = 0.0;

			for (int window_index = 0; window_index < contig->num_windows(); ++window_index) {
				mean += read_counts[window_index];
			}

			mean /= contig->num_windows();

			double variance = 0.0;

			for (int window_index = 0; window_index < contig->num_windows(); ++window_index) {
				variance += (read_counts[window_index] - mean) * (read_counts[window_index] - mean);
			}

			variance /= contig->num_windows();

			vmrs.push_back(variance / mean);
		}
	}

	std::sort(vmrs.begin(), vmrs.end());

	return vmrs[vmrs.size() / 2];
}