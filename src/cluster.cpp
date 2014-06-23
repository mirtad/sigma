#include <cstdlib>

#include "cluster.h"

#include "sigma.h"

Cluster::Cluster(Contig* contig) {
	num_contigs_ = 1;
	contigs_ = new Contig*[num_contigs_];

	contigs_[0] = contig;

	length_ = contig->modified_length();

	sum_read_counts_ = new int[Sigma::num_samples];
	arrival_rates_ = new double[Sigma::num_samples];

	for (int sample_index = 0; sample_index < Sigma::num_samples; ++sample_index) {
		sum_read_counts_[sample_index] = contig->sum_read_counts()[sample_index];
		arrival_rates_[sample_index] = sum_read_counts_[sample_index] / (double) length_;
	}

	child1_ = NULL;
	child2_ = NULL;

	score_ = 0.0;
	model_score_ = 0.0;
	modeled_ = false;
	connected_ = false;

	contig->set_cluster(this);
}

Cluster::Cluster(Cluster* child1, Cluster* child2) {
	num_contigs_ = child1->num_contigs_ + child2->num_contigs_;
	contigs_ = new Contig*[num_contigs_];

	for (int contig_index = 0; contig_index < child1->num_contigs(); ++contig_index) {
		contigs_[contig_index] = child1->contigs_[contig_index];
	}

	int contig_index_offset = child1->num_contigs();

	for (int contig_index = 0; contig_index < child2->num_contigs(); ++contig_index) {
		contigs_[contig_index_offset + contig_index] = child2->contigs_[contig_index];
	}

	length_ = child1->length_ + child2->length_;

	sum_read_counts_ = new int[Sigma::num_samples];
	arrival_rates_ = new double[Sigma::num_samples];

	for (int sample_index = 0; sample_index < Sigma::num_samples; ++sample_index) {
		sum_read_counts_[sample_index] = child1->sum_read_counts_[sample_index] + child2->sum_read_counts_[sample_index];
		arrival_rates_[sample_index] = sum_read_counts_[sample_index] / (double) length_;
	}

	child1_ = child1;
	child2_ = child2;

	score_ = 0.0;
	model_score_ = 0.0;
	modeled_ = false;
	connected_ = false;

	for (int contig_index = 0; contig_index < num_contigs_; ++contig_index) {
		contigs_[contig_index]->set_cluster(this);
	}

	delete[] child1->contigs_;
	delete[] child2->contigs_;
}

Cluster::~Cluster() {
	delete[] sum_read_counts_;
	delete[] arrival_rates_;
}

Contig** Cluster::contigs() const { return contigs_; }
void Cluster::set_contigs(Contig** contigs) { contigs_ = contigs; }

int Cluster::num_contigs() const { return num_contigs_; }
int Cluster::length() const { return length_; }
int* Cluster::sum_read_counts() const { return sum_read_counts_; }
double* Cluster::arrival_rates() const { return arrival_rates_; }

Cluster* Cluster::child1() const { return child1_; }
Cluster* Cluster::child2() const { return child2_; }

double Cluster::score() const { return score_; }
double Cluster::model_score() const { return model_score_; }
bool Cluster::modeled() const { return modeled_; }
bool Cluster::connected() const { return connected_; }

void Cluster::set_score(double score) { score_ = score; }
void Cluster::set_model_score(double model_score) { model_score_ = model_score; }
void Cluster::set_modeled(bool modeled) { modeled_ = modeled; }
void Cluster::set_connected(bool connected) { connected_ = connected; }