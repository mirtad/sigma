#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "cluster_graph.h"

#include "sigma.h"

ClusterGraph::ClusterGraph(ContigMap* contigs, EdgeQueue* edges) {
	num_contigs_ = (int) contigs->size();
	num_windows_ = 0;

	for (auto it = contigs->begin(); it != contigs->end(); ++it) {
		Contig* contig = (*it).second;

		num_windows_ += contig->num_windows();

		roots_.insert(new Cluster(contig));
	}

	while (!edges->empty()) {
		Edge edge = edges->top();

		Cluster* cluster1 = edge.contig1()->cluster();
		Cluster* cluster2 = edge.contig2()->cluster();

		if (cluster1 != cluster2) {
			roots_.insert(new Cluster(cluster1, cluster2));

			roots_.erase(cluster1);
			roots_.erase(cluster2);
		}

		edges->pop();
	}

	updateClusters();
}

void ClusterGraph::updateClusters() {
	ClusterStack clusters;

	for (auto it = roots_.begin(); it != roots_.end(); ++it) {
		clusters.push(*it);
	}

	while (!clusters.empty()) {
		Cluster* cluster = clusters.top();
		clusters.pop();

		if (cluster->num_contigs() > 1) {
			cluster->child1()->set_contigs(cluster->contigs());
			cluster->child2()->set_contigs(cluster->contigs() + cluster->child1()->num_contigs());

			clusters.push(cluster->child1());
			clusters.push(cluster->child2());
		}
	}
}

ClusterGraph::~ClusterGraph() {
	ClusterStack clusters;

	for (auto it = roots_.begin(); it != roots_.end(); ++it) {
		Cluster* cluster = *it;

		for (int contig_index = 0; contig_index < cluster->num_contigs(); ++contig_index) {
			delete cluster->contigs()[contig_index];
		}

		delete[] cluster->contigs();

		clusters.push(cluster);
	}

	while (!clusters.empty()) {
		Cluster* cluster = clusters.top();
		clusters.pop();

		if (cluster->num_contigs() > 1) {
			clusters.push(cluster->child1());
			clusters.push(cluster->child2());
		}

		delete cluster;
	}
}

ClusterSet* ClusterGraph::roots() { return &roots_; }

void ClusterGraph::computeScores(const ProbabilityDistribution* prob_dist) {
	ClusterStack clusters;

	for (auto it = roots_.begin(); it != roots_.end(); ++it) {
		clusters.push(*it);
	}

	while (!clusters.empty()) {
		Cluster* cluster = clusters.top();
		clusters.pop();

		if (cluster->num_contigs() > 1) {
			clusters.push(cluster->child1());
			clusters.push(cluster->child2());
		}

		computeClusterScore(cluster, prob_dist);
	}
}

void ClusterGraph::computeModels() {
	ClusterStack clusters;

	for (auto it = roots_.begin(); it != roots_.end(); ++it) {
		clusters.push(*it);
	}

	while (!clusters.empty()) {
		Cluster* cluster = clusters.top();

		if (cluster->num_contigs() == 1 || (cluster->child1()->modeled() && cluster->child2()->modeled())) {
			clusters.pop();
			computeClusterModel(cluster);
		} else {
			clusters.push(cluster->child1());
			clusters.push(cluster->child2());
		}
	}

	for (auto it = roots_.begin(); it != roots_.end(); ++it) {
		clusters.push(*it);
	}

	while (!clusters.empty()) {
		Cluster* cluster = clusters.top();
		clusters.pop();

		if (cluster->connected()) {
			for (int contig_index = 0; contig_index < cluster->num_contigs(); ++contig_index) {
				cluster->contigs()[contig_index]->set_cluster(cluster);
			}
		} else {
			clusters.push(cluster->child1());
			clusters.push(cluster->child2());
		}
	}
}

void ClusterGraph::computeClusterScore(Cluster* cluster, const ProbabilityDistribution* prob_dist) {
	double score = 0;

	for (int sample_index = 0; sample_index < Sigma::num_samples; ++sample_index) {
		double mean_read_count = 0.0;

		if (Sigma::contig_window_len > 0) {
			mean_read_count = cluster->arrival_rates()[sample_index] * Sigma::contig_window_len;
		}

		for (int contig_index = 0; contig_index < cluster->num_contigs(); ++contig_index) {
			Contig* contig = cluster->contigs()[contig_index];

			if (Sigma::contig_window_len == 0) {
				mean_read_count = cluster->arrival_rates()[sample_index] * contig->modified_length();
				score += prob_dist->logpf(mean_read_count, contig->sum_read_counts()[sample_index]);
			} else {
				for (int window_index = 0; window_index < contig->num_windows(); ++window_index) {
					score += prob_dist->logpf(mean_read_count, contig->read_counts()[sample_index][window_index]);
				}
			}
		}
	}

	score -= 0.5 * Sigma::num_samples * log(num_windows_);

	cluster->set_score(score);
}

void ClusterGraph::computeClusterModel(Cluster* cluster) {
	if (cluster->num_contigs() == 1) {
		cluster->set_model_score(cluster->score());
		cluster->set_connected(true);
	} else {
		double connected_score = cluster->score();
		double disconnected_score = cluster->child1()->model_score() + cluster->child2()->model_score();

		if (connected_score >= disconnected_score) {
			cluster->set_model_score(connected_score);
			cluster->set_connected(true);
		} else {
			cluster->set_model_score(disconnected_score);
			cluster->set_connected(false);
		}
	}

	cluster->set_modeled(true);
}

void ClusterGraph::saveClusters(const char* clusters_file_path) {
	FILE* clusters_fp = fopen(clusters_file_path, "w");

	if (clusters_fp != NULL) {
		int cluster_id = 1;

		ClusterStack clusters;

		for (auto it = roots_.begin(); it != roots_.end(); ++it) {
			clusters.push(*it);
		}

		while (!clusters.empty()) {
			Cluster* cluster = clusters.top();
			clusters.pop();

			if (cluster->connected()) {
				for (int contig_index = 0; contig_index < cluster->num_contigs(); ++contig_index) {
					Contig* contig = cluster->contigs()[contig_index];

					fprintf(clusters_fp, "%s\t%d\t%d\t%f\n",
							contig->id().c_str(), cluster_id, contig->sum_read_counts()[0], cluster->arrival_rates()[0]);
				}

				cluster_id++;
			} else {
				clusters.push(cluster->child1());
				clusters.push(cluster->child2());
			}
		}

		fclose(clusters_fp);
	} else {
		fprintf(stderr, "Error opening file: %s\n", clusters_file_path);
		exit(EXIT_FAILURE);
	}
}