#ifndef CLUSTER_GRAPH_H_
#define CLUSTER_GRAPH_H_

#include "contig.h"
#include "edge.h"
#include "cluster.h"
#include "probability_distribution.h"

/**
 * @brief A class for representing a graph of hierarchical clustering trees.
 *
 * Represents a greph of hierarchical clustering trees builit along the scaffold
 * edges based on contig arrival rate information.
 */
class ClusterGraph {
public:
	/**
	 * @brief Constructs a graph of hierarchical clustering trees.
	 *
	 * @param contigs	contigs
	 * @param edges		edges
	 */
	ClusterGraph(ContigMap* contigs, EdgeQueue* edges);

	~ClusterGraph(); /**< Default destructor. */

	/**
	 * @brief Getter for roots of hierarchical clustering trees.
	 *
	 * @return roots of hierarchical clustering trees
	 */
	ClusterSet* roots();

	/**
	 * @brief Computes scores for all clusters based on given probability distribution.
	 *
	 * @param prob_dist		probability distribution
	 */
	void computeScores(const ProbabilityDistribution* prob_dist);

	/**
	 * @brief Computes models for all clustering trees which maximize BIC.
	 */
	void computeModels();

	/**
	 * @brief Saves final clusters to a file.
	 *
	 * @param clusters_file_path	path to file for saving final clusters
	 */
	void saveClusters(const char* clusters_file_path);

private:
	/**
	 * @brief Updates contig array pointers of all clusters.
	 */
	void updateClusters();

	/**
	 * @brief Computes score for the cluster based on given probability distribution.
	 *
	 * @param cluster		cluster
	 * @param prob_dist		probability distribution
	 */
	void computeClusterScore(Cluster* cluster, const ProbabilityDistribution* prob_dist);

	/**
	 * @brief Computes model for the cluster which maximizes BIC.
	 *
	 * @param cluster	cluster
	 */
	void computeClusterModel(Cluster* cluster);

	int num_contigs_; /**< Number of contigs. */
	int num_windows_; /**< Number of windows. */
	ClusterSet roots_; /**< Roots of hierarchical clustering trees. */
};

#endif // CLUSTER_GRAPH_H_