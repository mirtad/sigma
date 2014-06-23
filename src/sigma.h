#ifndef SIGMA_H_
#define SIGMA_H_

#include <string>
#include <vector>
#include <unordered_map>

/** A map for storing key-value configuration parameters. */
typedef std::unordered_map<std::string, std::string> ParamsMap;


/**
 * @brief Main configuration class.
 *
 * This class is used for reading the configuration file and storing
 * configuration parameters.
 */
class Sigma {
public:
	/**
	 * @brief Reads the configuration file.
	 *
	 * @param config_file	path to configuration file
	 */
	static void readConfigFile(char* config_file);

	static std::string contigs_file_type; /**< Type of contigs file. */

	static std::string contigs_file; /**< Path to contigs file. */
	static std::vector<std::string> mapping_files; /**< Paths to mapping files. */
	static std::vector<std::string> edges_files; /**< Paths to edges files. */

	static std::string sigma_contigs_file; /**< Path to Sigma contigs file. */

	static std::string output_dir; /**< Path to output directory. */
	static std::vector<std::string> skipped_edges_files; /**< Paths to skipped edges files. */
	static std::vector<std::string> filtered_edges_files; /**< Paths to filtered edges files. */
	static std::string clusters_file;  /**< Path to clusters file. */

	static int num_samples; /**< Number of samples. */
	
	static int contig_len_thr; /**< Threshold on contig length. */
	static int contig_edge_len; /**< Contig edge length. */
	static int contig_window_len; /**< Contig window length. */

	static std::string pdist_type; /**< Type of read count probability distribution. */

	static double vmr; /**< Variance to mean ratio for negative binomial distribution. */

private:
	/**
	 * @brief Configures all parameters from the given map.
	 *
	 * @param params	map containing key-value configuration parameters
	 */
	static void configure(ParamsMap* params);

	/**
	 * @brief Returns int value for given key.
	 * 
	 * @param params	map containing key-value configuration parameters
	 * @param key		key
	 * @return int value for given key, or -1 if the key does not exist
	 */
	static int getIntValue(ParamsMap* params, std::string key);

	/**
	 * @brief Returns double value for given key.
	 *
	 * @param params	map containing key-value configuration parameters
	 * @param key		key
	 * @return double value for given key, or -1.0 if the key does not exist
	 */
	static double getDoubleValue(ParamsMap* params, std::string key);

	/**
	 * @brief Returns string value for given key.
	 *
	 * @param params	map containing key-value configuration parameters
	 * @param key		key
	 * @return string value for given key, or "-" if the key does not exist
	 */
	static std::string getStringValue(ParamsMap* params, std::string key);

	/**
	 * @brief Returns vector of string values for given key.
	 *
	 * @param params	map containing key-value configuration parameters
	 * @param key		key
	 * @return vector of string values for given key, or an empty vector if the key does not exist
	 */
	static std::vector<std::string> getVectorValue(ParamsMap* params, std::string key);
};

#endif // SIGMA_H_