#ifndef PROBABILITY_DISTRIBUTION_H_
#define PROBABILITY_DISTRIBUTION_H_

/**
 * @brief An interface for probability distributions.
 *
 * This class provides an interface for probability distributions for
 * modeling contig read count (and possibly other) statistics.
 */
class ProbabilityDistribution {
public:
	virtual ~ProbabilityDistribution(); /**< A virtual destructor. */

	/**
	 * @brief Computes log of pmf/pdf for given mean and value.
	 *
	 * Computes the natural logarithm of probability mass/density function
	 * for given mean and value.
	 *
	 * @param mean		mean of the distribution
	 * @param value		value for which log of pmf/pdf is computed
	 * @return log of pmf/pdf for given mean and value
	 */
	virtual double logpf(double mean, double value) const = 0;
};


/**
 * @brief An implementation of Poisson probability distribution.
 *
 * Implements <a href="http://en.wikipedia.org/wiki/Poisson_distribution">
 * Poisson probability distribution</a>.
 */
class PoissonDistribution : public ProbabilityDistribution {
public:
	PoissonDistribution(); /**< An empty constructor. */

	/**
	 * @brief Computes log of pmf for given mean and value.
	 *
	 * @param mean		mean of the distribution
	 * @param value		value for which log of pmf is computed
	 * @return log of pmf for given mean and value
	 */
	double logpf(double mean, double value) const;
};


/**
 * @brief An implementation of negative binomial probability distribution.
 *
 * Implements <a href="http://en.wikipedia.org/wiki/Negative_binomial_distribution">
 * negative binomial probability distribution</a>.
 */
class NegativeBinomialDistribution : public ProbabilityDistribution {
public:
	/**
	 * @brief Constructs a negative binomial distribution for given vmr value.
	 *
	 * @param vmr	variance-to-mean ratio
	 */
	NegativeBinomialDistribution(double vmr);

	/**
	 * @brief Computes log of pmf for given mean and value.
	 *
	 * @param mean		mean of the distribution
	 * @param value		value for which log of pmf is computed
	 * @return log of pmf for given mean and value
	 */
	double logpf(double mean, double value) const;

private:
	const double log_p; /**< Log probability of success. */
	const double log_1mp; /**< Log probability of failure. */
	const double xo1mx; /**< Multiplier for computing mean number of failures from mean number of successes. */
};

#endif // PROBABILITY_DISTRIBUTION_H_