#ifndef CONTIG_READER_H_
#define CONTIG_READER_H_

#include "contig.h"

/**
 * @brief An interface for contigs file readers.
 *
 * Provides an interface for reading contig information from the output of an
 * arbitrary assembly tool.
 */
class ContigReader {
public:
	virtual ~ContigReader(); /**< A virtual destructor. */

	/**
	 * @brief Reads contig information from given contigs file.
	 *
	 * @param contigs_file	path to contigs file
	 * @param contigs		map with contig information
	 */
	virtual void read(const char* contigs_file, ContigMap* contigs) = 0;
};


/**
 * @brief SOAPdenovo contigs file reader.
 *
 * Enables reading contig information from the output of
 * <a href="http://soap.genomics.org.cn/soapdenovo.html">SOAPdenovo</a>
 * assembler.
 */
class SOAPdenovoReader : public ContigReader {
public:
	SOAPdenovoReader(); /**< An empty constructor. */

	/**
	 * @brief Reads contig information from SOAPdenovo contigs file.
	 *
	 * @copydetails ContigReader::read(const char*, ContigMap*)
	 */
	void read(const char* contigs_file, ContigMap* contigs);
};


/**
 * @brief Velvet contigs file reader.
 *
 * Enables reading contig information from the output of
 * <a href="https://www.ebi.ac.uk/~zerbino/velvet/">Velvet</a>
 * assembler.
 */
class VelvetReader : public ContigReader {
public:
	VelvetReader(); /**< An empty constructor. */

	/**
	 * @brief Reads contig information from Velvet contigs file.
	 * 
	 * @copydetails ContigReader::read(const char*, ContigMap*)
	 */
	void read(const char* contigs_file, ContigMap* contigs);
};

#endif // CONTIG_READER_H_