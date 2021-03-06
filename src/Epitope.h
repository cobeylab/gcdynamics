#ifndef __bcellmodel__Epitope__
#define __bcellmodel__Epitope__

#include <memory>
#include <unordered_map>
#include "BCell.h"
#include "DatabaseManager.h"

class Antigen;
class PlasmaCell;

class Epitope
{
public:
	Epitope(
		SimParameters & params,
		Antigen * antigenPtr,
		uint32_t neighborSeed,
		uint32_t energySeed,
		zppsim::rng_t & rng, int64_t row, int64_t col, double energyMean,
		zppdb::Database * dbPtr,
		zppdb::Table<AntigenNeighborRow> * neighborTablePtr, zppdb::Table<AntigenEnergyRow> * energyTablePtr
	);
	void mutate(zppsim::rng_t & rng);
	
	double getEnergy(
		BCell const & cell
	);
	double getEnergy(
		BCell const & cell, uint16_t locus
	);
	
	void verifyNeighbors();
	void verifyEnergy();
	
	std::vector<uint8_t> const energySeedBytes;
	
	uint16_t const nLoci;
	int64_t const nInteractions;
	int64_t const alphabetSize;
	std::string const alphabet;
	
	double const pEnergyMutation;
	double const pNeighborMutation;
private:
	void writeNeighborsToDatabase();
	void writeNeighborToDatabase(
		uint32_t locus, uint32_t neighborIndex, uint32_t neighborLocus
	);
	void writeEnergyToDatabase(
		uint32_t locus, Sequence const & seq, double energy
	);
	
	int64_t row;
	int64_t col;
	double energyMean;
	
	std::vector<std::vector<uint32_t>> neighbors;
	std::vector<std::vector<Sequence>> activeNeighborSeqs;
//	std::vector<std::unordered_map<Sequence, double, HashSequence>> energyMaps;
	
	Antigen * antigenPtr;
	zppdb::Database * dbPtr;
	zppdb::Table<AntigenNeighborRow> * neighborTablePtr;
	zppdb::Table<AntigenEnergyRow> * energyTablePtr;
	std::normal_distribution<double> energyDist;
	
	std::unordered_map<uint32_t, double> energyCache;
};

#endif /* defined(__bcellmodel__Epitope__) */
