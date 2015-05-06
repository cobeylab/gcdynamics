/**
 Header for class representing host (human individual).
 */

#ifndef __bcellmodel__Host__
#define __bcellmodel__Host__

#include <memory>
#include <vector>
#include "BCell.h"
#include "Database.hpp"
#include "NaiveCellGenerator.h"
#include "Simulation.h"

class SimParameters;
class DatabaseManager;
class Antigen;
class GC;

/**
 Class representing host.
*/
class Host
{
friend class GC;
friend class NaiveCellGenerator;
public:
	Host(Simulation & sim);
	void runInfection(Antigen & antigen, zppsim::rng_t & rng);
	uint32_t getNextCellId();
	
	double effectiveLogAffinity(int64_t epRow, int64_t epCol, BCell const & bCell, zppsim::rng_t & rng);
	double maxEffectiveLogAffinity(BCell const & bCell, zppsim::rng_t & rng);
	
	double meanMaxAffinityMemory(Antigen & antigen, zppsim::rng_t & rng);
	
	int64_t const epitopeRows;
	int64_t const epitopeCols;
	
	SequenceCache seqCache;
private:
	/*** PRIVATE FIELDS ***/
	
	Simulation * simPtr;
	SimParameters * p;
	
	uint32_t nextCellId;
	
	double time;
	DatabaseManager * dbmPtr;
	Antigen * antigenPtr;
	
	bool gcStarted;
	int64_t gcRound;
	
	NaiveCellGenerator naiveCellGen;
	std::vector<std::vector<double>> logAffinitiesMasking;
	
	int64_t nNaiveBCells;
	std::vector<std::unique_ptr<BCell>> naiveBCells;
	std::vector<std::unique_ptr<BCell>> memoryBCells;
	std::vector<std::unique_ptr<BCell>> plasmaCells;
	
	/*** PRIVATE FUNCTIONS ***/
	
	std::vector<std::unique_ptr<GC>> createGCsFromSeeds(
		Array<Array<String>> & seedStrs, zppsim::rng_t & rng
	);
	
	std::vector<std::unique_ptr<GC>> createGCsOldMethod(
		Antigen const & antigen, zppsim::rng_t & rng
	);
	std::vector<std::unique_ptr<GC>> createGCsNewMethod(
		Antigen const & antigen, zppsim::rng_t & rng
	);
	
	BCell * createNaiveCell(
		Antigen const & antigen, int64_t epRow, int64_t epCol,
		double logAffMin, double logAffMax, zppsim::rng_t & rng
	);
	
	double logAffinityQuantile(double quantile, int64_t epRow, int64_t epCol);
	double logAffinityCDF(double logAff, int64_t epRow, int64_t epCol);
	
	double getLogAffinityMasking(
		std::vector<std::vector<std::vector<double>>> const & logAffsPlasma,
		std::vector<double> logSumAffsPlasma,
		Antigen const & antigen, int64_t i, int64_t j, zppsim::rng_t & rng
	);
	
	double meanMaxAffinity(Antigen & antigen, std::vector<std::unique_ptr<BCell>> & cells, zppsim::rng_t & rng);
	
	void initializeNaiveCells(zppsim::rng_t & rng);
	
	void writeNaiveBCellsToDatabase();
	void writeNaiveBCellsToDatabase(std::vector<std::unique_ptr<BCell>> const & cells);
	
	void writeBCellsToDatabase(GC const & gc, std::string const & type, std::vector<std::unique_ptr<BCell>> const & bCells);
	void writeBCellsToDatabase(GC const & gc, std::string const & type, std::vector<BCell *> const & bCells);
	void writeBCellToDatabase(GC const & gc, std::string const & type, BCell const & bCell);
	
	void writePlasmaConcentrationToDatabase();
	
	void writeAffinityToDatabase(std::string const & type, BCell & cell, zppsim::rng_t & rng);
	void writeAffinityToDatabase(std::string const & type, std::vector<std::unique_ptr<BCell>> const & cells, zppsim::rng_t & rng);
	void writeAffinityStatsToDatabase(std::string const & type, std::vector<std::unique_ptr<BCell>> const & bCells, zppsim::rng_t & rng);
	
	void writeGCStatsToDatabase(int64_t nAlive);
	
	void writeGCSeedCandidatesToDatabase(std::vector<BCell *> const & seedBCells);
	void writeGCSeedsToDatabase(GC const & gc, std::vector<BCell *> const & initBCells);
	
	void stepTime(double delay);
	
	std::vector<BCell *> identifyGCSeedCandidates(
		Antigen const & antigen, zppsim::rng_t & rng
	);
};

#endif /* defined(__bcellmodel__Host__) */
