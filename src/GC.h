/**
 Header for class representing germinal center (location of B-cell affinity maturation).
*/

#ifndef __bcellmodel__GC__
#define __bcellmodel__GC__

#include <memory>
#include <vector>
#include <random>

#include "Host.h"
#include "BCell.h"
#include "PlasmaCell.h"

class DatabaseManager;

/**
 Class representing germinal center.
*/
class GC
{
friend class Host;

public:
	GC(Host & host, std::vector<BCell *> const & initCells, int64_t id, zppsim::rng_t & rng);
	void runRound(
		zppsim::rng_t & rng, Host & host, Antigen & antigen,
		std::vector<BCell *> & gcfdCellsForMemoryExport,
		std::vector<BCell *> & gcfdCellsForPlasmaExport,
		double pMemExport
	);
	int64_t const id;
	bool isAlive();
	
	void writeGCFDCellsToDatabase(Antigen const & antigen, int64_t round);
	void writeEffectiveAffinityToDatabase(
		Host & host,
		Antigen const & antigen,
		int64_t round,
		zppsim::rng_t & rng
	);
private:
	Host * hostPtr;
	DatabaseManager * dbmPtr;
	SimParameters * p;
	
	bool alive;
	std::vector<std::unique_ptr<BCell>> cells;
	
	std::vector<BCell *> sampleNewCells(int64_t newCells, int64_t maxOffspring, Host & host, zppsim::rng_t & rng);
	std::vector<BCell *> sampleNewCellsNullModel(int64_t newCells, int64_t maxOffspring, Host & host, zppsim::rng_t & rng);
};

#endif /* defined(__bcellmodel__GC__) */
