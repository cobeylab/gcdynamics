#ifndef __bcellmodel__Antigen__
#define __bcellmodel__Antigen__

#include <memory>
#include <vector>

#include "Epitope.h"

class Simulation;
class DatabaseManager;
class BCell;
class PlasmaCell;

class Antigen
{
public:
	Antigen(
		Simulation & sim, SimParameters & params, int64_t id,
		zppsim::rng_t & rng,
		zppdb::Database * dbPtr,
		zppdb::Table<AntigenNeighborRow> * neighborTablePtr,
		zppdb::Table<AntigenEnergyRow> * energyTablePtr
	);
	
	double getLogAffinity(
		BCell const & bCell, int64_t i, int64_t j,
		zppsim::rng_t & rng
	) const;
//	double getLogAffinity(PlasmaCell const & pCell, int64_t i, int64_t j) const;
	double getMaxLogAffinity(
		BCell const & bCell, int64_t & iOut, int64_t & jOut,
		zppsim::rng_t & rng
	) const;
	double getEnergy(BCell const & bCell, int64_t i, int64_t j, zppsim::rng_t & rng) const;
	std::vector<std::vector<double>> getEnergies(BCell const & bCell, zppsim::rng_t & rng) const;
	
	void mutate(zppsim::rng_t & rng);
	
	int64_t id;
	
	int64_t const epitopeRows;
	int64_t const epitopeCols;
	int64_t const nLoci;
	
	double const affinityA;
	double const affinityB;
	
	double const quantity;
	
private:
	zppdb::Database * dbPtr;
	zppdb::Table<AntigenNeighborRow> * neighborTablePtr;
	zppdb::Table<AntigenEnergyRow> * energyTablePtr;
	
	std::vector<std::unique_ptr<Epitope>> epitopesVec;
	std::vector<std::vector<Epitope *>> epitopes;
};

#endif /* defined(__bcellmodel__Antigen__) */
