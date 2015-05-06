/**
 Implementation for class representing germinal center (location of B-cell affinity maturation).
*/

#include "GC.h"
#include "BCell.h"
#include "zppsim_random.hpp"
#include "DatabaseManager.h"

#include <cassert>
#include <iostream>
#include <unordered_map>
#include <functional>
#include "zppsim_random.hpp"

using namespace std;
using namespace zppdb;
using namespace zppsim;

static vector<double> normalizeExp(vector<double> const & logVals)
{
	vector<double> vals = logVals;
	double maxLogVal = -std::numeric_limits<double>::infinity();
	for(double logVal : logVals) {
		if(logVal > maxLogVal) {
			maxLogVal = logVal;
		}
	}
	if(maxLogVal > -std::numeric_limits<double>::infinity()) {
		for(int64_t i = 0; i < vals.size(); i++) {
			if(!std::isinf(vals[i])) {
				vals[i] -= maxLogVal;
			}
		}
	}
	for(int64_t i = 0; i < vals.size(); i++) {
		vals[i] = exp(vals[i]);
	}
	return vals;
}

GC::GC(Host & host, vector<BCell *> const & initCells, int64_t id, zppsim::rng_t & rng) :
	hostPtr(&host), dbmPtr(host.dbmPtr),
	p(host.p), id(id), alive(true)
{
	cerr << "GC " << id << "..." << endl;
	
	cells.reserve(p->gcs.nTotalCells);
	for(BCell * cellPtr : initCells) {
		cerr << "  " << cellPtr->toString(p->sequences.alphabet) << endl;
		
		BCell * newCellPtr = new BCell(host.getNextCellId(), *cellPtr);
		host.writeBCellToDatabase(*this, "gcfd_seed", *newCellPtr);
		cells.push_back(unique_ptr<BCell>(newCellPtr));
	}
}

void GC::runRound(
	zppsim::rng_t & rng, Host & host, Antigen & antigen,
		std::vector<BCell *> & gcfdCellsForMemoryExport,
		std::vector<BCell *> & gcfdCellsForPlasmaExport,
		double pMemExport
) {
	assert(alive);
	assert(pMemExport >= 0.0 && pMemExport <= 1.0);
	
	int64_t maxOffspring = 1 << p->gcs.nDivisionsPerRound;
	if(p->gcs.nDivisionsPerRound == 1) {
		assert(maxOffspring == 2);
	}
	else if(p->gcs.nDivisionsPerRound == 2) {
		assert(maxOffspring == 4);
	}
	
	// Calculate number of cells for this round
	int64_t nCells;
	if(p->gcs.exponentialGrowthActive) {
		nCells = std::min(
			int64_t(cells.size() * maxOffspring),
			int64_t(p->gcs.nTotalCells)
		);
	}
	else {
		nCells = p->gcs.nTotalCells;
	}
//	cerr << "# cells for round: " << nCells << endl;
	
	// New generation of cells
	std::vector<BCell *> newCells;
	if(p->gcs.useNullReplicationModel.present() && p->gcs.useNullReplicationModel) {
		newCells = sampleNewCellsNullModel(nCells, maxOffspring, host, rng);
	}
	else {
		newCells = sampleNewCells(nCells, maxOffspring, host, rng);
	}
	
	// Replace old cells with new cells
	cells.clear();
	for(BCell * newCell : newCells) {
		cells.emplace_back(newCell);
		if(p->gcs.shouldTerminateBelowThreshold) {
			double cellEffLogAff = host.maxEffectiveLogAffinity(*cells.back(), rng);
			if(cellEffLogAff < p->gcs.logAffinityTerminationThreshold) {
//				cerr << "CELL BELOW THRESHOLD" << endl;
				cells.pop_back();
			}
		}
	}
	
	// Terminate if the maximum affinity across all cells, across all epitopes, is less than the termination threshold.
	if(p->gcs.shouldTerminateBelowThreshold && cells.size() == 0) {
//		cerr << "NO NEW CELLS!" << endl;
		fprintf(stderr, "GC %llu terminating (all cells below threshold)\n", id);
		alive = false;
		return;
	}
	
	// Choose cells for memory export based on pMemExport (which is calculated based on the round number in Host.cpp).
	for(int64_t i : zppsim::drawMultipleBernoulli(rng, int64_t(cells.size()), pMemExport)) {
		gcfdCellsForMemoryExport.push_back(cells[i].get());
	}

	// Export plasma cells that are >= PLASMA_EXPORT_LOG_AFFINITY_THRESHOLD with probability P_PLASMA_EXPORT
	bernoulli_distribution flipCoin(p->gcs.pPlasmaExport);
	for(auto & cellUPtr : cells) {
		if(host.maxEffectiveLogAffinity(*cellUPtr, rng)
			>= p->gcs.plasmaExportLogAffinityThreshold
		) {
			if(flipCoin(rng)) {
				gcfdCellsForPlasmaExport.push_back(cellUPtr.get());
			}
		}
	}
}

std::vector<BCell *> GC::sampleNewCells(int64_t nCells, int64_t maxOffspring, Host & host, zppsim::rng_t & rng)
{
	std::vector<BCell *> newCells;
	newCells.reserve(p->gcs.nTotalCells);
	
	// Offspring counts (capped at 2^nDivisionsPerRound) and relative prob.
	vector<int64_t> offspringCounts(cells.size());
	vector<double> parentLogProbs(cells.size());
	for(int64_t i = 0; i < cells.size(); i++) {
		offspringCounts[i] = 0;
		parentLogProbs[i] = host.maxEffectiveLogAffinity(*cells[i], rng);
//		cerr << "parent log prob " << i << ": " << parentLogProbs[i] << endl;
	}
	
	vector<double> parentProbs = normalizeExp(parentLogProbs);
	for(int64_t i = 0; i < cells.size(); i++) {
//		cerr << "parent rel. prob " << i << ": " << parentProbs[i] << endl;
	}
	
	// Create new cells by sampling parent cells proportional to their affinity,
	// capping the number of offspring per cell at 2^nDivisionsPerRound
	for(int64_t i = 0; i < nCells; i++) {
		int64_t index = sampleDiscreteLinearSearch(rng, parentProbs);
//		cerr << "parent index " << index << endl;
		
		BCell * newCellPtr = new BCell(
			host.getNextCellId(), *cells[index], rng, p->gcs.pMutation,
			host.seqCache
		);
//		host.writeBCellToDatabase(*this, "gcfd", *newCellPtr);
//		host.writeAffinityToDatabase("gcfd", *newCellPtr, rng);
		newCells.push_back(newCellPtr);
		
		offspringCounts[index] += 1;
//		cerr << "offspringCounts " << index << " = " << offspringCounts[index] << endl;
		if(offspringCounts[index] == maxOffspring) {
			parentProbs[index] = 0.0;
		}
	}
	
	return newCells;
}

std::vector<BCell *> GC::sampleNewCellsNullModel(int64_t nCells, int64_t maxOffspring, Host & host, zppsim::rng_t & rng)
{
	std::vector<BCell *> newCells;
	newCells.reserve(p->gcs.nTotalCells);
	
	// Offspring counts (capped at 2^nDivisionsPerRound) and relative prob.
	vector<int64_t> offspringCounts(cells.size());
	for(int64_t i = 0; i < cells.size(); i++) {
		offspringCounts[i] = 0;
	}
	
	// Create new cells by sampling parent cells proportional to their affinity,
	// capping the number of offspring per cell at 2^nDivisionsPerRound
	uniform_int_distribution<int64_t> parentIndexDist(0, cells.size() - 1);
	for(int64_t i = 0; i < nCells; i++) {
		int64_t parentIndex;
		do {
			parentIndex = parentIndexDist(rng);
		}
		while(offspringCounts[parentIndex] == maxOffspring);
		
		BCell * newCellPtr = new BCell(
			host.getNextCellId(), *cells[parentIndex], rng, p->gcs.pMutation,
			host.seqCache
		);
//		host.writeBCellToDatabase(*this, "gcfd", *newCellPtr);
//		host.writeAffinityToDatabase("gcfd", *newCellPtr, rng);
		newCells.push_back(newCellPtr);
		
		offspringCounts[parentIndex] += 1;
	}
	
	return newCells;
}

bool GC::isAlive()
{
	return alive;
}

void GC::writeGCFDCellsToDatabase(Antigen const & antigen, int64_t round)
{
	if(!dbmPtr->tableExists(dbmPtr->gcfdCellsTable)) {
		return;
	}
	
	GCFDCellRow row;
	row.infection_id = antigen.id;
	row.gc_id = id;
	row.gc_round = round;
	for(auto & cellUPtr : cells) {
		row.cell_id = cellUPtr->getId();
		dbmPtr->insert(dbmPtr->TABLE(gcfdCells), row);
	}
}

void GC::writeEffectiveAffinityToDatabase(
	Host & host,
	Antigen const & antigen,
	int64_t round,
	zppsim::rng_t & rng
) {
	if(!dbmPtr->tableExists(dbmPtr->effectiveAffinityTable)) {
		return;
	}
	
	EffectiveAffinityRow row;
	row.infection_id = antigen.id;
	if(round == -1) {
		row.gc_round.setNull();
	}
	else {
		row.gc_round = round;
	}
	for(auto & cellUPtr : cells) {
		row.cell_id = cellUPtr->getId();
		for(int64_t i = 0; i < p->epitopes.rows; i++) {
			row.epitope_row = i;
			for(int64_t j = 0; j < p->epitopes.cols; j++) {
				row.epitope_column = j;
				double logAffinity = host.effectiveLogAffinity(i, j, *cellUPtr, rng);
				row.log_affinity = logAffinity;
				dbmPtr->insert(dbmPtr->TABLE(effectiveAffinity), row);
			}
		}
	}
}
