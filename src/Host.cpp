/**
	Implementation for class representing host (human individual).
*/

#include "Host.h"

#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <cassert>

#include "shared.h"
#include "Simulation.h"
#include "SimParameters.h"
#include "Antigen.h"
#include "GC.h"
#include "DatabaseManager.h"
#include "Simulation.h"
#include "NaiveCellGenerator.h"
#include <cassert>

#include <boost/math/distributions/normal.hpp>

using namespace std;
using namespace zppdb;

static double max(vector<vector<double>> const & mat) {
	double maxVal = -std::numeric_limits<double>::infinity();
	for(auto & vec : mat) {
		for(auto & val : vec) {
			if(val > maxVal) {
				maxVal = val;
			}
		}
	}
	return maxVal;
}

static double logSumExp(vector<vector<double>> const & mat)
{
	double maxLogVal = max(mat);
	double sumExp = 0.0;
	for(auto & vec : mat) {
		for(double logVal : vec) {
			sumExp += exp(logVal - maxLogVal);
		}
	}
	return maxLogVal + log(sumExp);
}

Host::Host(Simulation & sim) :
	epitopeRows(sim.p->epitopes.rows),
	epitopeCols(sim.p->epitopes.cols),
	simPtr(&sim), time(0.0), p(sim.p), dbmPtr(&sim.dbm),
	antigenPtr(nullptr), nextCellId(1),
	nNaiveBCells(sim.p->nBCellsNaive),
	naiveCellGen(*p, sim.db, dbmPtr->naiveAllelesTable, sim.rng),
	logAffinitiesMasking(epitopeRows, std::vector<double>(epitopeCols))
{
	clock_t clockStart = clock();
	fprintf(stderr, "Starting host initialization...\n");
	
	if(p->gcs.useOldSeedMethod) {
		// Create random B cells
		initializeNaiveCells(sim.rng);
	}
	
	clock_t clockEnd = clock();
	fprintf(stderr, "Host initialization finished (elapsed time: %.3f s)\n", elapsed(clockStart, clockEnd));
	
	writeNaiveBCellsToDatabase();
}

void Host::initializeNaiveCells(zppsim::rng_t & rng)
{
	naiveBCells.reserve(p->nBCellsNaive);
	
	for(int64_t i = 0; i < p->nBCellsNaive; i++) {
//		cerr << "GENERATING CELL..." << '\n';
		naiveBCells.push_back(unique_ptr<BCell>(naiveCellGen.generateCell(getNextCellId(), rng, seqCache)));
//		cerr << "SEQUENCE: " << naiveBCells[i]->toString(p->alphabet) << '\n';
	}
}

uint32_t Host::getNextCellId()
{
	assert(nextCellId < std::numeric_limits<uint32_t>::max());
	return nextCellId++;
}

void Host::runInfection(Antigen & antigen, zppsim::rng_t & rng)
{
	simPtr->db.beginTransaction();
	
	gcStarted = false;
	
	antigenPtr = &antigen;
	
	// Step time: reduce concentration of plasma cells, discard low-concentration cells
	if(antigen.id > 0) {
		assert(p->infections[antigen.id - 1].delay > p->gcs.nRounds * p->gcs.roundDuration);
		stepTime(p->infections[antigen.id - 1].delay - p->gcs.nRounds * p->gcs.roundDuration);
	}
	
	// Precalculate log-affinity of every plasma cell for every epitope
	cerr << "Precalculating log-affinity of every plasma cell" << endl;
	vector<vector<vector<double>>> logAffPlasma;
	for(int64_t ip = 0; ip < plasmaCells.size(); ip++) {
		PlasmaCell * cellPtr = (PlasmaCell *)plasmaCells[ip].get();
		logAffPlasma.emplace_back();
		for(int64_t i = 0; i < epitopeRows; i++) {
			logAffPlasma[ip].emplace_back();
			for(int64_t j = 0; j < epitopeCols; j++) {
				logAffPlasma[ip][i].push_back(antigen.getLogAffinity(*cellPtr, i, j, rng));
			}
		}
	}
	
	// Precalculate log-sum-affinity of every plasma cell across all epitopes
	cerr << "Adding affinities of plasma cells" << endl;
	vector<double> logSumAffPlasma;
	for(int64_t ip = 0; ip < plasmaCells.size(); ip++) {
		logSumAffPlasma.push_back(logSumExp(logAffPlasma[ip]));
	}
	
	// Calculate log-affinity of masking antibody
	for(int64_t i = 0; i < p->epitopes.rows; i++) {
		for(int64_t j = 0; j < p->epitopes.cols; j++) {
			logAffinitiesMasking[i][j] = getLogAffinityMasking(logAffPlasma, logSumAffPlasma, antigen, i, j, rng);
		}
	}
	
	// Write naive affinities
	writeAffinityToDatabase("naive", naiveBCells, rng);
	writeAffinityStatsToDatabase("naive", naiveBCells, rng);
	
	// Create GCs
	vector<unique_ptr<GC>> gcs;
	if(p->gcs.useOldSeedMethod) {
		gcs = createGCsOldMethod(antigen, rng);
	}
	else {
		gcs = createGCsNewMethod(antigen, rng);
	}
	simPtr->db.commitWithRetry(0.1, 100, std::cerr);
	
	// Perform rounds of SHM, selection, and memory cell creation
	gcStarted = true;
	gcRound = 0;
	while(gcRound < p->gcs.nRounds) {
		simPtr->db.beginTransaction();
		
		if(gcRound > 0) {
			// Update effective log-affinity functions
			for(int64_t i = 0; i < p->epitopes.rows; i++) {
				for(int64_t j = 0; j < p->epitopes.cols; j++) {
					logAffinitiesMasking[i][j] = getLogAffinityMasking(logAffPlasma, logSumAffPlasma, antigen, i, j, rng);
				}
			}
		}
		
		double pMemExpFirst = p->gcs.pMemoryExportFirstRound;
		double pMemExpLast = p->gcs.pMemoryExportLastRound;
		double pMemExport = pMemExpFirst +
			(pMemExpLast - pMemExpFirst) * (1.0 - double(gcRound) / double(p->gcs.nRounds - 1));
		
		clock_t clockStart = clock();
		fprintf(stderr, "Starting round %lld...\n", gcRound);
		
		simPtr->db.commitWithRetry(0.1, 100, std::cerr);
		
		int64_t nAlive = 0;
		for(int64_t i = 0; i < p->gcs.nGCs; i++) {
			simPtr->db.beginTransaction();
			if(gcs[i]->isAlive()) {
				nAlive++;
				vector<BCell *> gcfdCellsForMemoryExport;
				vector<BCell *> gcfdCellsForPlasmaExport;
				gcs[i]->runRound(
					rng, *this, antigen,
					gcfdCellsForMemoryExport, gcfdCellsForPlasmaExport,
					pMemExport
				);
				
				gcs[i]->writeGCFDCellsToDatabase(antigen, gcRound);
				gcs[i]->writeEffectiveAffinityToDatabase(
					*this, antigen, gcRound, rng
				);
				
				if(gcs[i]->isAlive()) {
					// Copy memory and plasma cells marked for export
					vector<BCell *> newMemoryCells;
					for(BCell * cellPtr : gcfdCellsForMemoryExport) {
						BCell * newCellPtr = new BCell(getNextCellId(), *cellPtr);
						memoryBCells.emplace_back(newCellPtr);
						newMemoryCells.push_back(newCellPtr);
					}
					vector<BCell *> newPlasmaCells;
					for(BCell * cellPtr : gcfdCellsForPlasmaExport) {
						PlasmaCell * newCellPtr = new PlasmaCell(
							getNextCellId(),
							*cellPtr,
							antigen,
							p->shortLivedPlasma.initialLogC,
							p->longLivedPlasma.initialLogC
						);
						plasmaCells.emplace_back(newCellPtr);
						newPlasmaCells.push_back(newCellPtr);
						
						// Add log-affinity, log-sum-affinity values for new plasma cells
						logAffPlasma.emplace_back();
						for(int64_t i = 0; i < epitopeRows; i++) {
							logAffPlasma.back().emplace_back();
							for(int64_t j = 0; j < epitopeCols; j++) {
								logAffPlasma.back()[i].push_back(antigen.getLogAffinity(*newCellPtr, i, j, rng));
							}
						}
						logSumAffPlasma.push_back(logSumExp(logAffPlasma.back()));
						
					}
					
					writeBCellsToDatabase(*gcs[i], "memory", newMemoryCells);
					writeBCellsToDatabase(*gcs[i], "plasma", newPlasmaCells);
				}
			}
			
			simPtr->db.commitWithRetry(0.1, 100, std::cerr);
		}
		simPtr->db.beginTransaction();
		writeGCStatsToDatabase(nAlive);
		writeAffinityStatsToDatabase("memory", memoryBCells, rng);
		writeAffinityStatsToDatabase("plasma", plasmaCells, rng);
		clock_t clockEnd = clock();
		fprintf(stderr, "Round %lld complete (elapsed time: %.3f s).\n", gcRound, elapsed(clockStart, clockEnd));
		
		writePlasmaConcentrationToDatabase();
		
//		verifyPerGCRound(gcs);
		
		stepTime(p->gcs.roundDuration);
		
		seqCache.writeToDatabase(simPtr->db, dbmPtr->sequencesTable, simPtr->p->sequences.alphabet);
		seqCache.collectGarbage();
		
		simPtr->db.commitWithRetry(0.1, 100, std::cerr);
		
		gcRound++;
	}
	simPtr->db.beginTransaction();
	writeAffinityToDatabase("memory", memoryBCells, rng);
	writeAffinityToDatabase("plasma", plasmaCells, rng);
	simPtr->db.commitWithRetry(0.1, 100, std::cerr);
}

vector<unique_ptr<GC>> Host::createGCsOldMethod(
	Antigen const & antigen, zppsim::rng_t & rng
) {
	// Select seed B cells and write to database
	vector<BCell *> seedBCells = identifyGCSeedCandidates(antigen, rng);
	writeGCSeedCandidatesToDatabase(seedBCells);
	
	// Create germinal centers from seed B cells
	uniform_int_distribution<int64_t> unifSeedBCell(0, int64_t(seedBCells.size() - 1));
	vector<unique_ptr<GC>> gcs(p->gcs.nGCs);
	for(int64_t i = 0; i < p->gcs.nGCs; i++) {
		vector<BCell *> initBCells;
		for(int64_t j = 0; j < p->gcs.nSeedCells; j++) {
			initBCells.push_back(seedBCells[unifSeedBCell(rng)]);
		}
		gcs[i] = unique_ptr<GC>(new GC(*this, initBCells, i, rng));
		writeGCSeedsToDatabase(*gcs[i], initBCells);
	}
	return gcs;
}

vector<unique_ptr<GC>> Host::createGCsNewMethod(
	Antigen const & antigen, zppsim::rng_t & rng
) {
	int64_t epRows = p->epitopes.rows;
	int64_t epCols = p->epitopes.cols;
	
	// Calculate log-affinity thresholds for each epitope
	vector<vector<double>> logAffLower(epRows, vector<double>(epCols));
	vector<vector<double>> logAffUpper(epRows, vector<double>(epCols));
	if(antigen.id == 0) {
		for(int64_t i = 0; i < epRows; i++) {
			for(int64_t j = 0; j < epCols; j++) {
				double logAffTh = logAffinityQuantile(1.0 - p->gcs.fNaiveSeed, i, j);
				logAffLower[i][j] = logAffTh;
				logAffUpper[i][j] = std::numeric_limits<double>::infinity();
				fprintf(stderr, "lal, lau = %lf, %lf\n", logAffLower[i][j], logAffUpper[i][j]);
			}
		}
	}
	else {
		for(int64_t i = 0; i < epRows; i++) {
			for(int64_t j = 0; j < epCols; j++) {
				double logAffLowerRaw = logAffinityQuantile(p->gcs.seedQuantileLower, i, j);
				double logAffUpperRaw = logAffinityQuantile(p->gcs.seedQuantileUpper, i, j);
				// C_Ab*F_Ab,j)*E_j / [Ag]
				if(logAffLowerRaw == -std::numeric_limits<double>::infinity()) {
//					logAffLower[i][j] = logAffLowerRaw;
					assert(false);
				}
				else {
//					logAffLower[i][j] = logAffLowerRaw - logAffinityOffsets[i][j];
					assert(false);
				}
				if(logAffUpperRaw == std::numeric_limits<double>::infinity()) {
//					logAffUpper[i][j] = logAffUpperRaw;
					assert(false);
				}
				else {
//					logAffUpper[i][j] = logAffUpperRaw - logAffinityOffsets[i][j];
					assert(false);
				}
				fprintf(stderr, "lal, lau = %lf, %lf\n", logAffLower[i][j], logAffUpper[i][j]);
			}
		}
	}
	
	// Calculate expected number of naive cells within bounds
	double nNaiveInBounds[epRows][epCols];
	double totalInBounds = 0.0;
	for(int64_t i = 0; i < epRows; i++) {
		for(int64_t j = 0; j < epCols; j++) {
			double quantileLower = logAffinityCDF(logAffLower[i][j], i, j);
			double quantileUpper = logAffinityCDF(logAffUpper[i][j], i, j);
			
			fprintf(stderr, "%lld, %lld: ql, qu = %lf, %lf\n", i, j, quantileLower, quantileUpper);
			
			nNaiveInBounds[i][j] = (quantileUpper - quantileLower) * nNaiveBCells;
			totalInBounds += nNaiveInBounds[i][j];
			
			if(antigen.id == 0) {
				assert(abs(quantileLower - (1.0 - p->gcs.fNaiveSeed)) < 0.0001);
			}
		}
	}
	
	// Assemble list of all memory cells in bounds (double-count if in bounds for multiple epitopes)
	vector<BCell *> memCellsInBounds;
	for(auto & memCellPtr : memoryBCells) {
		for(int64_t i = 0; i < epRows; i++) {
			for(int64_t j = 0; j < epCols; j++) {
				if(antigen.getLogAffinity(*memCellPtr, i, j, rng)) {
					memCellsInBounds.push_back(memCellPtr.get());
				}
			}
		}
	}
	totalInBounds += memCellsInBounds.size();
	
	// Create germinal centers from seed B cells
	vector<unique_ptr<GC>> gcs(p->gcs.nGCs);
	
	uniform_real_distribution<double> cellTypeUnifDist(0.0, totalInBounds);
	uniform_int_distribution<int64_t> memIndexDist(0, memCellsInBounds.size() - 1);
	
	for(int64_t gcId = 0; gcId < p->gcs.nGCs; gcId++) {
		vector<unique_ptr<BCell>> tmpNaiveCells;
		vector<BCell *> initBCells;
		for(int64_t cellNum = 0; cellNum < p->gcs.nSeedCells; cellNum++) {
			// Draw which kind of cell (naive within bounds for a particular epitope
			// or memory cell) using linear search
			double unifVal = cellTypeUnifDist(rng);
			double cumSum = 0.0;
			bool done = false;
			for(int64_t i = 0; !done && i < epRows; i++) {
				for(int64_t j = 0; !done && j < epCols; j++) {
					cumSum += nNaiveInBounds[i][j];
					if(unifVal <= cumSum) {
						fprintf(stderr, "generating naive cell within range for ep %lld, %lld\n", i, j);
						tmpNaiveCells.emplace_back(createNaiveCell(
							antigen, i, j, logAffLower[i][j], logAffUpper[i][j], rng
						));
						initBCells.push_back(tmpNaiveCells.back().get());
						done = true;
					}
					fprintf(stderr, "%lld, %lld: %lf expected naive cells in range\n", i, j, nNaiveInBounds[i][j]);
				}
			}
			if(!done) {
				assert(unifVal > cumSum);
				int64_t memIndex = memIndexDist(rng);
				fprintf(stderr, "using memory cell %lld\n", memIndex);
				initBCells.push_back(memCellsInBounds[memIndex]);
			}
		}
		writeNaiveBCellsToDatabase(tmpNaiveCells);
		gcs[gcId] = unique_ptr<GC>(new GC(*this, initBCells, gcId, rng));
		writeGCSeedsToDatabase(*gcs[gcId], initBCells);
	}
	return gcs;
}

BCell * Host::createNaiveCell(
	Antigen const & antigen, int64_t epRow, int64_t epCol,
	double logAffMin, double logAffMax, zppsim::rng_t & rng
) {
	BCell * bCellPtr;
	while(true) {
		bCellPtr = naiveCellGen.generateCell(getNextCellId(), rng, seqCache);
		double logAff = antigen.getLogAffinity(*bCellPtr, epRow, epCol, rng);
//		cerr << "log affinity: " << logAff << endl;
		if(logAff >= logAffMin && logAff <= logAffMax) {
			break;
		}
		delete bCellPtr;
	}
	return bCellPtr;
}

double Host::logAffinityQuantile(double quantile, int64_t epRow, int64_t epCol)
{
	if(quantile == 0.0) {
		return -std::numeric_limits<double>::infinity();
	}
	if(quantile == 1.0) {
		return std::numeric_limits<double>::infinity();
	}
	
	double a = p->epitopes.affinityA;
	double b = p->epitopes.affinityB;
	
	boost::math::normal norm(p->epitopes.energyMean[epRow][epCol], 1.0);
	double energyThreshold = boost::math::quantile(norm, 1.0 - quantile);
//	cerr << "energy threshold: " << energyThreshold << endl;
	
	return a - b * energyThreshold;
}

double Host::logAffinityCDF(double logAff, int64_t epRow, int64_t epCol)
{
	if(logAff == -std::numeric_limits<double>::infinity()) {
		return 0.0;
	}
	if(logAff == std::numeric_limits<double>::infinity()) {
		return 1.0;
	}
	
	double a = p->epitopes.affinityA;
	double b = p->epitopes.affinityB;
	
	// logAff = a - b * E
	double E = (a - logAff) / b;
	
	boost::math::normal norm(p->epitopes.energyMean[epRow][epCol], 1.0);
	return boost::math::cdf(boost::math::complement(norm, E));
}

void Host::stepTime(double delay)
{
	double longLivedLogFactor = delay * log(p->longLivedPlasma.decayRate);
	double shortLivedLogFactor = delay * log(p->shortLivedPlasma.decayRate);
	
	// Adjust concentrations and track removed plasma cells
	vector<int64_t> removedIndices;
	for(size_t i = 0; i < plasmaCells.size(); i++) {
		PlasmaCell * cellPtr = (PlasmaCell *)plasmaCells[i].get();
		assert(cellPtr != nullptr);
		
		double logCLongLived = cellPtr->getLogLongLivedConcentration();
		cellPtr->setLogLongLivedConcentration(logCLongLived + longLivedLogFactor);
		
		double logCShortLived = cellPtr->getLogShortLivedConcentration();
		cellPtr->setLogShortLivedConcentration(logCShortLived + shortLivedLogFactor);
		
		if(logCLongLived < p->longLivedPlasma.discardLogC
			&& logCShortLived < p->shortLivedPlasma.discardLogC
		) {
			removedIndices.push_back(i);
			plasmaCells[i] = nullptr;
		}
	}
	
	// Move existing cells from end into vacated positions
	int64_t oldSize = plasmaCells.size();
	
	int64_t nRemoved = 0;
	for(int64_t removedIndex : removedIndices) {
		// Remove any nulls at the end to begin with
		while(plasmaCells.back() == nullptr) {
			plasmaCells.pop_back();
			nRemoved++;
		}
		
		// We can't be exactly at the removed index because otherwise
		// it would have been removed by the nullptr for loop above
		assert(removedIndex != plasmaCells.size() - 1);
		
		// If we've removed any nulls greedily from the end, we'll be
		// done before we've swapped into all the removed indices
		if(removedIndex >= plasmaCells.size()) {
	 		break;
		}
		
		// But normally, we'll swap a non-null from the end into the removed
		// index position.
		plasmaCells[removedIndex] = std::move(plasmaCells.back());
		plasmaCells.pop_back();
		nRemoved++;
	}
	assert(nRemoved == removedIndices.size());
	
	/*for(int64_t removedIndex : removedIndices) {
		if(plasmaCells.back() != nullptr) {
			plasmaCells[removedIndex] = std::move(plasmaCells.back());
		}
		plasmaCells.pop_back();
	}*/

	assert(plasmaCells.size() == oldSize - removedIndices.size());
	
	time += delay;
	fprintf(stderr, "time advanced %lf\n", delay);
	fprintf(stderr, "plasma concentrations adjusted by log-factor of %lg (long-lived), %lg (short-lived)\n", longLivedLogFactor, shortLivedLogFactor);
	fprintf(stderr, "t = %lf\n", time);
}

vector<BCell *> Host::identifyGCSeedCandidates(
	Antigen const & antigen,
	zppsim::rng_t & rng
) {
	vector<BCell *> seedBCells;
	
	// Construct vector of naive and memory B cells and their affinities for all epitopes;
	// sort in descending order by affinity.
	// Exclude memory cells if gcs.seedExcludeMemory == true.
	vector<BCellDoublePair> baPairs;
	for(auto & bCellPtr : naiveBCells) {
		for(int64_t i = 0; i < epitopeRows; i++) {
			for(int64_t j = 0; j < epitopeCols; j++) {
				baPairs.emplace_back(
					bCellPtr.get(),
					effectiveLogAffinity(i, j, *bCellPtr, rng)
				);
			}
		}
	}
	
	if(!p->gcs.seedExcludeMemory) {
		for(auto & bCellPtr : memoryBCells) {
			for(int64_t i = 0; i < epitopeRows; i++) {
				for(int64_t j = 0; j < epitopeCols; j++) {
					baPairs.emplace_back(
						bCellPtr.get(),
						effectiveLogAffinity(i, j, *bCellPtr, rng)
					);
				}
			}
		}
	}
	sort(baPairs.begin(), baPairs.end(), std::greater<BCellDoublePair>());
	assert(baPairs.size() > 0);
	
	if(antigen.id == 0) {
		// Just take the top F_NAIVE_SEED fraction of cell/affinity pairs
		int64_t nNaiveSeed = int64_t(ceil(p->gcs.fNaiveSeed * baPairs.size()));
		assert(nNaiveSeed > 0);
		seedBCells.reserve(nNaiveSeed);
		for(int64_t i = 0; i < nNaiveSeed; i++) {
			seedBCells.push_back(baPairs[i].cellPtr);
		}
		assert(seedBCells.size() == nNaiveSeed);
	}
	else {
		// Take the quantile range GC_SEED_LOWER_QUANTILE, GC_SEED_UPPER_QUANTILE of cell/affinity pairs
		int64_t startIndex = min(
			int64_t(baPairs.size() - 1),
			int64_t(floor((1.0 - p->gcs.seedQuantileUpper) * baPairs.size()))
		);
		int64_t endIndex = max(
			startIndex + 1,
			int64_t(ceil((1.0 - p->gcs.seedQuantileLower) * baPairs.size()))
		);
		assert(startIndex >= 0);
		assert(endIndex > startIndex);
		assert(endIndex <= baPairs.size());
		seedBCells.reserve(endIndex - startIndex);
		for(int64_t i = startIndex; i < endIndex; i++) {
			seedBCells.push_back(baPairs[i].cellPtr);
		}
	}
	assert(seedBCells.size() > 0);
	return seedBCells;
}

double Host::getLogAffinityMasking(
	std::vector<std::vector<std::vector<double>>> const & logAffsPlasma, std::vector<double> logSumAffsPlasma,
	Antigen const & antigen, int64_t i, int64_t j, zppsim::rng_t & rng
) {
	if(plasmaCells.size() == 0) {
		return -std::numeric_limits<double>::infinity();
	}
	else {
//		double logAt = log(p->infections[antigen.id].antigenQuantity);
		double logFCeffMax = -std::numeric_limits<double>::infinity();
		double logFMasking = -std::numeric_limits<double>::infinity();
		
		// Find plasma cell with highest (affinity * effective concentration),
		// where effective concentration = affinity / (sum of affinities) * concentration
		for(int64_t ip = 0; ip < plasmaCells.size(); ip++) {
			PlasmaCell *cellPtr = (PlasmaCell *)plasmaCells[ip].get();
			double logFCeff = cellPtr->getLogConcentration()
				+ 2 * logAffsPlasma[ip][i][j]
				- logSumAffsPlasma[ip];
			if(logFCeff > logFCeffMax) {
				logFCeffMax = logFCeff;
				logFMasking = logAffsPlasma[ip][i][j];
			}
		}
		assert(logFCeffMax != -std::numeric_limits<double>::infinity());
		assert(logFMasking != -std::numeric_limits<double>::infinity());
		
		return logFMasking;
	}
}

void Host::writeNaiveBCellsToDatabase()
{
	writeNaiveBCellsToDatabase(naiveBCells);
}

void Host::writeNaiveBCellsToDatabase(vector<unique_ptr<BCell>> const & cells)
{
	if(!dbmPtr->tableExists(dbmPtr->bcellsTable)) {
		return;
	}
	
	BCellRow row;
	row.type = "naive";
	
	for(auto & cellUPtr : cells) {
		row.cell_id = cellUPtr->getId();
		row.sequence_id = cellUPtr->getSequenceId();
		dbmPtr->insert(dbmPtr->TABLE(bcells), row);
	}
}

void Host::writeBCellsToDatabase(
	GC const & gc, string const & type, vector<BCell *> const & cells
) {
	if(!dbmPtr->tableExists(dbmPtr->bcellsTable)) {
		return;
	}
	
	BCellRow row;
	row.type = type;
	row.infection_id = antigenPtr->id;
	row.gc_id = gc.id;
	if(gcStarted) {
		row.gc_round = gcRound;
	}
	else {
		row.gc_round.setNull();
	}
	for(auto cellPtr : cells) {
		row.cell_id = cellPtr->getId();
		row.parent_id = cellPtr->getParentId();
		row.sequence_id = cellPtr->getSequenceId();
		dbmPtr->insert(dbmPtr->TABLE(bcells), row);
	}
}

void Host::writeBCellToDatabase(GC const & gc, std::string const & type, BCell const & cell)
{
	if(!dbmPtr->tableExists(dbmPtr->bcellsTable)) {
		return;
	}
	
	BCellRow row;
	row.type = type;
	row.infection_id = antigenPtr->id;
	row.gc_id = gc.id;
	if(gcStarted) {
		row.gc_round = gcRound;
	}
	else {
		row.gc_round.setNull();
	}
	row.cell_id = cell.getId();
	row.parent_id = cell.getParentId();
	row.sequence_id = cell.getSequenceId();
	dbmPtr->insert(dbmPtr->TABLE(bcells), row);
}

void Host::writeAffinityToDatabase(std::string const & type, BCell & cell, zppsim::rng_t & rng)
{
	if(!p->affinityOutputTypesEnabled[type]) {
		return;
	}
	if(!dbmPtr->tableExists(dbmPtr->affinityTable)) {
		return;
	}
	
	AffinityRow row;
	row.type = type;
	row.infection_id = antigenPtr->id;
	row.cell_id = cell.getId();
	for(int64_t i = 0; i < p->epitopes.rows; i++) {
		row.epitope_row = i;
		for(int64_t j = 0; j < p->epitopes.cols; j++) {
			row.epitope_column = j;
			double energy = antigenPtr->getEnergy(cell, i, j, rng);
			double logAffinity = antigenPtr->getLogAffinity(cell, i, j, rng);
			row.energy = energy;
			row.log_affinity = logAffinity;
			dbmPtr->insert(dbmPtr->TABLE(affinity), row);
		}
	}
}

void Host::writeAffinityToDatabase(
	std::string const & type, vector<unique_ptr<BCell>> const & cells, zppsim::rng_t & rng
) {
	if(!p->affinityOutputTypesEnabled[type]) {
		return;
	}
	if(!dbmPtr->tableExists(dbmPtr->affinityTable)) {
		return;
	}
	
	AffinityRow row;
	row.type = type;
	row.infection_id = antigenPtr->id;
	for(auto & cellPtr : cells) {
		row.cell_id = cellPtr->getId();
		for(int64_t i = 0; i < p->epitopes.rows; i++) {
			row.epitope_row = i;
			for(int64_t j = 0; j < p->epitopes.cols; j++) {
				row.epitope_column = j;
				double energy = antigenPtr->getEnergy(*cellPtr, i, j, rng);
				double logAffinity = antigenPtr->getLogAffinity(*cellPtr, i, j, rng);
				row.energy = energy;
				row.log_affinity = logAffinity;
				dbmPtr->insert(dbmPtr->TABLE(affinity), row);
			}
		}
	}
}

void Host::writeGCStatsToDatabase(int64_t nAlive)
{
	GCStatsRow row;
	row.infection_id = antigenPtr->id;
	row.gc_round = gcRound;
	row.n_active_gcs = nAlive;
	dbmPtr->insert(dbmPtr->TABLE(gcStats), row);
}

void Host::writeAffinityStatsToDatabase(
	std::string const & type, vector<unique_ptr<BCell>> const & cells, zppsim::rng_t & rng
) {
	if(!p->affinityOutputTypesEnabled[type]) {
		return;
	}
	if(!dbmPtr->tableExists(dbmPtr->affinityStatsTable)) {
		return;
	}
	
	int64_t cols = p->epitopes.cols;
	int64_t rows = p->epitopes.rows;
	
	// Calculate sum, sum-of-squared affinities, matching cell counts for all epitopes
	double affSum[rows][cols];
	double affSumSq[rows][cols];
	double maxAffSum[rows][cols];
	double maxAffSumSq[rows][cols];
	int64_t counts[rows][cols];
	int64_t maxCounts[rows][cols];
	for(int64_t i = 0; i < p->epitopes.rows; i++) {
		for(int64_t j = 0; j < p->epitopes.cols; j++) {
			affSum[i][j] = 0.0;
			affSumSq[i][j] = 0.0;
			counts[i][j] = 0;
			maxAffSum[i][j] = 0.0;
			maxAffSumSq[i][j] = 0.0;
			maxCounts[i][j] = 0;
		}
	}
	
	// Log-affinity: matching for plasma; all cells for memory
	// TODO: check how to do log-affinity output
	for(auto & cellUPtr : cells) {
		for(int64_t i = 0; i < p->epitopes.rows; i++) {
			for(int64_t j = 0; j < p->epitopes.cols; j++) {
				counts[i][j]++;
				double logAffinity = antigenPtr->getLogAffinity(*cellUPtr, i, j, rng);
				affSum[i][j] += logAffinity;
				affSumSq[i][j] += logAffinity * logAffinity;
			}
		}
	}
	
	// Max-log-affinity: cells whose highest affinity is for each epitope
	for(auto & cellUPtr : cells) {
		int64_t row = std::numeric_limits<int64_t>::max();
		int64_t col = std::numeric_limits<int64_t>::max();
		double logAffinity = antigenPtr->getMaxLogAffinity(*cellUPtr, row, col, rng);
		
		assert(row != std::numeric_limits<int64_t>::max());
		assert(col != std::numeric_limits<int64_t>::max());
		
		maxAffSum[row][col] += logAffinity;
		maxAffSumSq[row][col] += logAffinity * logAffinity;
		maxCounts[row][col]++;
	}
	
	// Write mean, sd to database
	AffinityStatsRow row;
	row.type = type;
	row.infection_id = antigenPtr->id;
	if(gcStarted) {
		row.gc_round = gcRound;
	}
	else {
		row.gc_round.setNull();
	}
	for(int64_t i = 0; i < p->epitopes.rows; i++) {
		row.epitope_row = i;
		for(int64_t j = 0; j < p->epitopes.cols; j++) {
			row.epitope_column = j;
			
			row.n_cells = counts[i][j];
			if(counts[i][j] > 0) {
				double affMean = affSum[i][j] / counts[i][j];
				row.log_affinity_mean = affMean;
				
				double affSD = standardDeviationPopulation(affSum[i][j], affSumSq[i][j], uint32_t(counts[i][j]));
				row.log_affinity_sd = affSD;
			}
			else {
				row.log_affinity_mean.setNull();
				row.log_affinity_sd.setNull();
			}
			
			row.max_n_cells = maxCounts[i][j];
			if(maxCounts[i][j] > 0) {
				double maxAffMean = maxAffSum[i][j] / maxCounts[i][j];
				row.max_log_affinity_mean = maxAffMean;
				
				double maxAffSD = standardDeviationPopulation(maxAffSum[i][j], maxAffSumSq[i][j], uint32_t(maxCounts[i][j]));
				row.max_log_affinity_sd = maxAffSD;
			}
			else {
				row.max_log_affinity_mean.setNull();
				row.max_log_affinity_sd.setNull();
			}
			
			dbmPtr->insert(dbmPtr->TABLE(affinityStats), row);
		}
	}
}

void Host::writePlasmaConcentrationToDatabase()
{
	if(!dbmPtr->tableExists(dbmPtr->plasmaConcentrationsTable)) {
		return;
	}
	
	PlasmaConcRow row;
	row.infection_id = antigenPtr->id;
	row.gc_round = gcRound;
	for(auto & cellUniquePtr : plasmaCells) {
		PlasmaCell * cellPtr = (PlasmaCell *)cellUniquePtr.get();
		
		row.cell_id = cellPtr->getId();
		row.log_conc_short = cellPtr->getLogShortLivedConcentration();
		row.log_conc_long = cellPtr->getLogLongLivedConcentration();
		dbmPtr->insert(dbmPtr->TABLE(plasmaConcentrations), row);
	}
}

void Host::writeGCSeedCandidatesToDatabase(std::vector<BCell *> const & cells)
{
	if(!dbmPtr->tableExists(dbmPtr->gcSeedCandidatesTable)) {
		return;
	}
	
	GCSeedCandidateRow row;
	row.infection_id = antigenPtr->id;
	for(auto cellPtr : cells) {
		row.cell_id = cellPtr->getId();
		dbmPtr->insert(dbmPtr->TABLE(gcSeedCandidates), row);
	}
}

void Host::writeGCSeedsToDatabase(GC const & gc, std::vector<BCell *> const & cells)
{
	if(!dbmPtr->tableExists(dbmPtr->gcSeedsTable)) {
		return;
	}
	
	GCSeedRow row;
	row.infection_id = antigenPtr->id;
	row.gc_id = gc.id;
	for(auto cellPtr : cells) {
		row.cell_id = cellPtr->getId();
		dbmPtr->insert(dbmPtr->TABLE(gcSeeds), row);
	}
}

double Host::effectiveLogAffinity(int64_t epRow, int64_t epCol, BCell const & bCell, zppsim::rng_t & rng)
{
	double logFCell = antigenPtr->getLogAffinity(bCell, epRow, epCol, rng);
	double logFAntibody = logAffinitiesMasking[epRow][epCol];
	double FRatio = std::isinf(logFAntibody) ? 0.0 : exp(logFAntibody - logFCell);
	
	return logFCell - p->gcs.alphaMasking * FRatio;
}

double Host::maxEffectiveLogAffinity(BCell const & bCell, zppsim::rng_t & rng)
{
	double maxLA = -std::numeric_limits<double>::infinity();
	for(int64_t i = 0; i < epitopeRows; i++) {
		for(int64_t j = 0; j < epitopeCols; j++) {
			double LAij = effectiveLogAffinity(i, j, bCell, rng);
			if(LAij > maxLA) {
				maxLA = LAij;
			}
		}
	}
	return maxLA;
}

double Host::meanMaxAffinityMemory(Antigen & antigen, zppsim::rng_t & rng)
{
	return meanMaxAffinity(antigen, memoryBCells, rng);
}

double Host::meanMaxAffinity(Antigen & antigen, std::vector<std::unique_ptr<BCell>> & cells, zppsim::rng_t & rng)
{
	double sum = 0.0;
	for(auto & cell : cells) {
		int64_t i, j;
		sum += antigen.getMaxLogAffinity(*cell, i, j, rng);
	}
	return sum / cells.size();
}
