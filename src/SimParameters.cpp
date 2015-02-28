#include "SimParameters.h"
#include <cassert>
#include <limits>

using namespace std;

static void expandMatrix(Array<Array<Double>> & mat, int64_t rows, int64_t cols)
{
	// If the vector isn't big enough, its size will be expanded from 1
	// to match the required size.
	if(mat.size() < rows) {
		assert(mat.size() == 1);
		assert(mat[0].size() == 1);
		
		for(int64_t i = 0; i < rows; i++) {
			for(int64_t j = 0; j < cols; j++) {
				if(i == 0 && j == 0) {
					continue;
				}
				mat[i].push_back(double(mat[0][0]));
			}
		}
	}
	else {
		for(int64_t i = 0; i < rows; i++) {
			assert(mat[i].size() >= cols);
		}
	}
}

void verify(SimParameters & p)
{
	assert(p.randomSeed >= 0);
	assert(p.nBCellsNaive > 0);
	verify(p.epitopes);
	verify(p.gcs);
	verify(p.sequences);
	assert(p.infections.size() > 0);
	for(int64_t i = 0; i < p.infections.size(); i++) {
		verify(p.infections[i],
			i < p.infections.size() - 1 ? p.gcs.nRounds * p.gcs.roundDuration
				: -std::numeric_limits<double>::infinity()
		);
	}
	verify(p.shortLivedPlasma);
	verify(p.longLivedPlasma);
}

void verify(SegmentParams & p)
{
	assert(p.minLength > 0);
	assert(p.lengthDistribution.size() > 0);
	if(p.type == SegmentType::ALLELE()) {
		assert(p.alleleDistribution.size() > 0);
	}
}

void verify(ChainParams & p)
{
	assert(p.truncatedLength > 0);
	assert(p.segments.size() > 0);
	for(int64_t i = 0; i < p.segments.size(); i++) {
		verify(p.segments[i]);
	}
}

void verify(InfectionParams & p, double minDelay)
{
	assert(p.delay > minDelay);
	assert(p.antigenQuantity > 0.0);
}

void verify(EpitopeParams & p)
{
	assert(p.rows > 0);
	assert(p.cols > 0);
	assert(p.nInteractions >= 0);
	
	expandMatrix(p.pEnergyMutation, p.rows, p.cols);
	expandMatrix(p.pNeighborMutation, p.rows, p.cols);
	expandMatrix(p.energyMean, p.rows, p.cols);
	for(int64_t i = 0; i < p.rows; i++) {
		for(int64_t j = 0; j < p.cols; j++) {
			assert(p.pEnergyMutation[i][j] >= 0.0);
			assert(p.pEnergyMutation[i][j] <= 1.0);
			assert(p.pNeighborMutation[i][j] >= 0.0);
			assert(p.pNeighborMutation[i][j] <= 1.0);
			assert(p.energyMean[i][j] >= 0.0);
			assert(p.energyMean[i][j] <= 1.0);
		}
	}
}

void verify(SequenceParams & p)
{
	assert(p.alphabetSize > 1);
	assert(string(p.alphabet).size() >= p.alphabetSize);
	assert(p.chains.size() > 0);
//	nLoci = 0;
	for(int64_t i = 0; i < p.chains.size(); i++) {
		verify(p.chains[i]);
//		nLoci += chain.truncatedLength;
	}
}

void verify(GCParams & p)
{
	assert(p.nGCs >= 1);
	assert(p.nRounds >= 1);
	assert(p.roundDuration > 0.0);
	assert(p.fNaiveSeed > 0.0 && p.fNaiveSeed <= 1.0);
	assert(p.seedQuantileLower >= 0.0 && p.seedQuantileLower <= 1.0);
	assert(p.seedQuantileUpper > p.seedQuantileLower && p.seedQuantileUpper <= 1.0);
	assert(p.pMutation >= 0.0 && p.pMutation <= 1.0);
	assert(p.nSeedCells >= 1);
	assert(p.nTotalCells >= p.nSeedCells);
	assert(p.fWinners > 0.0 && p.fWinners <= 1.0);
	if(p.exponentialGrowthActive) {
		assert(p.nDivisionsPerRound > 0);
		assert(p.nDivisionsPerRound < 20);
	}
	assert(p.pPlasmaExport >= 0.0 && p.pPlasmaExport <= 1.0);
	assert(p.pMemoryExportFirstRound >= 0.0 && p.pMemoryExportFirstRound <= 1.0);
	assert(p.pMemoryExportLastRound >= 0.0 && p.pMemoryExportLastRound <= 1.0);
}

void verify(PlasmaParams & p)
{
	assert(p.discardLogC < p.initialLogC);
	assert(p.decayRate > 0.0 && p.decayRate <= 1.0);
}
