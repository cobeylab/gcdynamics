//
//  Antigen.cpp
//  bcellmodel
//
//  Created by Ed Baskerville on 10/18/13.
//  Copyright (c) 2013 Cobey Lab. All rights reserved.
//

#include "Antigen.h"
#include "Epitope.h"
#include "DatabaseManager.h"

#include "BCell.h"
#include "Simulation.h"
#include "SimParameters.h"

using namespace std;

Antigen::Antigen(
	Simulation & sim, SimParameters & params, int64_t id,
	zppsim::rng_t & rng,
	zppdb::Database * dbPtr, zppdb::Table<AntigenNeighborRow> * neighborTablePtr,
	zppdb::Table<AntigenEnergyRow> * energyTablePtr
) :
	id(id),
	epitopeRows(params.epitopes.rows),
	epitopeCols(params.epitopes.cols),
	nLoci(sim.nLoci),
	affinityA(params.epitopes.affinityA),
	affinityB(params.epitopes.affinityB),
	quantity(params.infections[id].antigenQuantity),
	dbPtr(dbPtr),
	neighborTablePtr(neighborTablePtr),
	energyTablePtr(energyTablePtr),
	epitopes(epitopeRows, vector<Epitope *>(epitopeCols, NULL))
{
	uniform_int_distribution<int64_t> energySeedGenerator;
	
	for(int64_t i = 0; i < epitopeRows; i++) {
		for(int64_t j = 0; j < epitopeCols; j++) {
			double energyMean = params.epitopes.energyMean[i][j];
			Epitope * epitopePtr = new Epitope(
				params,
				this,
				uint32_t(params.epitopes.neighborSeed[i][j]),
				uint32_t(params.epitopes.energySeed[i][j]),
				rng, i, j, energyMean,
				dbPtr, neighborTablePtr, energyTablePtr
			);
			
			epitopesVec.push_back(unique_ptr<Epitope>(epitopePtr));
			epitopes[i][j] = epitopePtr;
		}
	}
}

double Antigen::getLogAffinity(
	BCell const & bCell, int64_t i, int64_t j,
	zppsim::rng_t & rng
) const {
	return affinityA - affinityB * epitopes[i][j]->getEnergy(bCell);
}

double Antigen::getMaxLogAffinity(
	BCell const & bCell, int64_t & iOut, int64_t & jOut,
	zppsim::rng_t & rng
) const {
	double maxLogAffinity = -std::numeric_limits<double>::infinity();
	for(int64_t i = 0; i < epitopeRows; i++) {
		for(int64_t j = 0; j < epitopeCols; j++) {
			double logAffinity = getLogAffinity(bCell, i, j, rng);
			if(logAffinity > maxLogAffinity) {
				maxLogAffinity = logAffinity;
				iOut = (int64_t)i;
				jOut = (int64_t)j;
			}
		}
	}
	return maxLogAffinity;
}

/*double Antigen::getLogAffinity(PlasmaCell const & pCell, int64_t i, int64_t j) const
{
	return epitopes[i][j]->getLogAffinity(pCell);
}*/

double Antigen::getEnergy(BCell const & bCell, int64_t i, int64_t j, zppsim::rng_t & rng) const
{
	return epitopes[i][j]->getEnergy(bCell);
}

vector<vector<double>> Antigen::getEnergies(BCell const & bCell, zppsim::rng_t & rng) const
{
	vector<vector<double>> energies;
	for(int64_t i = 0; i < epitopeRows; i++) {
		energies.push_back(vector<double>(epitopeCols));
		for(int64_t j = 0; j < epitopeCols; j++) {
			energies[i][j] = epitopes[i][j]->getEnergy(bCell);
		}
	}
	return energies;
}

void Antigen::mutate(zppsim::rng_t & rng)
{
	id++;
	for(int64_t i = 0; i < epitopeRows; i++) {
		for(int64_t j = 0; j < epitopeCols; j++) {
			epitopes[i][j]->mutate(rng);
		}
	}
}
