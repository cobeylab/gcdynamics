//
//  Simulation.cpp
//  bcellmodel
//
//  Created by Ed Baskerville on 1/27/14.
//  Copyright (c) 2014 Cobey Lab. All rights reserved.
//

#include "Simulation.h"
#include "Host.h"
#include "Antigen.h"
#include <memory>
#include "shared.h"

using namespace std;
using namespace zppdb;
using namespace zppsim;

Simulation::Simulation(SimParameters * p) : p(p), db(p->dbFilename), dbm(db, p->dbTablesEnabled), nLoci(calculateNLoci(*p))
{
	// Generate random seed if set to 0
	if(!p->randomSeed.present() || p->randomSeed == 0) {
		random_device rd;
		uniform_int_distribution<uint32_t> ud(0, std::numeric_limits<uint32_t>::max());
		p->randomSeed = ud(rd);
	}
	
	// Write params to meta table
	db.beginTransaction();
	MetaRow row;
	row.key = "parameters";
	row.value = p->toJsonString();
	db.insert(dbm.metaTable, row);
	db.commitWithRetry(0.1, 100, std::cerr);
	
	rng = rng_t(int64_t(p->randomSeed));
}

int64_t Simulation::calculateNLoci(SimParameters & p)
{
	int64_t nLoci = 0;
	for(int64_t i = 0; i < p.sequences.chains.size(); i++) {
		nLoci += p.sequences.chains[i].truncatedLength;
	}
	return nLoci;
}

void Simulation::run()
{
	// Create a single host (seeded with random naive B cells)
	db.beginTransaction();
	unique_ptr<Host> host(new Host(*this));
	db.commitWithRetry(0.1, 100, std::cerr);
	
	// Perform nInfections rounds of infection
	unique_ptr<Antigen> antigen;
	for(int64_t infectionNumber = 0; infectionNumber < p->infections.size(); infectionNumber++) {
		cerr << "infections.size(): " << p->infections.size() << endl;
		
		db.beginTransaction();
		fprintf(stderr, "Starting infection %lld...\n", infectionNumber);
		clock_t clockStart = clock();
		
		// Generate infection: random initial state or modification from previous infection
		if(infectionNumber == 0) {
			antigen = unique_ptr<Antigen>(new Antigen(
				*this, *p, 0,
				rng, &db, &dbm.antigenNeighborsTable, &dbm.antigenEnergiesTable
			));
		}
		else {
			antigen->mutate(rng);
		}
		db.commitWithRetry(0.1, 100, std::cerr);
		
		// Process infection by host, writing cells out to DB in the process
		host->runInfection(*antigen, rng);
		
		clock_t clockEnd = clock();
		fprintf(stderr, "Infection %lld complete (elapsed time: %.3f s).\n", infectionNumber, elapsed(clockStart, clockEnd));
	}
	
	// Calculate test statistics
	db.beginTransaction();
	
	if(db.tableExists(dbm.testStatsTable)) {
		TestStatsRow row;
		
		row.name = "mean_max_affinity_memory";
		row.value = host->meanMaxAffinityMemory(*antigen, rng);
		db.insert(dbm.testStatsTable, row);
	}
	
	db.commitWithRetry(0.1, 100, std::cerr);
}
