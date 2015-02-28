//
//  Simulation.h
//  bcellmodel
//
//  Created by Ed Baskerville on 1/27/14.
//  Copyright (c) 2014 Cobey Lab. All rights reserved.
//

#ifndef __bcellmodel__Simulation__
#define __bcellmodel__Simulation__

#include "SimParameters.h"
#include "DatabaseManager.h"
#include "zppsim_random.hpp"

class Simulation
{
public:
	SimParameters * p;
	zppdb::Database db;
	DatabaseManager dbm;
	zppsim::rng_t rng;
	
	Simulation(SimParameters * p);
	void run();
	
	int64_t const nLoci;
private:
	int64_t calculateNLoci(SimParameters & params);
};

#endif /* defined(__bcellmodel__Simulation__) */
