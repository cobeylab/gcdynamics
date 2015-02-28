/**
	\mainpage bcellmodel
	
	Code documentation for B-cell affinity maturation model.
	
	Useful sections:
	- SimParameters_jsontypes.h: parameter documentation
	- DatabaseManager_dbtypes.h: database table documentation
*/

/**
	\file main.cpp
	\brief Globals and main function.
*/

#include <iostream>
#include <fstream>
#include <memory>
#include <sys/resource.h>
#include <ctime>
#include <cinttypes>
#include <sqlite3.h>
#include <sstream>

#include "SimParameters.h"
#include "Simulation.h"
#include "shared.h"

using namespace std;

/**
	\fn main
	\brief Loads parameters and runs simulation.
*/
int main(int argc, char const * argv[])
{
	try {
		// Start timer
		clock_t startClock = clock();
		time_t startTime = time(nullptr);
		fprintf(stderr, "Starting at %s", ctime(&startTime));
		
		// Load parameters
		char const * paramsFilename;
		if(argc > 1) {
			paramsFilename = argv[1];
		}
		else {
			paramsFilename = "parameters.json";
		}
		ifstream paramsStream(paramsFilename);
		if(!paramsStream.good()) {
			stringstream errStream;
			errStream << "Parameter file " << paramsFilename << " not present.";
			throw runtime_error(errStream.str());
		}
		
		JsonElement paramsJson = parseJson(paramsStream);
		paramsStream.close();
		
		SimParameters params;
		params.loadJson(paramsJson);
		
		// Create and run simulation
		unique_ptr<Simulation> simPtr = unique_ptr<Simulation>(new Simulation(&params));
		simPtr->run();
		
		// Calculate total time
		time_t endTime = time(nullptr);
		clock_t endClock = clock();
		fprintf(stderr, "Ending at %s", ctime(&endTime));
		fprintf(stderr, "Total elapsed time: %f\n", elapsed(startClock, endClock));
		
		// Output memory usage (in bytes)
		rusage resourceUsage;
		getrusage(RUSAGE_SELF, &resourceUsage);
		fprintf(stderr, "Memory usage: %ld\n", resourceUsage.ru_maxrss);
	}
	catch(std::runtime_error e) {
		fprintf(stderr, "Exception thrown.\n");
		cerr << e.what() << endl;
	}
	
	return 0;
}



