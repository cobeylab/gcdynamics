/**
	Implementation of shared utility functions.
*/

#include "shared.h"
#include <iostream>
#include <sstream>
#include <cassert>

using namespace std;

static uint64_t nextId = 0;

uint64_t getNextId()
{
	return nextId++;
}

float elapsed(clock_t clockStart, clock_t clockEnd)
{
	return float(clockEnd - clockStart) / CLOCKS_PER_SEC;
}

double meanFromFunction(std::function<double(uint32_t)> func, uint32_t size)
{
	assert(size > 0);
	
	double sum = 0.0;
	for(uint32_t i = 0; i < size; i++) {
		sum += func(i);
	}
	return sum / size;
}

double standardDeviationSample(double sum, double sumSq, uint32_t n)
{
	assert(n > 1);
	
	double mean = sum / n;
	double meanSq = sumSq / n;
	double var = double(n) / double(n-1) * (meanSq - mean*mean);
	
	return sqrt(var);
}

double standardDeviationSample(std::function<double(uint32_t index)> func, uint32_t size)
{
	assert(size > 1);
	
	double sum = 0.0;
	double sumSq = 0.0;
	
	for(uint32_t i = 0; i < size; i++) {
		double value = func(i);
		sum += value;
		sumSq +=  value * value;
	}
	
	return standardDeviationSample(sum, sumSq, size);
}

double standardDeviationSample(double mean, std::function<double(uint32_t)> func, uint32_t size)
{
	assert(size > 1);
	
	double sumSqDiff = 0.0;
	
	for(uint32_t i = 0; i < size; i++) {
		double diff = mean - func(i);
		sumSqDiff += diff * diff;
	}
	
	return sqrt(sumSqDiff / (size - 1));
}

double standardDeviationPopulation(double sum, double sumSq, uint32_t n)
{
	assert(n >= 1);
	
	double mean = sum / n;
	double meanSq = sumSq / n;
	double var = max(0.0, meanSq - mean*mean);
	
	return sqrt(var);
}

double standardDeviationPopulation(std::function<double(uint32_t index)> func, uint32_t size)
{
	assert(size >= 1);
	
	double sum = 0.0;
	double sumSq = 0.0;
	
	for(uint32_t i = 0; i < size; i++) {
		double value = func(i);
		sum += value;
		sumSq +=  value * value;
	}
	
	return standardDeviationPopulation(sum, sumSq, size);
}

double standardDeviationPopulation(double mean, std::function<double(uint32_t)> func, uint32_t size)
{
	assert(size >= 1);
	
	double sumSqDiff = 0.0;
	
	for(uint32_t i = 0; i < size; i++) {
		double diff = mean - func(i);
		sumSqDiff += diff * diff;
	}
	
	return sqrt(sumSqDiff / size);
}
