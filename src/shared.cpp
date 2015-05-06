/**
	Implementation of shared utility functions.
*/

#include "shared.h"
#include <iostream>
#include <sstream>
#include <cassert>
#include <openssl/sha.h>

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

std::vector<uint8_t> hostToNetworkBytes(uint32_t x)
{
	std::vector<uint8_t> byte_vec(4, 0);
	byte_vec[0] = (x & 0xFF000000U) >> 24;
	byte_vec[1] = (x & 0x00FF0000U) >> 16;
	byte_vec[2] = (x & 0x0000FF00U) >> 8;
	byte_vec[3] = (x & 0x000000FFU);
	return byte_vec;
}

uint32_t networkBytesToHostUInt32(uint8_t bytes[4])
{
	return (uint32_t(bytes[0]) << 24)
		+ (uint32_t(bytes[1]) << 16)
		+ (uint32_t(bytes[2]) << 8)
		+ uint32_t(bytes[3]);
}

double sha1Normal(std::vector<uint8_t> bytes)
{
	double uint32_max = std::numeric_limits<uint32_t>::max();
	uint8_t hash[SHA_DIGEST_LENGTH];
	SHA_CTX sha1;
	
	SHA1_Init(&sha1);
	SHA1_Update(&sha1, (void *)bytes.data(), bytes.size());
	SHA1_Final(hash, &sha1);
	
	double u = networkBytesToHostUInt32(hash) / uint32_max;
	double v = networkBytesToHostUInt32(hash + 4) / uint32_max;
	
	return sqrt(-2.0 * log(u)) * sin(2.0 * M_PI * v);
}
