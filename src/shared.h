/**
	\file shared.h
	\brief Shared types, macros, utility functions.
*/

#ifndef bcellmodel_shared_h
#define bcellmodel_shared_h

#include <cstdio>
#include <random>
#include <vector>
#include <functional>

float elapsed(clock_t clockStart, clock_t clockEnd);

double meanFromFunction(std::function<double(uint32_t)> func, uint32_t size);
double standardDeviationSample(double sum, double sumSq, uint32_t size);
double standardDeviationSample(std::function<double(uint32_t)> func, uint32_t size);
double standardDeviationSample(double mean, std::function<double(uint32_t)> func, uint32_t size);
double standardDeviationPopulation(double sum, double sumSq, uint32_t size);
double standardDeviationPopulation(std::function<double(uint32_t)> func, uint32_t size);
double standardDeviationPopulation(double mean, std::function<double(uint32_t)> func, uint32_t size);

#endif
