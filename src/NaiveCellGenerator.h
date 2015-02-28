//
//  NaiveCellGenerator.h
//  bcellmodel
//
//  Created by Ed Baskerville on 3/25/14.
//  Copyright (c) 2014 Cobey Lab. All rights reserved.
//

#ifndef __bcellmodel__NaiveCellGenerator__
#define __bcellmodel__NaiveCellGenerator__

#include <iostream>
#include "SimParameters.h"
#include "BCell.h"
#include "DatabaseManager.h"

class Segment
{
friend class NaiveCellGenerator;
public:
	virtual Sequence generateSequence(zppsim::rng_t & rng) = 0;
	virtual int64_t getSegmentType() = 0;
};

class AlleleSegment : public Segment
{
public:
friend class NaiveCellGenerator;
	AlleleSegment(
		zppsim::rng_t & rng, int64_t alphabetSize, uint16_t minLength,
		std::vector<double> const & frequencyDistVec, std::vector<double> const & lengthDistVec
	):
		alleleDist(frequencyDistVec.begin(), frequencyDistVec.end())
	{
		size_t nAlleles = frequencyDistVec.size();
//		std::cerr << "Creating " << nAlleles << " alleles for segment" << '\n';
		std::discrete_distribution<uint16_t> lengthDist(lengthDistVec.begin(), lengthDistVec.end());
		for(size_t i = 0; i < nAlleles; i++) {
			uint16_t length = minLength + lengthDist(rng);
			alleles.push_back(Sequence(alphabetSize, length, rng));
//			std::cerr << "allele " << i << ": " << alleles[i].toString(p->alphabet) << '\n';
		}
	}
	
	virtual Sequence generateSequence(zppsim::rng_t & rng)
	{
		return alleles[alleleDist(rng)];
	}
	
	virtual int64_t getSegmentType()
	{
		return SegmentType::ALLELE();
	}

private:
	std::discrete_distribution<size_t> alleleDist;
	std::vector<Sequence> alleles;
};

class RandomInsertSegment : public Segment
{
friend class NaiveCellGenerator;
public:
	RandomInsertSegment(int64_t alphabetSize, uint16_t minLength, std::vector<double> const & lengthDistVec):
		alphabetSize(alphabetSize), minLength(minLength), lengthDist(lengthDistVec.begin(), lengthDistVec.end())
	{
//		std::cerr << "Creating random insert segment with maximum length " << lengthDistVec.size() << '\n';
//		std::cerr << "lengthDistVec : " << std::endl;
//		for(double lengthDist : lengthDistVec) {
//			std::cerr << lengthDist << ", ";
//		}
//		std::cerr << std::endl;
	}
	
	virtual Sequence generateSequence(zppsim::rng_t & rng)
	{
		return Sequence(alphabetSize, minLength + lengthDist(rng), rng);
	}
	
	virtual int64_t getSegmentType()
	{
		return SegmentType::RANDOM_INSERT();
	}
	
private:
	uint8_t alphabetSize;
	uint16_t minLength;
	std::discrete_distribution<uint16_t> lengthDist;
};

class Chain
{
friend class NaiveCellGenerator;
public:
	Chain(zppsim::rng_t & rng, int64_t alphabetSize, ChainParams & chainInfo):
		truncatedLength(uint16_t(chainInfo.truncatedLength))
	{
		for(size_t i = 0; i < chainInfo.segments.size(); i++) {
			SegmentParams segInfo = chainInfo.segments[i];
//			std::cerr << "Creating segment " << i << '\n';
			if(segInfo.type == SegmentType::ALLELE()) {
				segments.push_back(std::unique_ptr<Segment>(
					new AlleleSegment(
						rng, alphabetSize, segInfo.minLength,
						segInfo.alleleDistribution.toDoubleVector(),
						segInfo.lengthDistribution.toDoubleVector()
					)
				));
			}
			else {
				uint16_t minLength = segInfo.minLength.present()
					? uint16_t(segInfo.minLength)
					: 0;
				
				segments.push_back(std::unique_ptr<Segment>(
					new RandomInsertSegment(alphabetSize, minLength, segInfo.lengthDistribution.toDoubleVector())
				));
			}
		}
	}
	
	Sequence generateSequence(zppsim::rng_t & rng)
	{
		std::vector<Sequence> sequences;
		for(size_t i = 0; i < segments.size(); i++) {
//			std::cerr << "Generating sequence for segment " << i << ":" << '\n';
			sequences.push_back(segments[i]->generateSequence(rng));
//			std::cerr << sequences[i].toString(p->alphabet) << '\n';
		}
		return Sequence(Sequence::concatenate(sequences), truncatedLength);
	}
private:
	uint16_t truncatedLength;
	std::vector<std::unique_ptr<Segment>> segments;
};

class NaiveCellGenerator
{
public:
	NaiveCellGenerator(SimParameters & p, zppdb::Database & db, zppdb::Table<NaiveAlleleRow> & table, zppsim::rng_t & rng)
	{
		for(size_t i = 0; i < p.sequences.chains.size(); i++) {
//			std::cerr << "Creating chain " << i << '\n';
			chains.emplace_back(rng, p.sequences.alphabetSize, p.sequences.chains[i]);
			
			Chain * chainPtr = &(chains[i]);
			for(size_t j = 0; j < chainPtr->segments.size(); j++) {
				Segment * segPtr = chainPtr->segments[j].get();
				if(segPtr->getSegmentType() == SegmentType::ALLELE()) {
					AlleleSegment * aSegPtr = (AlleleSegment *)segPtr;
					for(size_t k = 0; k < aSegPtr->alleles.size(); k++) {
						Sequence * seq = &aSegPtr->alleles[k];
						
						if(db.tableExists(table)) {
							NaiveAlleleRow row;
	//						std::cerr << i << ", " << j << ", " << k << '\n';
							row.chain_id = i;
							row.segment_id = j;
							row.allele_id = k;
							row.sequence = seq->toString(p.sequences.alphabet);
							db.insert(table, row);
						}
					}
				}
			}
		}
	}
	
	BCell * generateCell(uint32_t cellId, zppsim::rng_t & rng, SequenceCache & cache)
	{
		std::vector<Sequence> sequences;
		for(size_t i = 0; i < chains.size(); i++) {
//			std::cerr << "generating sequence for chain " << i << '\n';
			sequences.push_back(chains[i].generateSequence(rng));
//			std::cerr << "chain sequence: " << sequences[i].toString(p->alphabet) << '\n';
		}
		return new BCell(cellId, Sequence::concatenate(sequences), cache);
	}

private:
	std::vector<Chain> chains;
};

#endif
