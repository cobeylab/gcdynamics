#ifndef __bcellmodel__Sequence__
#define __bcellmodel__Sequence__

#include <memory>
#include <vector>
#include "zppsim_random.hpp"
#include <cassert>
#include <memory>
#include <random>
#include <unordered_map>
#include <sstream>
#include <iostream>
#include "DatabaseManager.h"

#define N_BITS_PER_LOCUS(alphabetSize) \
	(alphabetSize >> 6) ? 0 : \
	(alphabetSize >> 5) ? 6 : \
	(alphabetSize >> 4) ? 5 : \
	(alphabetSize >> 3) ? 4 : \
	(alphabetSize >> 2) ? 3 : \
	(alphabetSize >> 1) ? 2 : \
	(alphabetSize > 1) ? 1 : \
	0

class Sequence
{
public:
	static Sequence concatenate(std::vector<Sequence> const & sequences)
	{
		uint8_t alphabetSize = sequences[0].alphabetSize;
		uint16_t nLoci = 0;
		std::vector<uint8_t> catSeqVec;
		for(Sequence const & seq : sequences) {
			assert(seq.alphabetSize == alphabetSize);
			nLoci += seq.nLoci;
			std::vector<uint8_t> seqVec = seq.toVector();
			catSeqVec.insert(catSeqVec.end(), seqVec.begin(), seqVec.end());
		}
		assert(nLoci == catSeqVec.size());
		
		return Sequence(alphabetSize, catSeqVec);
	}
	
	Sequence(Sequence const & seq, uint32_t id) :
		alphabetSize(seq.alphabetSize),
		nLoci(seq.nLoci),
		nBitsPerLocus(seq.nBitsPerLocus),
		id(id),
		bits(seq.bits),
		hashValue(seq.hashValue)
	{
	}
	
	Sequence(Sequence const & seq, uint16_t nLociTrimmed):
		alphabetSize(seq.alphabetSize),
		nLoci(nLociTrimmed),
		nBitsPerLocus(seq.nBitsPerLocus),
		id(std::numeric_limits<uint32_t>::max()),
		bits(nLoci*nBitsPerLocus)
	{
		assert(nLoci <= seq.nLoci);
		for(uint16_t i = 0; i < nLoci; i++) {
			set(bits, i, seq.get(i));
		}
		cacheHash(bits);
	}
	
	Sequence(uint8_t alphabetSize, uint16_t nLoci) :
		alphabetSize(alphabetSize), nLoci(nLoci),
		nBitsPerLocus(N_BITS_PER_LOCUS(alphabetSize)),
		id(std::numeric_limits<uint32_t>::max()),
		bits(nLoci*nBitsPerLocus)
	{
		assert(nBitsPerLocus != 0);
		std::uniform_int_distribution<uint32_t> unif(0, alphabetSize - 1);
		for(uint32_t i = 0; i < nLoci; i++) {
			set(bits, i, 0);
		}
		cacheHash(bits);
	}
	
	Sequence(uint8_t alphabetSize, uint16_t nLoci, zppsim::rng_t & rng) :
		alphabetSize(alphabetSize), nLoci(nLoci),
		nBitsPerLocus(N_BITS_PER_LOCUS(alphabetSize)),
		id(std::numeric_limits<uint32_t>::max()),
		bits(nLoci*nBitsPerLocus)
	{
		assert(nBitsPerLocus != 0);
		std::uniform_int_distribution<uint32_t> unif(0, alphabetSize - 1);
		for(uint32_t i = 0; i < nLoci; i++) {
			set(bits, i, unif(rng));
		}
		cacheHash(bits);
	}
	
	Sequence(Sequence const & seq, zppsim::rng_t & rng, double pMutation) :
		alphabetSize(seq.alphabetSize), nLoci(seq.nLoci), nBitsPerLocus(seq.nBitsPerLocus),
		id(std::numeric_limits<uint32_t>::max()),
		bits(seq.bits)
	{
		assert(nBitsPerLocus != 0);
		
		// Generate random loci to be mutated and set each locus to a random new value
		std::vector<uint32_t> indexes = zppsim::drawMultipleBernoulli(rng, uint32_t(nLoci), pMutation);
		for(uint32_t index : indexes) {
			set(bits, index, zppsim::drawUniformIndexExcept(rng, alphabetSize, get(bits, index)));
		}
		cacheHash(bits);
	}
	
	Sequence(uint8_t alphabetSize, std::vector<uint8_t> seq) :
		alphabetSize(alphabetSize), nLoci(uint16_t(seq.size())),
		nBitsPerLocus(N_BITS_PER_LOCUS(alphabetSize)),
		id(std::numeric_limits<uint32_t>::max()),
		bits(nLoci*nBitsPerLocus)
	{
		assert(nBitsPerLocus != 0);
		for(uint32_t i = 0; i < nLoci; i++) {
			set(bits, i, seq[i]);
		}
		cacheHash(bits);
	}
	
	Sequence(
		uint8_t alphabetSize,
		std::string const & seqStr,
		std::unordered_map<char, uint8_t> & alphabetMap
	) :
		alphabetSize(alphabetSize), nLoci(uint16_t(seqStr.size())),
		nBitsPerLocus(N_BITS_PER_LOCUS(alphabetSize)),
		id(std::numeric_limits<uint32_t>::max()),
		bits(nLoci*nBitsPerLocus)
	{
		assert(nBitsPerLocus != 0);
		for(uint32_t i = 0; i < nLoci; i++) {
			set(bits, i, alphabetMap[seqStr[i]]);
		}
		cacheHash(bits);
	}
	
	uint8_t get(uint32_t const index) const
	{
		return get(bits, index);
	}
	
	std::vector<uint8_t> toVector() const
	{
		std::vector<uint8_t> vec;
		vec.reserve(nLoci);
		for(uint16_t i = 0; i < nLoci; i++) {
			vec.push_back(get(i));
		}
		return vec;
	}
	
	std::string toString(std::string const & alphabet) const
	{
		assert(alphabetSize < alphabet.size());
		
		std::string str(nLoci, alphabet[0]);
		for(uint32_t i = 0; i < nLoci; i++) {
			str[i] = alphabet[get(i)];
		}
		return str;
	}
	
	size_t hash() const
	{
		return hashValue;
	}
	
	bool equal_to(Sequence const & x) const
	{
		return alphabetSize == x.alphabetSize && nLoci == x.nLoci && bits == x.bits;
	}
	
	uint8_t const alphabetSize;
	uint16_t const nLoci;
	uint8_t const nBitsPerLocus;
	
	uint32_t const id;
	
private:
	std::vector<bool> bits;
	size_t hashValue;
	
	void cacheHash(std::vector<bool> & bits)
	{
		std::hash<std::vector<bool>> hasher;
		hashValue = hasher(bits);
	}
	
	void set(std::vector<bool> & bits, uint32_t const index, uint8_t value)
	{
		assert(value < alphabetSize);
		for(uint8_t i = 0; i < nBitsPerLocus; i++) {
			bits[index*nBitsPerLocus + i] = value & (1 << i);
		}
	}
	
	uint8_t get(std::vector<bool> const & bits, uint32_t const index) const
	{
		uint8_t value = 0;
		for(uint8_t i = 0; i < nBitsPerLocus; i++) {
			value |= bits[index*nBitsPerLocus + i] ? (1 << i) : 0;
		}
		return value;
	}
};

typedef std::shared_ptr<Sequence> SequencePtr;

struct HashSequence
{
	size_t operator()(Sequence const & seq) const
	{
		return seq.hash();
	}
};

bool operator==(Sequence const & s1, Sequence const & s2);

class SequenceCache
{
public:
	SequenceCache() :
		nextId(0),
		nextIdToWrite(0)
	{
	}
	
	SequencePtr get(Sequence const & seq)
	{
		auto itr = cache.find(seq);
		if(itr == cache.end()) {
			assert(nextId < std::numeric_limits<uint32_t>::max());
			cache[seq] = std::make_shared<Sequence>(seq, nextId++);
			idSeqMap[cache[seq]->id] = cache[seq];
		}
		return cache[seq];
	}
	
	void writeToDatabase(Database & db, Table<SequenceRow> & table, std::string const & alphabet)
	{
		if(!db.tableExists(table)) {
			return;
		}
		
		std::vector<uint32_t> idsToWrite;
		for(auto itr = idSeqMap.begin(); itr != idSeqMap.end(); ++itr) {
			if(itr->first >= nextIdToWrite) {
				idsToWrite.push_back(itr->first);
			}
		}
		std::sort(idsToWrite.begin(), idsToWrite.end());
		for(uint32_t id : idsToWrite) {
			SequencePtr seqPtr = idSeqMap[id];
			SequenceRow row;
			row.sequence_id = int64_t(id);
			row.sequence = seqPtr->toString(alphabet);
			db.insert(table, row);
		}
		if(!idsToWrite.empty()) {
			nextIdToWrite = idsToWrite.back() + 1;
		}
	}
	
	void collectGarbage()
	{
		for(auto itr = cache.begin(); itr != cache.end(); ++itr) {
			if(itr->second.use_count() == 2) {
				idSeqMap.erase(itr->second->id);
				itr = cache.erase(itr);
				if(itr == cache.end()) {
					break;
				}
			}
		}
	}
private:
	std::unordered_map<uint32_t, SequencePtr> idSeqMap;
	std::unordered_map<Sequence, SequencePtr, HashSequence> cache;
	uint32_t nextId;
	uint32_t nextIdToWrite;
};

#endif /* defined(__bcellmodel__Sequence__) */
