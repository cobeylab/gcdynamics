/**
 Header for class representing B cells.
*/

#ifndef __bcellmodel__BCell__
#define __bcellmodel__BCell__

#include <memory>
#include <vector>
#include <unordered_map>
#include "SimParameters.h"
#include "Sequence.h"

#define NO_PARENT std::numeric_limits<uint32_t>::max()

class BCell
{
friend class Epitope;
public:
	BCell(uint32_t id, Sequence const & seq, SequenceCache & cache) :
		id(id),
		parentId(NO_PARENT),
		seqPtr(cache.get(seq))
	{
	}
	
	BCell(uint32_t id, BCell const & parent) :
		id(id),
		parentId(parent.id),
		seqPtr(parent.seqPtr)
	{
	}
	
	BCell(uint32_t id, uint8_t alphabetSize, uint16_t nLoci,
		zppsim::rng_t & rng, SequenceCache & cache
	) :
		id(id),
		parentId(NO_PARENT),
		seqPtr(cache.get(Sequence(alphabetSize, nLoci, rng)))
	{
	}
	
	BCell(
		uint32_t id, BCell const & parent,
		zppsim::rng_t & rng, double pMutation,
		SequenceCache & cache
	) :
		id(id), parentId(parent.id),
		seqPtr(cache.get(Sequence(*(parent.seqPtr), rng, pMutation)))
	{
	}
	
	uint8_t get(uint32_t const locus) const
	{
		return seqPtr->get(locus);
	}
	
	std::string toString(std::string const & alphabet) const
	{
		return seqPtr->toString(alphabet);
	}
	
	uint32_t getId() const
	{
		return id;
	}
	
	uint32_t getParentId() const
	{
		return parentId;
	}
	
	uint32_t getSequenceId() const
	{
		return seqPtr->id;
	}

protected:
	uint32_t id;
	uint32_t parentId;
	std::shared_ptr<Sequence> seqPtr;
};

struct BCellDoublePair
{
	BCellDoublePair(BCell * cellPtr, double value) : cellPtr(cellPtr), value(value)
	{
	}
	
	BCell * cellPtr;
	double value;
	
	bool operator>(BCellDoublePair other) const
	{
		return value > other.value;
	}
};

#endif /* defined(__bcellmodel__BCell__) */
