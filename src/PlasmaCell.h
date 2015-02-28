#ifndef bcellmodel_PlasmaCell_h
#define bcellmodel_PlasmaCell_h

#include "BCell.h"
#include "Antigen.h"

class PlasmaCell : public BCell
{
public:
	PlasmaCell(uint32_t id, PlasmaCell const & parent)
		: BCell(id, parent),
		logCShortLived(parent.logCShortLived),
		logCLongLived(parent.logCLongLived)
	{
	}
	
	PlasmaCell(uint32_t id, BCell & parent, Antigen & antigen, double logCShortLived, double logCLongLived)
		: BCell(id, parent), logCShortLived(logCShortLived), logCLongLived(logCLongLived)
	{
	}
	
	double getLogShortLivedConcentration() const
	{
		return logCShortLived;
	}
	
	void setLogShortLivedConcentration(double logCShortLived)
	{
		this->logCShortLived = logCShortLived;
	}
	
	double getLogLongLivedConcentration() const
	{
		return logCLongLived;
	}
	
	void setLogLongLivedConcentration(double logCLongLived)
	{
		this->logCLongLived = logCLongLived;
	}
	
	double getLogConcentration() const
	{
		// log-sum-exp(logCLongLived, logCShortLived)
		double z = std::max(logCLongLived, logCShortLived);
		return z + log(exp(logCLongLived - z) + exp(logCShortLived - z));
	}
	
private:
	double logCShortLived;
	double logCLongLived;
};


#endif
