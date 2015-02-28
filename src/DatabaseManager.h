//
//  DatabaseManager.h
//  bcellmodel
//
//  Created by Ed Baskerville on 1/29/14.
//  Copyright (c) 2014 Cobey Lab. All rights reserved.
//

#ifndef __bcellmodel__DatabaseManager__
#define __bcellmodel__DatabaseManager__

#include <iostream>
#include <unordered_map>
#include <map>
#include <memory>
#include "zppdb.hpp"
#include "SimParameters.h"

using namespace zppdb;

#ifndef DOXYGEN_SHOULD_SKIP_THIS

/**
	\brief Row type for meta table.
*/
ZPPDB_DEFINE_ROW_TYPE(
	MetaRow,
	
	/**
		\brief meta-information key (`"parameters"`)
	*/
	((Text)(key))
	
	/**
		\brief value (e.g., JSON parameter text)
	*/
	((Text)(value))
)

/**
	\brief Row type for test statistics table.
*/
ZPPDB_DEFINE_ROW_TYPE(
	TestStatsRow,
	
	/**
		\brief name of statistic
	*/
	((Text)(name))

	/**
		\brief value of statistic
	*/
	((Real)(value))
)

/**
	\brief Row type for gcStats table.
*/
ZPPDB_DEFINE_ROW_TYPE(
	GCStatsRow,
	((Integer)(infection_id))
	((Integer)(gc_round))
	((Integer)(n_active_gcs))
)

/**
	\brief Sequence of alleles used to construct naive cell sequences.
*/
ZPPDB_DEFINE_ROW_TYPE(
	NaiveAlleleRow,
	
	/**
		\brief zero-indexed chain ID
	*/
	((Integer)(chain_id))
	
	/**
		\brief zero-indexed segment ID
	*/
	((Integer)(segment_id))
	
	/**
		\brief zero-indexed allele ID
	*/
	((Integer)(allele_id))
	
	/**
		\brief amino acid sequence
	*/
	((Text)(sequence))
)

/**
	\brief Row for bcells table, containing full complement of all cells that appear in the simulation.
	
	Each cell should appear exactly once in this table, not necessarily in
	ascending order of ID.
*/
ZPPDB_DEFINE_ROW_TYPE(
	BCellRow,
	
	/**
		\brief unique ID of cell
	*/
	((Integer)(cell_id))
	
	/**
		\brief unique ID of parent cell
	*/
	((Integer)(parent_id))
	
	/**
		\brief cell type (`"naive"`, `"gcfd"`, `"memory"`, `"plasma"`)
	*/
	((Text)(type))
	
	/**
		\brief ID of antigen-recognizing sequence
	*/
	((Integer)(sequence_id))
	
	/**
		\brief ID of infection in which cell was created
	*/
	((Integer)(infection_id))
	
	/**
		\brief ID of GC in which cell was created
	*/
	((Integer)(gc_id))
	
	/**
		\brief GC round during which cell was created
	*/
	((Integer)(gc_round))
)

/**
	\brief Row for sequence table, containing sequence IDs that are shared among cells
*/
ZPPDB_DEFINE_ROW_TYPE(
	SequenceRow,
	
	/**
		\brief Sequence ID
	*/
	((Integer)(sequence_id))

	/**
		\brief Sequence
	*/
	((Text)(sequence))
)

/**
	\brief Row for gcfdCells table (full population of cells in all GCs in every round).
	
	The population is written out after proliferation and selection.
*/
ZPPDB_DEFINE_ROW_TYPE(
	GCFDCellRow,
	
	/**
		\brief ID of active infection
	*/
	((Integer)(infection_id))
	
	/**
		\brief ID of germinal center
	*/
	((Integer)(gc_id))
	
	/**
		\brief GC round
	*/
	((Integer)(gc_round))
	
	/**
		\brief ID of cell
	*/
	((Integer)(cell_id))
)

/**
	\brief Row for gcSeedCandidates table (all cells that were chosen as potential seeds for GCs).
*/
ZPPDB_DEFINE_ROW_TYPE(
	GCSeedCandidateRow,
	
	/**
		\brief ID of infection
	*/
	((Integer)(infection_id))
	
	/**
		\brief cell ID
	*/
	((Integer)(cell_id))
)

/**
	\brief Cells that were chosen as seeds for GCs.
*/
ZPPDB_DEFINE_ROW_TYPE(
	GCSeedRow,
	
	/**
		\brief infection ID
	*/
	((Integer)(infection_id))
	
	/**
		\brief GC ID
	*/
	((Integer)(gc_id))
	
	/**
		\brief Cell ID
	*/
	((Integer)(cell_id))
)

/**
	\brief Row for affinity table (intrinsic affinity of cells to all epitopes, all infections).
	
	Naive cells are not included.
*/
ZPPDB_DEFINE_ROW_TYPE(
	AffinityRow,
	
	/**
		\brief type of cell
	*/
	((Text)(type))
	
	/**
		\brief ID of cell
	*/
	((Integer)(cell_id))
	
	/**
		\brief ID of infection
	*/
	((Integer)(infection_id))
	
	/**
		\brief row in epitope matrix
	*/
	((Integer)(epitope_row))
	
	/**
		\brief column in epitope matrix
	*/
	((Integer)(epitope_column))
	
	/**
		\brief binding energy
	*/
	((Real)(energy))
	
	/**
		\brief log affinity
	*/
	((Real)(log_affinity))
)

/**
	\var affinityStats
	\brief Mean/standard deviation of memory and plasma cell log-affinity over time.
*/
ZPPDB_DEFINE_ROW_TYPE(
	AffinityStatsRow,
	
	/**
		\brief type of cell (memory or plasma)
	*/
	((Text)(type))
	
	/**
		\brief ID of infection
	*/
	((Integer)(infection_id))
	
	/**
		\brief GC round
	*/
	((Integer)(gc_round))
	
	/**
		\brief number of cells in calculation; plasma cells are only included for the epitome they bind best to
	*/
	((Integer)(n_cells))
	
	/**
		\brief row in epitope matrix
	*/
	((Integer)(epitope_row))
	
	/**
		\brief column in epitope matrix
	*/
	((Integer)(epitope_column))
	
	/**
		\brief mean log affinity
	*/
	((Real)(log_affinity_mean))
	
	/**
		\brief standard deviation of log affinity
	*/
	((Real)(log_affinity_sd))
	
	/**
		\brief number of cells that have maximum affinity to this epitope
	*/
	((Integer)(max_n_cells))
	
	/**
		\brief mean log affinity among cells with maximum affinity to this epitope
	*/
	((Real)(max_log_affinity_mean))
	
	/**
		\brief sd log affinity among cells with maximum affinity to this epitope
	*/
	((Real)(max_log_affinity_sd))
)

/**
	\brief Row for effectiveAffinity table (effective affinity of cells within GCs).
*/
ZPPDB_DEFINE_ROW_TYPE(
	EffectiveAffinityRow,
	
	/**
		\brief ID of cell
	*/
	((Integer)(cell_id))
	
	/**
		\brief ID of infection
	*/
	((Integer)(infection_id))
	
	/**
		\brief GC round
	*/
	((Integer)(gc_round))
	
	/**
		\brief row of epitope in epitope matrix
	*/
	((Integer)(epitope_row))
	
	/**
		\brief column of epitope in epitope matrix
	*/
	((Integer)(epitope_column))
	
	/**
		\brief effective affinity to this epitope
	*/
	((Real)(log_affinity))
)

/**
	\brief Row for plasmaConcentrations table (concentration of plasma cells over time).
*/
ZPPDB_DEFINE_ROW_TYPE(
	PlasmaConcRow,
	
	/**
		\brief ID of cell
	*/
	((Integer)(cell_id))
	
	/**
		\brief ID of infection
	*/
	((Integer)(infection_id))
	
	/**
		\brief GC round
	*/
	((Integer)(gc_round))
	
	/**
		\brief short-lived plasma cell concentration
	*/
	((Real)(log_conc_short))
	
	/**
		\brief long-lived plasma cell concentration
	*/
	((Real)(log_conc_long))
)

/**
	\brief Row for antigenNeighbor table: the neighbors of each locus in the energy
	landscape for each epitope in each infection.
*/
ZPPDB_DEFINE_ROW_TYPE(
	AntigenNeighborRow,
	
	/**
		\brief ID of infection
	*/
	((Integer)(infection_id))
	
	/**
		\brief row in epitope matrix
	*/
	((Integer)(epitope_row))
	
	/**
		\brief column in epitope matrix
	*/
	((Integer)(epitope_column))
	
	/**
		\brief the locus
	*/
	((Integer)(locus))
	
	/**
		\brief the order of the neighbor
	*/
	((Integer)(neighbor_index))
	
	/**
		\brief the neighboring locus
	*/
	((Integer)(neighbor_locus))
)

/**
	\brief Row for antigenEnergy table (binding energies for neighbor sequences).
	
	Binding energies, for each epitope in each infection, of all possible
	"neighbor sequences" consisting of the value of a locus and the values
	of its neighbors.
*/
ZPPDB_DEFINE_ROW_TYPE(
	AntigenEnergyRow,
	
	/**
		\brief ID of infection
	*/
	((Integer)(infection_id))
	
	/**
		\brief row in epitope matrix
	*/
	((Integer)(epitope_row))
	
	/**
		\brief column in epitope matrix
	*/
	((Integer)(epitope_column))
	
	/**
		\brief the locus
	*/
	((Integer)(locus))
	
	/**
		\brief the sequence consisting of locus + neighbors
	*/
	((Text)(neighbor_sequence))
	
	/**
		\brief the binding energy
	*/
	((Real)(energy))
)

#define DB_TABLES \
((MetaRow)(meta)) \
((TestStatsRow)(testStats)) \
((GCStatsRow)(gcStats)) \
((NaiveAlleleRow)(naiveAlleles)) \
((BCellRow)(bcells)) \
((SequenceRow)(sequences)) \
((GCFDCellRow)(gcfdCells)) \
((GCSeedCandidateRow)(gcSeedCandidates)) \
((GCSeedRow)(gcSeeds)) \
((AffinityRow)(affinity)) \
((AffinityStatsRow)(affinityStats)) \
((EffectiveAffinityRow)(effectiveAffinity)) \
((PlasmaConcRow)(plasmaConcentrations)) \
((AntigenNeighborRow)(antigenNeighbors)) \
((AntigenEnergyRow)(antigenEnergies))

#define DECLARE_TABLE_ITERATION(r, data, tableSpec) \
	DECLARE_TABLE( \
		BOOST_PP_SEQ_ELEM(0, tableSpec), \
		BOOST_PP_SEQ_ELEM(1, tableSpec) \
	)
#define DECLARE_TABLE(RowType, tableName) Table<RowType> BOOST_PP_CAT(tableName, Table);

#define TABLE(tableName) tableName ## Table

#endif // #ifndef DOXYGEN_SHOULD_SKIP_THIS

class DatabaseManager
{
public:
	BOOST_PP_SEQ_FOR_EACH(DECLARE_TABLE_ITERATION, data, DB_TABLES)
		
	DatabaseManager(zppdb::Database & db, Map<Bool> & dbTablesEnabledParams);
	void beginTransaction();
	void commit();
	
	template<typename Row>
	bool tableExists(Table<Row> & table)
	{
		return dbPtr->tableExists(table);
	}
	
	template<typename Row>
	void insert(Table<Row> & table, Row & row)
	{
		if(dbPtr->tableExists(table)) {
			dbPtr->insert(table, row);
		}
	}
private:
	std::unordered_map<std::string, bool> dbTablesEnabled;
	zppdb::Database * dbPtr;
};

#endif /* defined(__bcellmodel__DatabaseManager__) */
