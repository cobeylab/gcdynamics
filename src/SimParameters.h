/**
	Model parameters. All values are loaded from external JSON file.
*/

#ifndef __bcellmodel__SimParameters__
#define __bcellmodel__SimParameters__

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#include <map>
#include "zppjson.hpp"

using namespace zppjson;

/**
	\brief JSON enum type defining sequence segment type.
*/
ZPPJSON_DEFINE_ENUM(
	SegmentType,
	
	/**
		\brief Allelic segment type.
	*/
	(ALLELE)
	
	/**
		\brief Random insert segment type.
	*/
	(RANDOM_INSERT)
)

/**
	\brief Parameters controlling a single segment of a sequence.
	
	Multiple segments are combined to form chains, which are combined to produce
	a complete sequence for a B cell.
*/
ZPPJSON_DEFINE_TYPE(
	SegmentParams,
	
	/**
		\brief Either ALLELE or RANDOM_INSERT.
	*/
	((SegmentType)(type))
	
	/**
		\brief The minimum length of the sequence.
		
		The minimum length of the sequence is used to determine the meaning of
		`lengthDistribution`, so that the weight in entry 0 of
		`lengthDistribution` corresponds to `minLength`; entry `i` corresponds
		to `minLength + i`.
	*/
	((Int64)(minLength))
	
	/**
		\brief Relative probabilities of different sequence lengths.
		
		These relative probabilites need not add up to 1; they are normalized
		when used. Entry `i` corresponds to length `minLength + i`.
	*/
	((Array<Double>)(lengthDistribution))
	
	/**
		\brief Relative probabilities of different alleles.
	*/
	((Array<Double>)(alleleDistribution))
)

/**
	\brief Parameters controlling a single chain of a sequence.
	
	Chains are made of segment; multiple chains are combined to form a sequence.
*/
ZPPJSON_DEFINE_TYPE(
	ChainParams,
	/**
		\brief The truncated total length of a chain in the sequence.
		
		During sequence generation, segments are attached end-to-end and then
		truncated at constant length `truncatedLength`.
	*/
	((Int64)(truncatedLength))
	
	/**
		\brief Parameters for segments contained in this chain.
	*/
	((Array<SegmentParams>)(segments))
)

/**
	\brief Parameters controlling a single infection.
*/
ZPPJSON_DEFINE_TYPE(
	InfectionParams,
	
	/**
		\brief Delay between the start of this infection and the start of the next.
		
		The delay must exceed `gcs.nRounds * gcs.roundDuration`, the
		time it takes for the completion of all GCs following an infection.
	*/
	((Double)(delay))
	
	/**
		\brief Quantity of antigen for this infection.
	*/
	((Double)(antigenQuantity))
)

/**
	\brief Parameters controlling all epitopes.
	
	Rather than having an object for each epitope, these parameters are
	implemented as parallel matrices that can contain either one entry,
	which will be used for all epitopes, or multiple entries, one for each
	epitope.
*/
ZPPJSON_DEFINE_TYPE(
	EpitopeParams,
	
	/**
		\brief The number of rows in the epitope lattice.
	*/
	((Int64)(rows))
	
	/**
		\brief The number of columns in the epitope lattice.
	*/
	((Int64)(cols))
	
	/**
		\brief The number of interacting loci, per locus.
		
		This number affects the interactions present in the energy landscape
		generated for each epitope.

	*/
	((Int64)(nInteractions))
	
	/**
		\brief The seed value used to uniquely determine the neighbors in the
		energy landscape.
		
		If not present or zero, a random value will be generated.
	*/
	((Array<Array<Int64>>)(neighborSeed))
	
	/**
		\brief The seed value used to uniquely determine the energies in the
		energy landscape.
		
		If not present or zero, a random value will be generated.
	*/
	((Array<Array<Int64>>)(energySeed))
	
	/**
		\brief The probability that an energy changes from the last antigen.
		
		This number affects the creation of epitopes during secondary
		infections, which have their energy landscape created as mutations of
		previous energy landscapes.
	*/
	((Array<Array<Double>>)(pEnergyMutation))
	
	/**
		\brief The probability that a neighbor changes from the last antigen.
		
		See `pEnergyMutation`.
	*/
	((Array<Array<Double>>)(pNeighborMutation))
	
	/**
		\brief Mean of energy distribution.
		
		The energy distribution is used to draw values for different
		locus+interacting-neighbor sequences. The standard deviation is 1.0.
	*/
	((Array<Array<Double>>)(energyMean))
	
	/**
		\brief Parameter `A` in energy-to-affinity mapping.
		
		Log-affinity is calculated as `affinityA - affinityB * energy`.
	*/
	((Double)(affinityA))
	
	/**
		\var affinityb
		\brief Parameter `B` in energy-to-affinity mapping.
		
		Log-affinity is calculated as `affinityA - affinityB * energy`.
	*/
	((Double)(affinityB))
)

/**
	\brief Parameters controlling B-cell sequences.
*/
ZPPJSON_DEFINE_TYPE(
	SequenceParams,
	
	/**
		\brief Number of amino acids in alphabet.
	*/
	((Int64)(alphabetSize))
	
	/**
		\brief Characters used to represent amino acid.
	*/
	((String)(alphabet))
	
	/**
		\brief Vector parameters controlling chains that make up sequence.
	*/
	((Array<ChainParams>)(chains))
)

/**
	\brief Parameters controlling germinal centers.
*/
ZPPJSON_DEFINE_TYPE(
	GCParams,
	
	/**
		\brief Number of germinal centers created in response to infections.
	*/
	((Int64)(nGCs))
	
	/**
		\brief Number of rounds of B-cell affinity maturation for GCs.
	*/
	((Int64)(nRounds))
	
	/**
		\brief Simulation time for a single round of affinity maturation.
	*/
	((Double)(roundDuration))
	
	/**
		\brief Whether or not to use a null model where affinity does not affect replication
	*/
	((Bool)(useNullReplicationModel))
	
	/**
		\brief Whether or not to use old GC seed method (generate lots of cells)
	*/
	((Bool)(useOldSeedMethod))
	
	/**
		\brief Fraction of naive cell/affinity pairs to seed GCs with.
		
		This parameter applies only to the first infection course, before a pool
		of naive + memory cells has been created. Each cell/affinity pair is
		calculated, and the top `fNaiveSeed` fraction of them, rounded up to the
		nearest integer, is used to construct a list of GC seed candidates.
	*/
	((Double)(fNaiveSeed))
	
	/**
		\brief Whether or not to exclude memory cells from seeding GCs.
		
		This parameter applies only to secondary infections. If `true`, memory
		cells are excluded from the choice of GC seeds, but secondary infections
		will be chosen from `[seedQuantileLower, seedQuantileUpper]` naive
		cells.
	*/
	((Bool)(seedExcludeMemory))
	
	/**
		\brief The lower-bound quantile for seeding GCs.
		
		This parameter applies only to secondary infectino courses. Each
		cell/affinity pair is calculated, and entries in
		`[seedQuantileLower, seedQuantileUpper]` are eligible to be GC seeds.
		At least one cell is always deemed to be eligible.
	*/
	((Double)(seedQuantileLower))
	
	/**
		/brief The upper-bound quantile for seeding GCs.
		
		See `seedQuantileLower`.
	*/
	((Double)(seedQuantileUpper))
	
	/**
		\brief Probability of somatic hypermutation, per B-cell locus.
	*/
	((Double)(pMutation))
	
	/**
		\brief Number of cells used to seed each GC.
	*/
	((Int64)(nSeedCells))
	
	/**
		\brief Seed cells for each GC: if present, used instead of randomized
		seed generation procedure.
	*/
	((Array<Array<String>>)(seedCells))
	
	/**
		\brief Total number of cells in each GC after replication.
	*/
	((Int64)(nTotalCells))
	
	/**
		\brief Fraction of cells chosen as winners after replication.
		
		Number of cells is rounded up to the nearest integer.
	*/
	((Double)(fWinners))
	
	/**
		\brief Whether or not cell populations should grow exponentially in GC.
		
		If true, during each GC round the number of cells in the GC after
		replication will be `nCells = min(nTotalCells, 2^nDivisionsPerRound)`.
		The number of cells will then be reduced to `ceil(fWinners * nCells)`
		via selection.
	*/
	((Bool)(exponentialGrowthActive))
	
	/**
		\brief Number of cell divisions per round.
		
		See exponentialGrowthActive for details.
	*/
	((Int64)(nDivisionsPerRound))
	
	
	/**
		\brief Whether or not GCs should be terminated when the highest
		log-affinity drops below a threshold.
	*/
	((Bool)(shouldTerminateBelowThreshold))
	
	/**
		\brief Threshold below which GCs should be terminated.
		
		This parameter only applies if `shouldTerminateBelowThreshold == true`.
	*/
	((Double)(logAffinityTerminationThreshold))
	
	/**
		\brief Log-affinity threshold above which cells are eligible for export
		as plasma cells.
	*/
	((Double)(plasmaExportLogAffinityThreshold))
	
	/**
		\brief Probability that an eligible cell will be exported as a plasma
		cell, per GC round.
	*/
	((Double)(pPlasmaExport))
	
	/**
		\brief The probability that a cell will be exported as a memory cell
		in the first GC round. Intermediate rounds will have probabilities
		linearly inteprolated between start and end rounds.
	*/
	((Double)(pMemoryExportFirstRound))
	
	/**
		\brief The probability that a cell will be exported as a memory cell
		in the last GC round. Intermediate rounds will have probabilities
		linearly inteprolated between start and end rounds.
	*/
	((Double)(pMemoryExportLastRound))
	
	/**
		\brief Coefficient for antibody masking.
		
		E_j = F_j exp(-alpha * FAb_j / F_j)
		
		where E_j is the effective affinity of the cell to epitope j;
		F_j is the intrinsic affinity of the cell to epitope j;
		and FAb_j is the intrinsic affinity of the masking antibody to epitope j.
	*/
	((Double)(alphaMasking))
)

/**
	\brief Parameters controlling plasma cells.
	
	Two instances of this class exist in SimParameters, one for short-lived
	and one for long-lived cells.
*/
ZPPJSON_DEFINE_TYPE(
	PlasmaParams,
	
	/**
		\brief Initial log-concentration of this type of plasma cell.
	*/
	((Double)(initialLogC))
	
	/**
		\brief The log-concentration at which this type of plasma cell
		is discarded from the simulation.
	*/
	((Double)(discardLogC))
	
	/**
		\brief The rate at which concentration decays.
		
		The concentration is multiplied by a factor of `decayRate` per unit
		time, so that the `C(1) == decayRate * C(0)`. That is, the log-
		concentration decreases linearly by `log(decayRate)` each unit of time.
	*/
	((Double)(decayRate))
)

/**
	\brief Top-level parameters class.
*/
ZPPJSON_DEFINE_TYPE(
	SimParameters,
	
	/**
		\brief Seed for random number generator.
		
		If absent or set to 0, a seed will be generated using `std::random_device`.
	*/
	((Int64)(randomSeed))
	
	/**
		\brief Filename for the database.
		
		If set to the empty string, no database will be created.
	*/
	((String)(dbFilename))
	
	/**
		\brief Dictionary of which database tables are enabled.
	*/
	((Map<Bool>)(dbTablesEnabled))
	
	/**
		\brief Dictionary of which cell types have affinity output enabled.
	*/
	((Map<Bool>)(affinityOutputTypesEnabled))
	
	/**
		\brief Number of naive cells to create.
	*/
	((Int64)(nBCellsNaive))
	
	/**
		\brief Parameters controlling sequences.
	*/
	((SequenceParams)(sequences))
	
	/**
		\brief Parameters controlling epitopes
	*/
	((EpitopeParams)(epitopes))
	
	/**
		\brief Parameters controlling germinal centers.
	*/
	((GCParams)(gcs))
	
	/**
		\brief Vector of parameters for infections, one per infection.
	*/
	((Array<InfectionParams>)(infections))
	
	/**
		\brief Parameters controlling short-lived plasma cells.
	*/
	((PlasmaParams)(shortLivedPlasma))
	
	/**
		\brief Parameters controlling long-lived plasma cells.
	*/
	((PlasmaParams)(longLivedPlasma))
)

#endif // #ifndef DOXYGEN_SHOULD_SKIP_THIS

/*** VERIFICATION METHODS ***/

void verify(SimParameters & p);
void verify(InfectionParams & p, double minDelay);
void verify(EpitopeParams & p);
void verify(PlasmaParams & p);
void verify(GCParams & p);
void verify(SequenceParams & p);
void verify(ChainParams & p);
void verify(SegmentParams & p);

#endif // #ifndef __bcellmodel__SimParameters__
