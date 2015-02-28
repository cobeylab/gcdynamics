#include <iostream>
#include <unordered_map>
#include <map>
#include "DatabaseManager.h"


#define INIT_TABLE_ITERATION(r, data, tableSpec) \
	INIT_TABLE(BOOST_PP_SEQ_ELEM(1, tableSpec))
#define INIT_TABLE(tableName) BOOST_PP_CAT(tableName, Table) ( BOOST_PP_STRINGIZE(tableName) ),

#define CREATE_TABLE_ITERATION(r, data, tableSpec) \
	CREATE_TABLE(BOOST_PP_SEQ_ELEM(1, tableSpec))
#define CREATE_TABLE(tableName) if(dbTablesEnabled[BOOST_PP_STRINGIZE(tableName)]) { \
	dbPtr->createTable(BOOST_PP_CAT(tableName, Table)); \
}

DatabaseManager::DatabaseManager(zppdb::Database & db, Map<Bool> & dbTablesEnabledParams) :
	BOOST_PP_SEQ_FOR_EACH(INIT_TABLE_ITERATION, data, DB_TABLES)
	dbPtr(&db)
{
	for(std::string key : dbTablesEnabledParams.keys()) {
		dbTablesEnabled[key] = dbTablesEnabledParams[key];
	}
	
	BOOST_PP_SEQ_FOR_EACH(CREATE_TABLE_ITERATION, data, DB_TABLES)
}
	
void DatabaseManager::beginTransaction()
{
	dbPtr->beginTransaction();
}

void DatabaseManager::commit()
{
	dbPtr->commitWithRetry(0.1, 100, std::cerr);
}
