//
//  Sequence.cpp
//  bcellmodel
//
//  Created by Ed Baskerville on 2/3/14.
//  Copyright (c) 2014 Cobey Lab. All rights reserved.
//

#include "Sequence.h"

using namespace std;

//unordered_map<vector<bool>, shared_ptr<vector<bool>>> Sequence::bitStore = unordered_map<vector<bool>, shared_ptr<vector<bool>>>();

bool operator==(Sequence const & s1, Sequence const & s2)
{
	return s1.equal_to(s2);
}
