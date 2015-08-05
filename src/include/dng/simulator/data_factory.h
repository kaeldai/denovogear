/*
 * Copyright (c) 2014 Reed A. Cartwright
 * Authors:  Reed A. Cartwright <reed@cartwrig.ht>, Kael Dai <kdai1@asu.edu>
 *
 * This file is part of DeNovoGear.
 *
 * DeNovoGear is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef DNG_SIM_DATA_FACTORY_H
#define DNG_SIM_DATA_FACTORY_H


#include "dng/hts/bcf.h"
#include "dng/hts/bam.h"
#include <vector>

#include <set>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <dng/io/ped.h>
#include <boost/format.hpp>
#include <random>
#include <math.h>
#include <dng/simulator/bam.h>
#include <dng/simulator/ped_member.h>
#include <dng/simulator/models.h>
#include <array>
#include <boost/algorithm/string.hpp>
/*
#define ID_UNKNOWN -1
#define IS_UNKNOWN(mem) ((mem)->mid == ID_UNKNOWN)
#define IS_EMPTY(mem) (mem == nullptr)
*/

//#define PED_NOPARENT "0"




namespace dng {
namespace sim {
    

typedef std::array<double, 16> genotype_dist;

// Format of output data file. 'None' indicates only PED file is generated.
//enum DataFormat { VCF, BCF, SAM, BAM, CRAM, None };

//static char nt2char[] = {'A', 'C', 'G', 'T', 'N', 'R', 'U'};

typedef std::unique_ptr<SimBuilder> sim_ptr;

#include <memory>

static std::string model_dng = "DNG";
static std::string model_test1 = "Test1";

class Factory {
public:


	static sim_ptr newInstance(std::string &model) {
		if(model == model_dng) {
			return sim_ptr(new DNGModel());
		}
		else if(model == model_test1) {
			return sim_ptr(new Test1());
		}
		else{
			throw std::runtime_error("Unknown model type '" + model + "'. ");
		}

	}

private:
	Factory() {

	}

};


/*
class Factory1 {
public:
private:

};
*/
  
} // namespace sim
} // namespace dng

#endif
