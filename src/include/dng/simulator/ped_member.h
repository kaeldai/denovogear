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

#ifndef DNG_SIM_PED_MEMBER_H_
#define DNG_SIM_PED_MEMBER_H_

#include <unordered_map>
#include <string>
#include <vector>
#include <boost/format.hpp>
#include <unordered_map>
#include "dng/io/ped.h"

#define ID_UNKNOWN -1
#define IS_UNKNOWN(mem) ((mem)->mid == ID_UNKNOWN)
#define IS_EMPTY(mem) (mem == nullptr)

#define PED_NOPARENT "0"

namespace dng {
namespace sim {



typedef size_t member_id;
typedef size_t family_id;
typedef dng::io::Pedigree::Gender Gender;

enum Base : uint8_t { A = 0, C, G, T, N, REF, UNASSIGNED };

static char nt2char[] = {'A', 'C', 'G', 'T', 'N', 'R', 'U'};
static Base index2Base[] = {A, C, G, T, N};

typedef std::unordered_map<size_t, Base> reference_map;



/**
 * Stores information about a sample/library. Including
 */
struct Sample {
  //std::string id;
  std::string lb;
  std::string sm;
  
  //std::string name; // library name
	size_t depth; // number of reeds per site.
	reference_map dna[2]; // list of somatic mutation diffs between the reference and sample

  
	Base get_nt(int chrom, size_t site) {
		reference_map &rmap = dna[chrom];
		reference_map::const_iterator nt = rmap.find(site);
		if(nt ==  rmap.end()) {
			return REF;
		}
		else {
			return nt->second;
		}
	}

  std::string id_vcf() {
    if(lb == sm)
      return lb;
    else
      return sm + ":" + lb;
  }

  std::string id_sam() {
    if(lb == sm)
      return lb;
    else 
      return sm + "-" + lb;
  }
};

struct Member {
	member_id mid;
	std::string name; // name that shows up in pedigree and data file
	family_id fid;
	Gender sex;
	Member *mom;
	Member *dad;
	bool hasDNA; // Used when generating the pedigree contigs
  
  //std::shared_ptr<Member> mom1;
  
	// Used to keep track of differences between member's DNA and the reference. We need one for each chromatid and somatic vs gametic
	reference_map gametes[2];
	std::vector<Sample> libraries;

	Member(member_id mid_, family_id fid_, Gender sex_) : mid(mid_), fid(fid_), sex(sex_) {
		name = (boost::format("NA%04d") % mid).str();
		mom = nullptr;
		dad = nullptr;
		hasDNA = false;
	}

	std::string get_family_name() {
		return (boost::format("F%03d") % fid).str();
	}

	std::string get_dad_id() {
		if(dad == nullptr)
			return 0;
		else
			return dad->name;
	}

	std::string get_mom_id() {
		if(mom == nullptr)
			return std::string("0");
		else
			return mom->name;
	}

	/**
	* Change the nucleotide at 'site' for chromosome pair number 'chrom_num' to 'base'.
	*/
	void update_gamete_dna(int chrom_num, size_t site, Base base) {
		if(chrom_num >= 2)
			throw std::runtime_error("Attempting to set base from more than diploidy DNA.");

		gametes[chrom_num].insert({site, base});
	}

	/**
	* Makes copies of mom and dad's dna and saves into member. Selection of chromosome pair is random
	*/
	void inherit_dna(Member *mom, Member *dad) {
		// randomly select one of each pararents contigs and copy to child
		int mom_chrom_num = rand() % 2;
		int dad_chrom_num = rand() % 2;
		int placement = rand() % 2;

		gametes[placement] = reference_map(mom->gametes[mom_chrom_num]);
		gametes[(placement+1)%2] = reference_map(dad->gametes[dad_chrom_num]);
	}

		/**
		* Get the nucleotide base at specified site, from either chromatid 0 or 1.
		*/
	Base get_gamete_nt(size_t chromatid, size_t site) {
		if(chromatid >= 2)
			throw std::runtime_error("Attempting to get nuclitide from chromatid " + std::to_string(chromatid) + ". Members only defined diploidy DNA.");

		reference_map &chrom = gametes[chromatid];
		reference_map::const_iterator nt = chrom.find(site);
		if(nt == chrom.end()) {
			return REF;
		}
		else {
			return nt->second;
		}
	}

	void add_library(const std::string &lb, size_t depth) {
		// TODO: Check that library with same name doesn't already exists

		Sample lib;
		/*
		if(lb != name) {
		  lib.id = name + "-" + lb;
		}
		else {
		  lib.id = name; 
		}
		*/
		lib.lb = lb;
		lib.sm = name;
		//lib.name = std::string(name);
		lib.depth = depth;
		lib.dna[0] = gametes[0];
		lib.dna[1] = gametes[1];
		libraries.push_back(lib);
	}

	Base get_lib_dna(size_t l_indx, size_t chrom, size_t site) {
		reference_map &rmap = libraries[l_indx].dna[chrom];
		reference_map::const_iterator nt = rmap.find(site);
		if(nt ==  rmap.end()) {
			return REF;
		}
		else {
			return nt->second;
		}
	}

	void update_lib_dna(size_t l_indx, size_t chrom, size_t site, Base base) {
		libraries[l_indx].dna[chrom].insert({site, base});
	}
};



} // namespace sim
} // namespace dng



#endif /* PED_MEMBER_H_ */
