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
#define MEMBER_NA ""


#define PED_NOPARENT "0"

namespace dng {
namespace sim {



  //typedef size_t member_id;
//typedef size_t family_id;
typedef dng::io::Pedigree::Gender Gender;

enum Base : uint8_t { A = 0, C, G, T, N, REF, UNASSIGNED };

static char nt2char[] = {'A', 'C', 'G', 'T', 'N', 'R', 'U'};
static Base index2Base[] = {A, C, G, T, N};

typedef std::unordered_map<size_t, Base> reference_map;

struct Member;
typedef std::shared_ptr<Member> member_ptr;
 

/**
 * Stores information about a sample library. Including VCF/SAM header tags, allele depth size, and list of variant locations.
 */
struct Library {
  std::string lb; // SAM LB tag 
  std::string sm; // SAM SM tag
  std::string id; // SAM ID tag
  std::string sample_name; // VCF genotype header
  
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
};

struct Member {
  std::string id; // identifier for member of pedigree
  std::string family_id; // family unit member belongs to (first col in PED file)
  Gender sex;
  member_ptr mom_ptr; // reference to member's mom
  member_ptr dad_ptr; // reference to member's dad

  bool hasDNA; // Used when generating the pedigree contigs
  
  // Used to keep track of differences between member's DNA and the reference. One for each chromatid and somatic vs gametic
  reference_map gametes[2];
  std::vector<Library> libraries;

  Member(const std::string &id_, const std::string &fam_id_, Gender sex_) : id(id_), family_id(fam_id_), sex(sex_) {
    mom_ptr = nullptr;
    dad_ptr = nullptr;
    hasDNA = false;
  }

  
  // Change the nucleotide at 'site' for chromosome pair number 'chrom_num' to 'base'.
  void update_gamete_dna(int chrom_num, size_t site, Base base) {
    if(chrom_num >= 2)
      throw std::runtime_error("Attempting to set base from more than diploidy DNA.");
    
    gametes[chrom_num].insert({site, base});
  }

  
  // Makes copies of mom and dad's dna and saves into member. Selection of chromosome pair is random
  void inherit_dna(member_ptr mom, member_ptr dad) {
    // randomly select one of each pararents contigs and copy to child
    int mom_chrom_num = rand() % 2;
    int dad_chrom_num = rand() % 2;
    int placement = rand() % 2;
    
    gametes[placement] = reference_map(mom->gametes[mom_chrom_num]);
    gametes[(placement+1)%2] = reference_map(dad->gametes[dad_chrom_num]);
  }

  // Get the nucleotide base at specified site, from either chromatid 0 or 1.
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
  
  void add_library(const std::string &lib_id, size_t depth) {
    // Check that library with same name doesn't already exists
    for(Library l : libraries) {
      if(l.lb == lib_id) {
	std::cerr << "WARNING. Member " + id + " already has a library with lb tag " + lib_id + ". Skipping!" << std::endl;
	return;
      }
    }
    
    Library lib;
    lib.lb = lib_id;
    lib.sm = id;
    if(lib_id == id) {
      lib.id = id;
      lib.sample_name = id;
    }
    else {
      lib.id = lib.sm + "-" + lib.lb;
      lib.sample_name = lib.sm + ":" + lib.lb;
    }
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
