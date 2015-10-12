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

#ifndef DNG_SIM_SIMULATOR_H_
#define DNG_SIM_SIMULATOR_H_

#include <dng/hts/bcf.h>
#include <dng/simulator/ped_member.h>
#include <dng/simulator/bam.h>
#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include "htslib/faidx.h"

namespace dng {
namespace sim {

// Output format type
enum SeqFormat { SAM, BAM, VCF, BCF };

// How to divided the output data files
enum DataScheme {
	SINGLE_FILE, // All data gets put into a single file
	PER_MEMBER, // One data file per member - with potentially multiple libraries per file
	PER_LIBRARY // each library sample is put into a separate file
};

std::array<std::pair<Base, Base>, 16> genotypes = {{{A, A}, {A, C}, {A, G}, {A, T},
						    {C, A}, {C, C}, {C, G}, {C, T},
						    {G, A}, {G, C}, {G, G}, {G, T},
						    {T, A}, {T, C}, {T, G}, {T, T}}};
 
typedef std::map<std::string, member_ptr> pedigree_map;
 

class SimBuilder {

public:

  
  void AddTrio(const std::string &child, const std::string &mom, const std::string &dad, const std::string &family = "F1") { 
    Gender gen = Gender::Unknown;
    pedigree_map::iterator c = pedigree.find(child);
    if(c != pedigree.end())
      gen = (c->second)->sex;
    
    AddTrio(child, gen, mom, dad, family);
  }


  
	void AddTrio(const std::string &child, Gender sex, const std::string &mom, const std::string &dad, const std::string &family = "F1") {
		member_ptr c = buildMember(child, sex, family);
		member_ptr m = buildMember(mom, Gender::Female, family);
		member_ptr d = buildMember(dad, Gender::Male, family);

		if(c == nullptr) {
			// Throw error, Can't have a trio without a child
			throw std::runtime_error("Attempting to create a trio without a child.");
		}

		if(m != nullptr) {
			if(c->mom_ptr != nullptr && c->mom_ptr->id != m->id) {
				throw std::runtime_error("Failed to create child-mother relationship between " + c->id + " and " + m->id +
										 ". The child already has mother " + c->mom_ptr->id);
			} else {
				c->mom_ptr = m;
				gametic_nodes_++;
			}
		}

		if(d != nullptr) {
			if(c->dad_ptr != nullptr && c->dad_ptr->id != d->id) {
				throw std::runtime_error("Failed to create child-father relationship between " + c->id + " and " + d->id +
										 ". The child already has mother " + c->dad_ptr->id);
			}
			else {
				c->dad_ptr = d;
				gametic_nodes_++;
			}
		}
	}

  /**
   *  Adds library to all members of the pedigree. LB and ID tags are the same as the sample 
   *   depth: Read depth
   */ 
  void SetDefaultLibraries(size_t depth = 10) {
    for(pedigree_map::const_iterator itr = pedigree.begin(); itr != pedigree.end(); ++itr) {
      AddLibrary(itr->first, itr->first, depth);
    }
  }

  /**
   * Add library of LB=lib_id to all members of the pedigree.
   *  lib_id: library name, unique for each member
   *  depth: Read depth
   */
  void SetDefaultLibraries(const std::string &lib_id, size_t depth = 10) {
    for(pedigree_map::const_iterator itr = pedigree.begin(); itr != pedigree.end(); ++itr) {
      AddLibrary(itr->first, lib_id, depth);
    }
  }

  /**
   * Add library, LB = lib_id, to member of the pedigree.
   *  member_id: id of pedigree member
   *  lib_id: name of library, must be unique for each member
   *  detph: read depth
   */
  void AddLibrary(const std::string &member_id, const std::string &lib_id, size_t depth = 10) {
      
    pedigree_map::const_iterator itr = pedigree.find(member_id);
    if(itr == pedigree.end()) {
      std::runtime_error("Member " + member_id + " does not exist in the pedigree, cannot add library.");
    }
    else {
      itr->second->add_library(lib_id, depth);
    } 
           
  }
	
  virtual void setParameter(std::string &name, int val) {
    if(name == "sam_seq_len") {
      sam_seq_len_ = val;
    }
    else {
      std::cerr << "Unrecognized parameter '" + name + "' (integer). Skipping!" << std::endl;
    }
  }


  virtual void setParameter(std::string &name, double val) {
    std::cerr << "Unrecognized parameter '" + name + "' (double). Skipping!" << std::endl;
  }

  virtual void setParameter(std::string &name, std::vector<double> val) {
    std::cerr << "Unrecognized parameter '" + name + "' (double[]). Skipping!" << std::endl;
  }


  // Publish a PED file containing all the pedigree information.
  template<typename Stream>
  void publishPed(Stream &output, bool newick_formatted = false) {

    for(pedigree_map::const_iterator itr = pedigree.begin(); itr != pedigree.end(); ++itr) {
      Member *m = itr->second.get();
      std::string mom = (m->mom_ptr == nullptr ? PED_NOPARENT : m->mom_ptr->id);
      std::string dad = (m->dad_ptr == nullptr ? PED_NOPARENT : m->dad_ptr->id);
      
      output << m->family_id << "\t"
	     << m->id << "\t"
	     << dad << "\t"
	     << mom << "\t"
	     << getGender(m->sex) << "\t"
	     // << m->family_id + m->id
	     << std::endl;
    }
  }

  // Publish sequence or variant pedigree data to a file
  //  fname = output file name
  //  format = currently supports BAM, SAM, BCF, VCF
  //  scheme = number of files to put write to: SINGLE_FILE or PER_MEMBER
  //
  void publishData(const std::string &fname, SeqFormat format, DataScheme scheme = SINGLE_FILE) {
    createSeqData();
    
    if(scheme == SINGLE_FILE) {
      // Create a single file
      switch(format) {
      case BAM: publishDataSAM(fname.c_str(), "wb", GetPedigree()); break;
      case SAM: publishDataSAM(fname.c_str(), "w", GetPedigree()); break;
      case BCF: publishDataVCF(fname.c_str(), "wb", GetPedigree()); break;
      case VCF: publishDataVCF(fname.c_str(), "w", GetPedigree()); break;
      }
    }
    else if(scheme == PER_MEMBER){
      // creating multiple files, turn fname.bam --> fname_{MID}.bam
      std::string base = fname;
      std::string postfix = "";
      size_t indx = fname.find_last_of(".");
      if(indx != std::string::npos) {
	base = fname.substr(0, indx);
	postfix = fname.substr(indx);
      }
      
      // create a file for each member in the pedigree
      for(member_ptr m : GetPedigree()) {
	std::string file = base + "_" +  m->id + postfix;
	std::vector<member_ptr> mems = {m};
	switch(format) {
	case BAM: publishDataSAM(fname.c_str(), "wb", mems); break;
	case SAM: publishDataSAM(fname.c_str(), "w", mems); break;
	case BCF: publishDataVCF(fname.c_str(), "wb", mems); break;
	case VCF: publishDataVCF(fname.c_str(), "w", mems); break;
	}
      }
    }
  }

  // Set the reference that will be used to build the pedigree gametic/somatic sequences.
  //  seq = A contiguous nucleotide sequence of A/C/T/G/N's
  //  chrom = name of contig/chromosome.
  //  start_pos = start location of sequence along contig
  virtual void setReference(std::string &seqence, std::string &chrom, size_t start_pos) {
    chrom_ = std::string(chrom);
    start_pos_ = start_pos;
    
    reference.clear(); 
    for(int nt : seqence) {
      reference.push_back(char2base(nt));
    }
    std::cout << reference.size() << std::endl;
  }


  // Set the reference that will be used to build the pedigree based on a subsequence (range)
  // of sequenced data taken from a fasta file.
  virtual void setReference(std::string &fasta, std::string &range) {
    fai_build(fasta.c_str());
    //std::string fai_file = fasta + ".fai";
    faidx_t *fai = fai_load(fasta.c_str());
    int len;
    char *char_ref = fai_fetch(fai, range.c_str(), &len);
    std::cout << len << std::endl;
    reference.clear();
    for(int a = 0; a < len; a++) {
      //std::cout << char_ref[a];
      reference.push_back(char2base(char_ref[a]));
    }
    delete char_ref;
    
    size_t chrom_indx = range.find(":");
    if(chrom_indx == std::string::npos) {
      chrom_ = std::string(range);
    }
    else {
      chrom_ = std::string(range.substr(0, chrom_indx));
    }
  }
  

  // Creates the library (i.e. somatic samples) sequences.
  virtual void createSeqData() = 0;

  
  // Mimics the reading of a base by a sequencer, returning a collection of definite nt's
  // depending on the read depth. Implementation will be dependent on the model.
  virtual std::vector<Base> baseCall(Library &lib, size_t site) = 0;
  

  // Returns reference to pedigree member
  member_ptr GetMember(const std::string &member_id) {
    pedigree_map::iterator itr = pedigree.find(member_id);
    if(itr == pedigree.end())
      return nullptr;
    else
      return itr->second;
  }

  // Returns a collection of all the members in the pedigree.
  // No specific order guarenteed.
  std::vector<member_ptr> GetPedigree() {
    std::vector<member_ptr> mems;
    mems.reserve(pedigree.size());
    for(pedigree_map::iterator itr = pedigree.begin(); itr != pedigree.end(); ++itr)
      mems.push_back(itr->second);
    
    return mems;
  }
  
  std::size_t gametic_nodes() {
	  return gametic_nodes_;
  }

    
  virtual ~SimBuilder() {  }
  
protected:

  // Returns a member if member_id already exists in the pedigree. If it doesn't create the new member given it's
  // id, sex, and family information.
  // If Member does alreay exists method will check that the sex and family_id are not different.
  member_ptr buildMember(const std::string &member_id, Gender sex, const std::string &family_id) {
    if(member_id == "") {
      return nullptr;
    }
    
    pedigree_map::iterator itr = pedigree.find(member_id);
    if(itr == pedigree.end()) {
      // Create a brand new member if no matching id is in the peidgree
      member_ptr m = std::make_shared<Member>(member_id, family_id, sex);
      pedigree.insert(std::pair<std::string, member_ptr>(member_id, m));
      return m;
    }
    else {
      // Member already exists, check for clashes.
      member_ptr m = itr->second;
      if(m->sex != sex) {
	throw std::runtime_error("Attempting to redifine the sex of pedigree member " + m->id + ". exiting!");
      }

      if(m->family_id != family_id) {
	std::cerr << "WARNING. Pedigree member " << m->id
		  << " has already been defined as a member of family " << m->family_id
		  << ". Skipping attempt to add to second family " << family_id << std::endl;
      }
      return m;
    }
  }

 
  // Convert Gender into a format used by .ped file
  int getGender(Gender sex) {
    int ret;
    switch(sex){
    case Gender::Male: ret = 1; break;
    case Gender::Female: ret = 2; break;
    default: ret = 0; break;
    }
    return ret;
  }

  virtual void publishDataSAM(const char *file, const char *mode, const std::vector<member_ptr> &mems) {

    // SAM Header
    std::stringstream hdr_txt;
    uint32_t chrom_len_ = static_cast<uint32_t>(2*start_pos_ + reference.size()); // take a guess of full contig length
    hdr_txt << "@HD\tVN:0.1\tSO:unknown\tGO:none" << std::endl;
    hdr_txt << "@SQ\tSN" << chrom_ << "\tLN:" << chrom_len_ << std::endl;
    for(member_ptr m : mems) {
      for(Library &lib : m->libraries) {
	hdr_txt << "@RG" << "\t"
		<< "ID:" << lib.id << "\t"
		<< "LB:" << lib.lb << "\t"
		<< "SM:" << lib.sm << std::endl;
      }
    }
    dng::sim::BAMFile out(file, mode);
    out.set_header(hdr_txt.str().c_str(), hdr_txt.str().size());      

    // Iterate through member of mems 
    size_t contig_len_ = 30;
    for(member_ptr member : mems) {
      for(int l = 0; l < member->libraries.size(); l++) {
	Library &library = member->libraries[l];

	for(int pos = 0; pos < reference.size(); pos += contig_len_) {
	  size_t l_ccontig = (pos + contig_len_ < reference.size() ? contig_len_ : reference.size() - pos);

	  // TODO: Figure out a better data structure for storing the results. Row base look-ups are not efficient.
	  std::vector<std::vector<Base>> reads;
	  for(int next = 0; next < l_ccontig; next++) {
	    reads.push_back(baseCall(library, next+pos));
	  }

	  for(int d = 0; d < library.depth; d++) {
	    dng::sim::BAMRec rec = out.init_rec();
	    rec.set_qname(chrom_);
	    rec.add_cigar(BAM_CMATCH, l_ccontig);

	    std::string seq;
	    for(int site = 0; site < l_ccontig; site++) {
	      seq += nt2char[reads[site][d]];
	    }

	    rec.set_seq(seq);

	    std::string qual = "*";
	    rec.set_qual(qual);

	    std::string aux = member->id + "-" + member->libraries[l].id;
	    char sm[2] = {'S', 'M'};
	    rec.add_aux(sm, 'Z', aux);

	    out.write_record(rec);
	  }
	}
      }
    }
    out.save();
  }

  virtual void publishDataVCF(const char *file, const char *mode, const std::vector<member_ptr> &mems) {

    // Create VCF header
    hts::bcf::File out(file, mode);
    out.AddHeaderMetadata("##fileformat-VCFv4.2");
    out.AddHeaderMetadata("##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic Depths\">");
    uint32_t chrom_len_ = static_cast<uint32_t>(2*start_pos_ + reference.size()); // take a guess of full contig length
      out.AddContig(chrom_.c_str(), chrom_len_);
    for(member_ptr m : mems) {
      for(Library &lib : m->libraries) {
	out.AddSample(lib.sample_name.c_str());
      }
    }

    out.WriteHeader();    
    
    // Go through each position in the reference and build a line in the VCF file
    size_t l_reference = reference.size();
    for(size_t pos  = 0; pos < l_reference; pos++) {
      Base ref = reference[pos];


      int gt_totals[4] = {0, 0, 0, 0}; // total num of A,C,G,Ts across all the libraries
      int total_reads = 0;
      std::vector<uint32_t> gtcounts;//(n_libs*4); // individual A,C,G, and T counts for each library
      for(size_t m = 0; m < mems.size(); m++) {
	member_ptr mem = mems[m];
	for(size_t l = 0; l < mem->libraries.size(); l++) {
	  std::vector<Base> reads = baseCall(mem->libraries[l], pos);
	  int gt_tmp[4] = {0, 0, 0, 0};
	  for(int r : reads) {
	    gt_tmp[r]++;
	    gt_totals[r]++;
	  }
	  gtcounts.push_back(gt_tmp[0]);
	  gtcounts.push_back(gt_tmp[1]);
	  gtcounts.push_back(gt_tmp[2]);
	  gtcounts.push_back(gt_tmp[3]);
	}
      }
      // allele order should start with the ref, and include only those nucleotides which have a read
      std::vector<int> allele_order_map;
      allele_order_map.push_back(ref);
      for(int a = 1; a < 4; a++) {
	int index = (ref + a)%4;
	if(gt_totals[index] > 0) {
	  //std::cout << "HERE at site " << pos << std::endl;
	  allele_order_map.push_back(index);
	}
      }
      if(allele_order_map.size() == 1)
	continue;

      hts::bcf::Variant rec = out.InitVariant();
      rec.target(chrom_.c_str());
      rec.position(pos);

      std::string allele_order_str;
      allele_order_str = Index2Allele(ref);
      for(int indx = 1; indx < allele_order_map.size(); indx++) {
	allele_order_str += std::string(",") + Index2Allele(allele_order_map[indx]);
      }
      rec.alleles(allele_order_str);

      std::vector<int32_t> allele_depths;
      for(int lindx = 0; lindx+3 < gtcounts.size(); lindx += 4) {
	for(int a : allele_order_map) {
	  int gt_loc = lindx + a;
	  allele_depths.push_back(gtcounts[gt_loc]);
	}
      }

      rec.samples("AD", allele_depths);
      out.WriteRecord(rec);
    }
  }
  

  

protected:
  static inline Base char2base(char c) {
    Base ret;
    switch(c) {
    case 'A':
    case 'a': ret = A; break;
    case 'C':
    case 'c': ret = C; break;
    case 'G':
    case 'g': ret = G; break;
    case 'T':
    case 't': ret = T; break;
    case 'N':
    case 'n': ret = N; break;
    default: ret = UNASSIGNED;
    }
    return ret;
  }

  // TODO: Change to Base
  char Index2Allele(int indx) {
    char ret;
    switch(indx) {
    case 0: ret = 'A'; break;
    case 1: ret = 'C'; break;
    case 2: ret = 'G'; break;
    case 3: ret = 'T'; break;
    }

    return ret;
  }
  


protected:
  pedigree_map pedigree; // Map of all members of the pedigree, indexed by id
  
  std::vector<Base> reference; // Reference to contig used to create pedigree sequences
  std::string chrom_; // name of chromosome/contig
  size_t start_pos_; // start location on contig
  int sam_seq_len_;

  int gametic_nodes_ = 0;

};


} // namespace sim
} // namespace dng




#endif /* SIMULATOR_H_ */

  
