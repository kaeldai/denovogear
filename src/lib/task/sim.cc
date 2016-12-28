/*
 * Copyright (c) 2014-2016 Reed A. Cartwright
 * Copyright (c) 2015-2016 Kael Dai
 * Authors:  Reed A. Cartwright <reed@cartwrig.ht>
 *           Kael Dai <kdai1@asu.edu>
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

#include <iostream>
#include <boost/algorithm/string.hpp>

#include <dng/task/sim.h>
#include <dng/sim/pedigree.h>
#include <dng/sim/models.h>
#include <dng/utility.h>

#include <htslib/faidx.h>
#include <htslib/hts.h>

using namespace sim;

typedef std::pair<std::string, std::size_t> loci; 

std::vector<loci> parse_loci(std::string &input_list) {
  std::vector<loci> ret_list;
  std::vector<std::string> loci_list;
  if(!input_list.empty()) {
    boost::split(loci_list, input_list, boost::is_any_of(","));
    for(std::string &str : loci_list) {
      std::size_t cpos = str.rfind(":");
      if(cpos == 0 || cpos == str.size()-1) {
	throw std::runtime_error("Invalid mutation position \"" + str +
				 "\". Must be in ID[:loc] format.");
      }

      
      if(cpos != std::string::npos) {
	std::string id = str.substr(0, cpos);
	char *endptr;
	std::size_t pos = std::strtoul(str.substr(cpos+1).c_str(), &endptr, 0);
	if(*endptr != 0 || str[cpos+1] == '-') {
	  throw std::runtime_error("Invalid mutation position \"" + str +
				   "\". loci location must be a non-negative integer.");
	}
	ret_list.emplace_back(id, pos);
      }
      else {
	ret_list.emplace_back(str, 0);
      }
    }    
  }

  return ret_list;
}


int dng::task::Sim::operator()(Sim::argument_type &arg) {

  // Must have specified a pedigree and region
  if(arg.ped.empty()) {
    throw std::runtime_error("No pedigree file.");
  }

  if(arg.region.empty()) {
    throw std::runtime_error("No region specified.");
  }

  // founder NT freqs.
  std::array<double, 4> freqs = dng::utility::parse_nuc_freqs(arg.nuc_freqs);
  
  gamma_t model_a(arg.gamma[0]);
  gamma_t model_b(arg.gamma[0]);
  

  
  // read the pedigree file
  ped::parse_pedigree(arg.ped);
  std::vector<ped::Member> family = ped::family;


  std::cout << "> Reading reference." << std::endl;
    
  // get information about contig
  size_t pos_begin = 0;
  size_t pos_end = 0;
  std::string contig_name;
  parse_region(arg.region, contig_name, pos_begin, pos_end);

  std::string base_fname = arg.prefix;
  if(base_fname.empty()) {
    base_fname = ped::family_name;
  }
  
  // read in reference contig
  int contig_len = 0;
  const faidx_t *contig_ref = fai_load(arg.input[0].c_str());
  char *contig_str = fai_fetch(contig_ref, arg.region.c_str(), &contig_len);

  // Create a reference
  std::vector<Base> reference(contig_len);  
  for(size_t site = 0; site < contig_len; ++site) {
    reference[site] = char2base(contig_str[site]);
  }

  // Check the values of --germline-mutation
  std::vector<std::vector<std::size_t>> germline_mutations(family.size());
  for(std::pair<std::string, std::size_t> &loci : parse_loci(arg.germline_mutations)) {
    std::string id = loci.first;
    std::size_t pos = loci.second;
  
    // make sure loci position is valid
    if(loci.second >= contig_len) {
      throw std::runtime_error("Germline mutation position \"" + id + ":" + std::to_string(pos) +
			       "\" is out of range. Please specify location between [0," + std::to_string(contig_len) + ")");
    }

    // make sure member is in the family
    
    std::size_t mem_loc = 0;
    for( ; mem_loc < family.size(); ++mem_loc) {
      if(id == family[mem_loc].id) {
	// TODO: May want to allow germline mutation on founder by forcing loci to be different from reference?
	if(family[mem_loc].is_founder) {
	  throw std::runtime_error("Invalid germline mutation \"" + id + ":" + std::to_string(pos) +
				   "\". Pedigree member " + id + " is a founder.");	  
	}	
	break;
      }
    }
    if(mem_loc == family.size()) {
      throw std::runtime_error("Invalid germline mutation \"" + id + ":" + std::to_string(pos) +
			       "\". No pedigree member " + id + ".");
    }
    else {
      germline_mutations[mem_loc].push_back(pos);
    }
  }

    // Check the values of --somatic-mutation
  std::vector<std::vector<std::size_t>> somatic_mutations(family.size());
  for(std::pair<std::string, std::size_t> &loci : parse_loci(arg.somatic_mutations)) {
    std::string id = loci.first;
    std::size_t pos = loci.second;
  
    // make sure loci position is valid
    if(loci.second >= contig_len) {
      throw std::runtime_error("Somatic mutation position \"" + id + ":" + std::to_string(pos) +
			       "\" is out of range. Please specify location between [0," + std::to_string(contig_len) + ")");
    }

    // make sure member is in the family    
    std::size_t mem_loc = 0;
    for( ; mem_loc < family.size(); ++mem_loc) {
      if(id == family[mem_loc].id) {
	break;
      }
    }
    if(mem_loc == family.size()) {
      throw std::runtime_error("Invalid somatic mutation \"" + id + ":" + std::to_string(pos) +
			       "\". No pedigree member " + id + ".");
    }
    else {
      somatic_mutations[mem_loc].push_back(pos);
    }
  }


  std::vector<std::vector<Genotype>> gametic_dna(family.size());
  std::vector<std::vector<Genotype>> somatic_dna(family.size());
	     
  model::DNGModel model(arg.ref_weight, arg.theta, arg.nuc_freqs, arg.mu, arg.mu_somatic, arg.model_a, arg.model_b);
  
  // A model that will always force a germline mutation
  model::DNGModel1GM model_g1(arg.ref_weight, arg.theta, arg.nuc_freqs, arg.mu_somatic, arg.model_a, arg.model_b);
  
  // A model that will always force a somatic mutation
  model::DNGModel1SM model_s1(arg.ref_weight, arg.theta, arg.nuc_freqs, arg.mu, arg.model_a, arg.model_b);
  
  size_t seed = (arg.seed > 0 ? seed : time(0));
  model.set_seed(seed);
  model_g1.set_seed(seed);
  model_s1.set_seed(seed);



  /*
  std::cout << "> Creating germline DNA for founders." << std::endl;
  
  // Use the population priors to intialize founder dna
  size_t remaining = family.size();
  for(int mem = 0; mem < family.size() && family[mem].is_founder; ++mem) {
    ped::Member &m = family[mem]; // reference to member
    std::vector<Genotype> &mdna = gametic_dna[mem]; // reference to member's dna
    for(size_t site = 0; site < contig_len; site++) {
      Base ref = reference[site];
      Genotype gt = model.createGameticDNA(m, ref);
      mdna.push_back(gt);
      ped::collect_founder_stats(m, ref, gt);
    }
    m.dna_set = true;
    --remaining;
  }


  std::cout << "> Creating germline DNA for children." << std::endl;
  
  // Create gametic DNA for non-founders
  while(remaining) {
    for(size_t mem; mem < family.size(); ++mem) {
      ped::Member &m = family[mem];
      if(m.dna_set) {
	// member already has gametic DNA
	continue;
      }
      ped::Member &dad = get_dad(m);
      ped::Member &mom = get_mom(m);
      if(!dad.dna_set || !mom.dna_set) {
	// one of member's parent's doesn't have their dna, come back later
	continue;
      }

      // member doesn't have DNA, but both parents do. use a germline tranformation
      // to create new DNA for member.
      std::vector<Genotype> &cdna = gametic_dna[mem]; // child's dna
      std::vector<Genotype> &ddna = gametic_dna[m.dpos]; // dad's dna
      std::vector<Genotype> &mdna = gametic_dna[m.mpos]; // mom's dna
      for(size_t site = 0; site < contig_len; ++site) {
	Base ref = reference[site];
	Genotype dgt = ddna[site];
	Genotype mgt = mdna[site];
	Genotype cgt = model.getGermlineTransition(m, ref, mgt, dgt);
	cdna.push_back(cgt);
      }

      // Force user-defined mutations at specific locations.
      // NOTE: For a typically case with a large contig and only a few forced mutations. It should be faster to
      //       force mutation afterwards than having to check at every position.
      if(germline_mutations[mem].size() > 0) {
	for(size_t site : germline_mutations[mem]) {
	  Base ref = reference[site];
	  Genotype dgt = ddna[site];
	  Genotype mgt = mdna[site];
	  Genotype cgt = model_g1.getGermlineTransition(m, ref, mgt, dgt);
	  cdna[site] = cgt;
	}
      }    

      // collect statistics about non-founder germline DNA
      for(size_t site = 0; site < contig_len; ++site) {
	ped::collect_germline_stats(m, cdna[site], mdna[site], ddna[site]);
      }
      
      m.dna_set = true;
      --remaining;
    }
  }

  std::cout << "> Creating somatic DNA." << std::endl;
  
  // Create somatic mutation
  for(size_t mem = 0; mem < family.size(); ++mem) {
    ped::Member &m = family[mem];
    std::vector<Genotype> &dna = gametic_dna[mem];

    for(size_t site = 0; site < contig_len; ++site) {
      Base ref = reference[site];
      Genotype gametic_gt = dna[site];
      Genotype somatic_gt = model.getSomaticTransition(m, ref, gametic_gt);
      somatic_dna[mem].push_back(somatic_gt);
    }

    // Changes sites where user if forcing a mutation
    if(somatic_mutations[mem].size() > 0) {
      for(size_t site : somatic_mutations[mem]) {
	Base ref = reference[site];
	Genotype gametic_gt = dna[site];
	somatic_dna[mem][site] = model_s1.getSomaticTransition(m, ref, gametic_gt);
      }
    }

    // Collect statistics about somatic mutations
    for(size_t site = 0; site < contig_len; ++site) {
      ped::collect_somatic_stats(m, dna[site], somatic_dna[mem][site]);
    }
  }

  
  // Save germline and somatic genotypes in a vcf file.
  io::VCFOutput vcf_file((base_fname + ".vcf").c_str());
  vcf_file.addContig(contig_name.c_str(), contig_len);
  for(ped::Member &m : family) {
   vcf_file.addSample((std::string("GL-") + m.id).c_str());
   vcf_file.addSample((std::string("SL-") + m.id).c_str());
  }
  vcf_file.writeHeader();
  std::vector<Genotype> gts(family.size()*2);
  for(size_t site = 0; site < contig_len; ++site) {
    Base ref = reference[site];
    vcf_file.newSite(contig_name, (pos_begin + site), ref);
    for(size_t mem = 0; mem < family.size(); ++mem) {
      gts[2*mem] = gametic_dna[mem][site];
      gts[2*mem+1] = somatic_dna[mem][site];
    }
    vcf_file.addGenotypes(gts);
    vcf_file.writeSite();
  }
  vcf_file.close();
  
  
  std::cout << "> Creating read calls." << std::endl;
  
  // Create library reads and output to file file
  io::TadOutput tad_file((base_fname + ".tad").c_str());
  tad_file.addContig(contig_name, contig_len);
  std::vector<std::string> lb_names;
  for(ped::Member &m : family) {
    lb_names.emplace_back(m.id);
    tad_file.addLibrary(m.id);
  }
  tad_file.writeHeader();

  
  std::vector<reads_list> read_depths(family.size());
  for(size_t site = 0; site < contig_len; ++site) {
    Base ref = reference[site];
    size_t pos = pos_begin + site + 1;
    tad_file.newSite(contig_name, pos, ref);
    for(size_t mem = 0; mem < family.size(); ++mem) {
      model.call(arg.read_depths, ref, somatic_dna[mem][site], read_depths[mem]); 
    }
    tad_file.addCalls(read_depths); 
    tad_file.writeSite();
  }
  tad_file.close();

  io::CSVOutput statsfile((base_fname + ".csv").c_str());
  statsfile.addStats(family);
  statsfile.close();

  */

  return EXIT_SUCCESS;
}
