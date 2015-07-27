#ifndef DNG_SIM_DATA_FACTORY_H
#define DNG_SIM_DATA_FACTORY_H

#include "htslib/sam.h"
#include "dng/hts/bcf.h"
#include <vector>
#include <unordered_map>
#include <set>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <dng/io/ped.h>
#include <boost/format.hpp>
#include <random>
#include <math.h>



#define ID_UNKNOWN -1
#define IS_UNKNOWN(mem) ((mem)->mid == ID_UNKNOWN)
#define IS_EMPTY(mem) (mem == nullptr)

#define PED_NOPARENT "0"


namespace dng {
  namespace sim {
    
    typedef size_t member_id;
    typedef size_t family_id;
    typedef dng::io::Pedigree::Gender Gender;

    typedef std::array<double, 16> genotype_dist;

    // Format of output data file. 'None' indicates only PED file is generated.
    enum DataFormat { VCF, BCF, SAM, BAM, CRAM, None };

    enum Base : uint8_t { A, C, G, T, N, REF, UNASSIGNED };

    static char nt2char[] = {'A', 'C', 'G', 'T', 'N', 'R', 'U'};


    // The number of data files generated.
    // TODO: Currenlty only one large file implemented
    enum FileScheme { 
      ONE_FILE, // One file containing all samples and libraries
      ONE_PER_SAMPLE, // Each sample gets its own library
      ONE_PER_LIB // Each sample:library pair gets a seperate file
    };


    typedef std::unordered_map<size_t, Base> reference_map;

    struct SampleLibrary {
      std::string name;
      size_t depth;
      reference_map dna[2];
      
      //SampleLibrary(std::string &name_, size_t depth_) : name(name_), depth 

    };

    struct Member {
      member_id mid; // interal usage
      std::string name; // name that shows up in pedigree and data file
      family_id fid;
      Gender sex;
      Member *mom; 
      Member *dad;
      bool hasDNA; // Used when generating the pedigree contigs

      // Used to keep track of differences between member's DNA and the reference. We need one for each chromatid and somatic vs gametic
      reference_map gametes[2];
      std::vector<SampleLibrary> libraries;
      
      Member(member_id mid_, family_id fid_, Gender sex_) : mid(mid_), fid(fid_), sex(sex_) {
	name = (boost::format("NA%04d") % mid).str();
	mom = nullptr;
	dad = nullptr;
	hasDNA = false;
      }
      
      std::string family_name() {
	return (boost::format("F%03d") % fid).str();
      }

      std::string dad_id() {
	if(dad == nullptr)
	  return 0;
	else 
	  return dad->name;
      }

      std::string mom_id() {
	if(mom == nullptr)
	  return std::string("0");
	else 
	  return mom->name;
      }

      /**
       * Change the nucleotide at 'site' for chromosome pair number 'chrom_num' to 'base'.
       */
      void updateDNA(int chrom_num, size_t site, Base base) {
	if(chrom_num >= 2)
	  throw std::runtime_error("Attempting to set base from more than diploidy DNA.");
	
	gametes[chrom_num].insert({site, base});
      }

      /**
       * Makes copies of mom and dad's dna and saves into member. Selection of chromosome pair is random
       */
      void inheritDNA(Member *mom, Member *dad) {
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
      Base getNuc(size_t chromatid, size_t site) {
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

      void addLibrary(std::string &name, size_t depth) {
	// TODO: Check that library with same name doesn't already exists

	SampleLibrary lib;
	lib.name = std::string(name);
	lib.depth = depth;
	libraries.push_back(lib);
      }

      Base getLibDNA(size_t l_indx, size_t chrom, size_t site) {
	reference_map &rmap = libraries[l_indx].dna[chrom];
	reference_map::const_iterator nt = rmap.find(site);
	if(nt ==  rmap.end()) {
	  return REF;
	}
	else {
	  return nt->second;
	}
      }

      void updateSample(size_t l_indx, size_t chrom, size_t site, Base base) {
	libraries[l_indx].dna[chrom].insert({site, base});
      }
    };

 
    class Factory {
    public:
      // static Factory *newInstance(Model=...)
      // AddTrio()
      // buildPedigree(stream);
      // buildData(vector<streams>, Type=BCF/VCF/BAM/SAM)  
      // getLB(), getSM(), getGroup() will be different depending on if BCF or BAM format

      static Factory *newInstance(DataFormat df = DataFormat::None) {
	return new Factory(df);
      }

      /*
      static SAMBuilder *newSAMBuilder() {
	SAMBuilder *sam = new SAMBuilder();
	return sam;
      }
      */

      Member *AddTrio(Member *child, Member *mom, Member *dad) {
	return AddTrio(child, (child == nullptr ? Gender::Unknown : child->sex), mom, dad);
      }

      Member *AddTrio(Member *child, Gender sex, Member *mom, Member *dad) {
	// Make sure the trio includes a child
	if(!IS_EMPTY(child) && IS_UNKNOWN(child)) {
	  // TODO: Attempting to add a childless mom and dad to pedigree. Should this be valid? Currently just skipping
	  return nullptr;
	}

	// TODO: Check that if either child, mom or dad are not NULL, then their genders are as to be expected 


	// Determine family units
	std::set<family_id> families;
	for(Member *m : {child, mom, dad}) {
	  if(!IS_EMPTY(m) && m->fid != 0)
	    families.insert(m->fid);
	}
	
	family_id family;
	if(families.empty()) {
	  // No one in the trio belongs to another another family unit, create a new fam
	  family = get_famid();
	}
	else if(families.size() == 1) {
	  // Only one of the members of the trio belong to another family, add the other two
	  family = *families.begin();
	}
	else {
	  // Merging of two different families by one trio
	  std::set<family_id>::iterator itr = families.begin();
	  for(itr++; itr != families.end(); itr++) {
	    //merge_familes(families[a], family);
	  }
	}
	  
	// If either mom, dad, or child are not specified create new entries 
	if(IS_EMPTY(mom)){
	  mom = new Member(get_uid(), family, Gender::Female);
	  members.push_back(mom);
	}
	if(IS_EMPTY(dad)){
	  dad = new Member(get_uid(), family, Gender::Male);
	  members.push_back(dad);
	}
	if(IS_EMPTY(child)) {
	  child = new Member(get_uid(), family, sex);
	  members.push_back(child);
	}
	
	child->mom = (IS_UNKNOWN(mom) ? nullptr : mom);
	child->dad = (IS_UNKNOWN(dad) ? nullptr : dad);

	return child;           
      }

      void AddLibrary(Member *m, std::string &libname, size_t depth) {


      }


      template<typename Stream>
      void publishPed(Stream &output) {
	
	for(int a = 0; a < members.size(); a++) {
	  Member *mem = members[a];
	  
	  // Create a unique family name
	  //char family_name[5];
	  //sprintf(family_name, "F%03d", mem->fam);  

	  // Get the names of the mom and dad, 
	  std::string mom = (mem->mom == nullptr ? PED_NOPARENT : mem->mom->name);
	  std::string dad = (mem->dad == nullptr ? PED_NOPARENT : mem->dad->name);
	  
	  output << mem->family_name() << "\t"
		 << mem->name << "\t"
		 << dad << "\t"
		 << mom << "\t"
		 << get_gender(mem->sex) << "\t"
		 << (mem->family_name() + mem->name)  << std::endl; // Not exactly sure what's this is for, but parsed out by ped.h! 
	  
	}
      }

      
      
      //genotype_dist genotype_probs(char ref, double theta, std::array<double, 4> nuc_freq, std::array<double,4> prior) {
      genotype_dist genotype_weights(double theta, std::array<double, 4> nuc_freq, std::array<double,4> prior) {
	std::array<double, 16> weights;
	//Eigen::Array<double, 16, 1> weights;
	//dng::GenotypeArray ret{10};
	//double weights[16];
	double alpha[4] = {
	  theta *nuc_freq[0] + prior[0], theta *nuc_freq[1] + prior[1],
	  theta *nuc_freq[2] + prior[2], theta *nuc_freq[3] + prior[3]
	};
	double alpha_sum = alpha[0] + alpha[1] + alpha[2] + alpha[3];
	double alpha_sum2 = (alpha_sum + alpha_sum*alpha_sum);

	// chromatid pattern will matter in heritence, so when we randomly assign a parent a heterozygous site A/T != T/A 
	weights[0]  = (alpha[0] + alpha[0]*alpha[0]) / alpha_sum2; // AA
	weights[1]  = alpha[0]*(alpha[1]) / alpha_sum2; // AC
	weights[2]  = alpha[0]*(alpha[2]) / alpha_sum2; // AG
	weights[3]  = alpha[0]*(alpha[3]) / alpha_sum2; // AT

	weights[4]  = weights[1]; // CA
	weights[5]  = (alpha[1] + alpha[1]*alpha[1]) / alpha_sum2; // CC
	weights[6]  = alpha[1]*(alpha[2]) / alpha_sum2; // CG
	weights[7]  = alpha[1]*(alpha[3]) / alpha_sum2; // CT

	weights[8]  = weights[2]; // GA 
	weights[9]  = weights[6]; // GC 
	weights[10] = (alpha[2] + alpha[2]*alpha[2]) / alpha_sum2; // GG
	weights[11] = alpha[2]*(alpha[3]) / alpha_sum2; // GT 

	weights[12] = weights[3]; // TA 
	weights[13] = weights[7]; // TC 
	weights[14] = weights[11]; // TG
	weights[15] = (alpha[3] + alpha[3]*alpha[3]) / alpha_sum2; // TT
	
	return weights;
      }

      int Allele2Index(char c) {
	int ret = 4;
	switch(c) {
	case 'a':
	case 'A': ret = 0; break;
	case 'c':
	case 'C': ret = 1; break;
	case 'g':
	case 'G': ret = 2; break;
	case 't':
	case 'T': ret = 3; break;
	}
	return ret;
      }
      
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

      
      std::array<std::pair<Base, Base>, 16> genotypes = {{{A, A}, {A, C}, {A, G}, {A, T},
							  {C, A}, {C, C}, {C, G}, {C, T},
							  {G, A}, {G, C}, {G, G}, {G, T},
							  {T, A}, {T, C}, {T, G}, {T, T}}};


	/*
      std::pair<Base, Base> index2Genotype(int index) {

	std::vector<std::pair<char, char>> genotypes = {{'A','A'}, {'A','C'}, {'A','G'}, {'A','T'},
							{'C','A'}, {'C','C'}, {'C','G'}, {'C','T'},
							{'G','A'}, {'G','C'}, {'G','G'}, {'G','T'},
							{'T','A'}, {'T','C'}, {'T','G'}, {'T','T'}};
	return genotypes[index];
      }
	*/

      void createPriorsDist() {
	genotype_dist ref_weights[] = { genotype_weights(theta_, nuc_freqs_, {ref_weight_, 0, 0, 0}),
					genotype_weights(theta_, nuc_freqs_, {0, ref_weight_, 0, 0}),
					genotype_weights(theta_, nuc_freqs_, {0, 0, ref_weight_, 0}),
					genotype_weights(theta_, nuc_freqs_, {0, 0, 0, ref_weight_})};
	double interval[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};

	pop_priors_dists.emplace_back(std::begin(interval), std::end(interval), ref_weights[0].begin());
	pop_priors_dists.emplace_back(std::begin(interval), std::end(interval), ref_weights[1].begin());
	pop_priors_dists.emplace_back(std::begin(interval), std::end(interval), ref_weights[2].begin());
	pop_priors_dists.emplace_back(std::begin(interval), std::end(interval), ref_weights[3].begin());
      }

      /*
      void createFounderGenotypes() {
	genotype_dist ref_weights[] = { genotype_weights(theta_, nuc_freqs_, {ref_weight_, 0, 0, 0}),
					genotype_weights(theta_, nuc_freqs_, {0, ref_weight_, 0, 0}),
					genotype_weights(theta_, nuc_freqs_, {0, 0, ref_weight_, 0}),
					genotype_weights(theta_, nuc_freqs_, {0, 0, 0, ref_weight_})};

	std::vector<std::piecewise_constant_distribution<>> dists;
	double interval[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
	dists.emplace_back(std::begin(interval), std::end(interval), ref_weights[0].begin());
	dists.emplace_back(std::begin(interval), std::end(interval), ref_weights[1].begin());
	dists.emplace_back(std::begin(interval), std::end(interval), ref_weights[2].begin());
	dists.emplace_back(std::begin(interval), std::end(interval), ref_weights[3].begin());

	std::cout << "           ";
	for(int a = 0; a*10 < reference.size(); a++) {
	  std::cout << a << "        ";
	}
	std::cout << std::endl;
	//std::cout << "reference: " << reference << std::endl;
	std::random_device rd;
	std::mt19937 gen(rd());
	for(int a = 0; a < members.size(); a++) {
	  Member *m = members[a];
	  if(m->mom == NULL && m->dad == NULL) {
	    for(int b = 0; b < reference.size(); b++) {
	      char ref = reference[b];
	      int ref_idx = Allele2Index(ref);
	      std::pair<Base, Base> genotype = genotypes[floor(dists[ref_idx](gen))];//index2Genotype(floor(dists[ref_idx](gen)));
	      if(genotype.first != ref) {
		m->updateDNA(1, b, genotype.first);
	      }
	      if(genotype.second != ref) {
		m->updateDNA(2, b, genotype.second);
	      }
	    }
	    m->hasDNA = true;
	  }
	  
	  // PRINTING
	  if(m->mom == NULL && m->dad == NULL)
	  {
	    std::cout << m->name << "   : ";
	    for(size_t a = 0; a < reference.size(); a++) {
	      std::cout << m->getNuc(0, a);
	    }
	    std::cout << std::endl;
	    
	    std::cout << "           ";
	    for(size_t a = 0; a < reference.size(); a++) {
	      std::cout << m->getNuc(1, a);
	    }
	    std::cout << std::endl;
	  }

	}
      }
      */

      Base transitions[4] = {G, A, T, C};
      Base transversions[4][2] = {{C, T}, {C, T}, {A, G}, {A, G}};
      
      /*
      inline Base transition(Base nuc){ 
	Base mutation;
	switch(nuc) {
	case A: mutation = G; break;
	case G: mutation = A; break;
	case C: mutation = T; break;
	case T: mutation = C; break;
	default:mutation = N;
	}
	return mutation;
      }
      
      inline char transversion(char nuc){
	char options[2][2] = {{'C', 'T'}, {'A','G'}};
	int r = rand() % 2;
	if(nuc == 'A' || nuc == 'G') {
	  return options[0][r];
	}
	else {
	  return options[1][r];
	  }
      }
      */

      /*
      void createChildrenDNA() {
	
	double interval[] = {0, 1, 2, 3};
	double probs[] =  {transitons_mut_, transversion_mut_, (1.0 - transitons_mut_ - transversion_mut_)};
	std::piecewise_constant_distribution<> dist(std::begin(interval), std::end(interval), std::begin(probs));
	std::random_device rd;
	std::mt19937 gen(rd());
	int remaining_members = 0;
	do {
	  for(int a = 0; a < members.size(); a++) {
	    Member *m = members[a];
	    if(m->hasDNA == true)
	      continue;
	    else if(!m->dad->hasDNA || !m->mom->hasDNA) {
	      remaining_members++;
	      continue;
	    }
	    else {
	      m->inheritDNA(m->mom, m->dad);
	      for(int b = 0; b < reference.size(); b++) {
		int mut_type = floor(dist(gen));

		if(mut_type == 2) {
		  continue;
		}
		else if(mut_type == 0) {
		  Base oldbase = m->getNuc(1, b);
		  if(oldbase == ' ')
		    oldbase = reference[b];
		  char newbase = transition(oldbase);
		  m->updateDNA(1, b, newbase);

		}
		else if(mut_type == 1) {
		  char oldbase = m->getNuc(1, b);
		  if(oldbase == ' ')
		    oldbase = reference[b];
		  char newbase = transversion(oldbase);
		  m->updateDNA(1, b, newbase);
		  //std::cout << "Transversion at " << b << std::endl;
		}

	      }
	      
	    }
	  }
	} while(remaining_members > 0);
      }
      */

      char mutation_options[4][3] = {{C, G, T}, {A, G, T}, {A, C, T}, {A, C, G}};

      void createLibraryMutations() {
	double interval[] = {0, 1, 2, 3};
	double probs[] =  {transitons_mut_, transversion_mut_, (1.0 - transitons_mut_ - transversion_mut_)};
	std::piecewise_constant_distribution<> dist(std::begin(interval), std::end(interval), std::begin(probs));
	std::random_device rd;
	std::mt19937 gen(rd());

	for(Member *m : members) {
	  // If user hasn't specified any libraries create a default one
	  if(m->libraries.size() == 0) {
	    std::string libname = std::string(m->name) + "-" + m->family_name();
	    m->addLibrary(libname, 10);
	  }

	  for(size_t l_indx = 0; l_indx < m->libraries.size(); l_indx++) {
	    for(size_t site = 0; site < reference.size(); site++) {
	      for(size_t chrom : {0, 1}) {
		int mut_type = floor(dist(gen));
		if(mut_type == 2) {
		  // No Mutation
		  continue;
		}
		else if(mut_type == 0) {
		  // Transition mutation
		  Base oldbase = m->getLibDNA(l_indx, chrom, site);
		  if(oldbase == REF)
		    oldbase = reference[site];
		  m->updateSample(l_indx, chrom, site, transitions[oldbase]);
		}
		else {
		  // Transversion
		  Base oldbase = m->getLibDNA(l_indx, chrom, site);
		  if(oldbase == REF)
		    oldbase = reference[site];
		  
		  int r = rand() % 2; // 2 possible choices for each transverison
		  m->updateSample(l_indx, chrom, site, transversions[oldbase][r]);
		  
		}
	      }
	    }
	  }
	}
      }


      void publishData() {
	//////////////////////////
	std::cout << "ref:   ";
	for(Base b : reference) {
	  std::cout << nt2char[b];
	}
	std::cout << std::endl;
      	//////////////////////////

	createPriorsDist();
	initGameteDNA();
	
	//////////////////////////
	for(Member *m : members) {
	  std::cout << m->name << " ";
	    for(size_t chrom : {0, 1}) {
	      for(int site = 0; site < reference.size(); site++) {
		Base base = m->getNuc(chrom, site);
		Base ref = reference[site];
		if(base == REF)
		  std::cout << ' ';
		else   
		  std::cout << nt2char[base];
	      }
	      std::cout << std::endl;
	      std::cout << "       ";
	    }
	  std::cout << std::endl;
	}
	//////////////////////////

	
	
	createLibraryMutations();
	//publishVCFData();

      }

      void createFounderDNA(Member *mem) {
	// Go through each site in the reference contig and randomly select a genotype based on the reference. 
	for(size_t site = 0; site < reference.size(); site++) {
	  Base ref = reference[site];
	  std::pair<Base, Base> gt = genotypes[floor(pop_priors_dists[ref](ran_generator))];

	  // Only need to keep track if site is different from the reference.
	  if(gt.first != ref) {
	    mem->updateDNA(0, site, gt.first);
	  }
	  if(gt.second != ref) {
	    mem->updateDNA(1, site, gt.second);
	  }
	}	    
      }

      void createChildDNA(Member *child, Member *parent) {
	// Temporarly create a parent
	Member *tmp_parent = new Member(-1, -1, Gender::Unknown);
	createFounderDNA(tmp_parent);

	// update child
	createChildDNA(child, tmp_parent, parent);

	// delete parent
	delete tmp_parent;
      }

      void createChildDNA(Member *child, Member *mom, Member *dad) {
	if(child == NULL || mom == NULL || dad == NULL)
	  // TODO: throw error
	  return;

	// Weighted probs of type of mutation; transtion, transversion, or none.
	double interval[] = {0, 1, 2, 3};
	double probs[] =  {transitons_mut_, transversion_mut_, (1.0 - transitons_mut_ - transversion_mut_)};
	std::piecewise_constant_distribution<> dist(std::begin(interval), std::end(interval), std::begin(probs));
	std::random_device rd;
	std::mt19937 gen(rd());

	child->inheritDNA(mom, dad);
	// Go through the DNA and 
	// TODO: Figure out a way to sample a few sites rather than iterate through the entire contig
	for(int site = 0; site < reference.size(); site++) {
	  for(size_t chrom : {0, 1}) {
	    int mut_type = floor(dist(gen));
	    if(mut_type == 2) {
	      // No Mutation
	      continue;
	    }
	    else if(mut_type == 0) {
	      // Transition mutation
	      Base oldbase = child->getNuc(chrom, site);
	      if(oldbase == REF)
		oldbase = reference[site];
	      child->updateDNA(chrom, site, transitions[oldbase]);
	    }
	    else {
	      // Transversion
	      Base oldbase = child->getNuc(chrom, site);
	      if(oldbase == REF)
		oldbase = reference[site];

	      int r = rand() % 2; // 2 possible choices for each transverison
	      child->updateDNA(chrom, site, transversions[oldbase][r]);
	      
	    }
	  }
	}
      }

      /**
       * Iterate through all the members of the pedigree creating DNA based on the reference and parential DNA.
       */
      void initGameteDNA() {
	// loop through all the members, since we don't assume anything about the order of the pedigree list we have to 
	// make mulitple loops through; first initializing the pedigree founders, then the founders' children, and then
	// their children, etc.
	size_t remaining = members.size();
	while(remaining > 0) {
	  for(Member *m : members) {
	    if(m->hasDNA) {
	      // Member already has gametic DNA, skip
	      continue;
	    }
	    else if(m->mom == NULL && m->dad == NULL) {
	      // Member is a pedigree founder - does not have parents
	      createFounderDNA(m);
	      m->hasDNA = true;
	      remaining--;
	    }
	    else if(!(m->mom == NULL) != !(m->dad == NULL)) {
	      // Special case where member has only one parent in pedigree and another is missing. Need to wait to parent has
	      // DNA filled.
	      if(m->mom != NULL && m->mom->hasDNA) {
		createChildDNA(m, m->mom);
		m->hasDNA = true;
		remaining--;
	      }
	      else if(m->dad != NULL && m->dad->hasDNA) {
		createChildDNA(m, m->dad);
		m->hasDNA = true;
		remaining--;
	      }
	      else {
		continue;
	      }
	    } 
	    else {
	      // A child with both parents in pedigree. Wait until parents have their DNA before initializing.
	      if(m->mom->hasDNA && m->dad->hasDNA) {
		createChildDNA(m, m->mom, m->dad);
		m->hasDNA = true;
		remaining--;
	      }
	    }	    
	  }
	}
      }



      std::array<size_t, 4> depth_count(Member *m, size_t pos) {
	char ref = reference[pos];
	char allele1 = m->getNuc(0, pos);
	if(allele1 == ' ')
	  allele1 = ref;

	char allele2 = m->getNuc(1, pos);
	if(allele2 == ' ')
	  allele2 = ref;

	double interval[] = {0, 1, 2, 3, 4};
	std::array<double, 4> probs;
	if(allele1 == allele2) {
	  double e = (1.0-homozygote_match_)/3.0;
	  probs = {e, e, e, e};
	  probs[Allele2Index(allele1)] = homozygote_match_;
	}
	else {
	  double e = (1.0-heterozygote_match_)/2.0;
	  probs = {e, e, e, e};
	  probs[Allele2Index(allele1)] = homozygote_match_/2.0;
	  probs[Allele2Index(allele2)] = homozygote_match_/2.0;
	}

	std::piecewise_constant_distribution<> dist(std::begin(interval), std::end(interval), std::begin(probs));
	std::random_device rd;
	std::mt19937 gen(rd());

	
	std::array<size_t, 4> ret = {0, 0, 0, 0};
	for(int a = 0; a < depth; a++) {
	  ret[floor(dist(gen))]++;
	}

	/*
	std::cout << "-----" << std::endl;
	std::cout << "depth = " << depth << std::endl;
	for(int a : ret) {
	  std::cout << a << ",";
	}
	std::cout << std::endl << "-----" << std::endl;
	*/

	return ret;
      }


      void publishVCFData(std::string &file) {
	hts::bcf::File out(file.c_str(), "w");
	out.AddHeaderMetadata("##fileformat=VCFv4.2");
	out.AddHeaderMetadata("##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">");
	out.AddContig("1", 243199373);
	for(int a = 0; a < members.size(); a++) {
	  Member *mem = members[a];
	  std::string sample = (boost::format("%s:LB-%s") % mem->name % mem->family_name()).str();
	  out.AddSample(sample.c_str());
	}
	out.WriteHeader();

	for(size_t pos  = 0; pos < reference.size(); pos++) {
	  //hts::bcf::Variant rec = out.InitVariant();
	  //rec.target("1");
	  //rec.position(pos);
	  char ref = reference[pos];
	  int ref_indx = Allele2Index(ref);
	  int gt_totals[4] = {0, 0, 0, 0};
	  int total_reads = 0;
	  std::vector<int> gtcounts(members.size()*4);
	  for(size_t m = 0; m < members.size(); m++) {
	    Member *mem = members[m];
	    std::array<size_t, 4> counts = depth_count(mem, pos);
	    for(int a : {0, 1, 2, 3}) {
	      gtcounts[m*4+a] = counts[a];
	      gt_totals[a] += counts[a];
	      total_reads += counts[a];
	    }
	  }

	  /*
	  std::cout << ">> ";
	  for(int a : gtcounts)
	    std::cout << a << ",";
	  std::cout << std::endl;
	  */

	  std::vector<int> allele_order_map;
	  allele_order_map.push_back(ref_indx);
	  for(int a = 1; a < 4; a++) {
	    int index = (ref_indx + a)%4;
	    if(gt_totals[index] > 0) {
	      //std::cout << "HERE at " << index << " having " << gt_totals[index] << "reads." << std::endl;
	      allele_order_map.push_back(index);
	    }
	  }

	  if(allele_order_map.size() == 1)
	    continue;


	  hts::bcf::Variant rec = out.InitVariant();
	  rec.target(chrom.c_str());
	  rec.position(pos);

	  std::string allele_order_str;
	  allele_order_str = ref;
	  for(int indx = 1; indx < allele_order_map.size(); indx++) {
	    allele_order_str += std::string(",") + Index2Allele(allele_order_map[indx]);
	  }

	  rec.info("LL", static_cast<float>(1.0));

	  rec.alleles(allele_order_str);

	  //std::cout << "HERE" << std::endl;

	  
	  std::vector<int32_t> allele_depths;
	  for(int m_indx = 0; m_indx < members.size(); m_indx++) {
	    for(int a : allele_order_map) {
	      int ad_indx = m_indx*4 + a;//allele_order_map[a];
	      allele_depths.push_back(gtcounts[ad_indx]);
	    }
	  }

	  rec.samples("AD", allele_depths);
	  out.WriteRecord(rec);
	}
      }

      std::string getContig(Member *m) {
	std::string contig;
	// randomly choose a contig to sample
	int chr = (rand() % 2);
	
	for(int pos = 0; pos < reference.size(); pos++) {
	  Base ref = reference[pos];
	  Base allele = m->getNuc(chr, pos);
	  if(allele == REF) 
	    allele = ref;

	  // need to add random mutation
	  
	  contig += nt2char[allele];  
	}
	return contig;
      }

      void publishSAMData(std::string &filename) {
	//std::cout << "----HERE-----" << std::endl;
	std::stringstream hdr_txt;
	//filename += ".sam";
	hdr_txt << "@HD\tVN:0.1\tSO:unknown\tGO:none" << std::endl;
	hdr_txt << "@SQ\tSN:1\tLN:249250621" << std::endl;
	for(int a = 0; a < members.size(); a++) {
	  Member *mem = members[a];
	  if(mem->libraries.size() == 0) {
	    std::string id = std::to_string(mem->mid);
	    //char name[8];
	    //sprintf(name, "NA%05d", id);
	    hdr_txt << "@RG\t"
		    << "ID:" << id << "\t"
		    << "LB:" << mem->name + "-" + mem->family_name() << "\t" //(std::string(name)+"-"+ring(id)) << "\t"
		    << "SM:" << mem->name << std::endl;
	  }
	  else {

	  }
	}
	std::string header_str = hdr_txt.str();

	bam_hdr_t *hdr = bam_hdr_init();
	hdr->text = (char *)malloc(sizeof(char)*(header_str.size()+1));
	strcpy(hdr->text, header_str.c_str());

	std::stringstream sam_file;
	for(int mindx = 0; mindx < members.size(); mindx++) {
	  for(int i = 0; i < depth; i++) {
	    Member *mem = members[mindx];
	    sam_file << chrom << "\t";
	    sam_file << 0 << "\t";
	    sam_file << 1 << "\t";
	    sam_file << 0 << "\t";
	    sam_file << 0 << "\t";
	    sam_file << reference.size() << "M" << "\t";
	    sam_file << "=" << "\t";
	    sam_file << 0 << "\t";
	    sam_file << 0 << "\t";
	    sam_file << getContig(mem) << "\t";
	    sam_file << "*" << "\t";
	    sam_file << "RG:Z:" << std::to_string(mem->mid);
	    sam_file << std::endl;
	    
	    //std::cout << getMContig(mem) << std::endl;

	    //bam1_t *rec = bam_init1();
	    //rec->core.tid = 20;
	    //rec->core.pos = 0;
	    
	  }
	}
	
	std::cout << hdr_txt.str();
	std::cout << sam_file.str();

	//sam_hdr_parse(header_str.size(), header_str.c_str());
	
	samFile *fp = sam_open(filename.c_str(), "w");
	//std::cout << "> " << hdr->text << std::endl;
	
	//bam1_t *line = bam_init1(); 
	//sam_write1(fp, hdr, line); 
	sam_hdr_write(fp, hdr);
	sam_close(fp);
	
      }


    private:
      Factory(DataFormat df) : dataformat_(df){ 
	std::string reference_str = "TTAATAGGGCGTTGCTGGCGGGCGTTGGGTGTGGCCCGCAGTCCTGGTTGAGGATTGCCC";
	setReference(reference_str);
      };

      int get_uid() { return id_list++; };

      family_id get_famid() { 
	family_id id = fam_list++;
	fam_rels[id] = std::vector<Member*>();
	return id;
      };

      
      static Base char2base(char c) {
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
      

      void setReference(std::string &fasta, std::string &range) {

      }

      /*
      static const Base char2base[256] = {
	N, N, N, N,  N, N, N, N,  N, N, N, N,  N, N, N, N,
	N, N, N, N,  N, N, N, N,  N, N, N, N,  N, N, N, N,
	N, N, N, N,  N, N, N, N,  N, N, N, N,  N, N, N, N,
	N, N, N, N,  N, N, N, N,  N, N, N, N,  N, N, N, N,

	N, A, N, C,  N, N, N, G,  N, N, N, N,  N, N, N, N,
	N, N, N, N,  T, N, N, N,  N, N, N, N,  N, N, N, N,
	N, A, N, C,  N, N, N, G,  N, N, N, N,  N, N, N, N,
	N, N, N, N,  T, N, N, N,  N, N, N, N,  N, N, N, N,
	
	N, N, N, N,  N, N, N, N,  N, N, N, N,  N, N, N, N,
	N, N, N, N,  N, N, N, N,  N, N, N, N,  N, N, N, N,
	N, N, N, N,  N, N, N, N,  N, N, N, N,  N, N, N, N,
	N, N, N, N,  N, N, N, N,  N, N, N, N,  N, N, N, N,
	
	N, N, N, N,  N, N, N, N,  N, N, N, N,  N, N, N, N,
	N, N, N, N,  N, N, N, N,  N, N, N, N,  N, N, N, N,
	N, N, N, N,  N, N, N, N,  N, N, N, N,  N, N, N, N,
	N, N, N, N,  N, N, N, N,  N, N, N, N,  N, N, N, N,
      };
      */


      void setReference(std::string &ref) {
	if(ref.size() == 0) 
	  return;
	
	reference.resize(0);
	for(char nt : ref) {
	  reference.push_back(char2base(nt));
	}
      }



      void merge_families(family_id source, family_id dest) {
	if(fam_rels.find(source) == fam_rels.end() || fam_rels.find(dest) == fam_rels.end())
	  return;
	else {
	  std::vector<Member*> old_fam = fam_rels.find(source)->second;
	  std::vector<Member*> new_fam = fam_rels.find(dest)->second;
	  for(Member *m : old_fam) {
	    new_fam.push_back(m);
	    m->fid = dest;
	  }
	  fam_rels.erase(source);
	}
      }
      
      

      size_t get_gender(Gender sex) {
	size_t ret;
	switch(sex){
	case Gender::Male: ret = 1; break;
	case Gender::Female: ret = 2; break;
	default: ret = 0; break;
	}
	return ret;
      }

      unsigned int id_list = 1;
      unsigned int fam_list = 1;
      std::unordered_map<family_id, std::vector<Member*>> fam_rels;
      std::vector<Member*> members;

      DataFormat dataformat_;
      // TODO: Either keep reference in file stream, or compact into 2 bits.
      //std::string reference;
      std::vector<Base> reference;

      double pop_prior_mutation = .1;
      double theta_ = 0.1;
      std::array<double, 4> nuc_freqs_ = {0.3, 0.2, 0.2, 0.3};
      double ref_weight_ = 1.0;
      double transitons_mut_ = 0.015;
      double transversion_mut_ = 0.005;

      double somatic_mutation_rate_ = 0.01;
      double homozygote_match_ = 0.98;
      double heterozygote_match_ = 0.99;
      std::vector<std::piecewise_constant_distribution<>> pop_priors_dists;
      std::mt19937 ran_generator;
      
      std::string chrom = "20";
      size_t depth = 15;
    };

  
  } // namespace sim
} // namespace dng


#endif
