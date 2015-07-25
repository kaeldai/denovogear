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

    // Used to keep track of differences between member's DNA and the reference
    typedef std::unordered_map<size_t, char> reference_map;

    typedef std::array<double, 16> genotype_dist;

    // Format of output data file. 'None' indicates only PED file is generated.
    enum DataFormat { VCF, BCF, SAM, BAM, CRAM, None };

    enum Nucleotide : uint8_t { A = 0, C, G, T, N, REF, UNASSIGNED };

    // The number of data files generated.
    // TODO: Currenlty only one large file implemented
    enum FileScheme { 
      ONE_FILE, // One file containing all samples and libraries
      ONE_PER_SAMPLE, // Each sample gets its own library
      ONE_PER_LIB // Each sample:library pair gets a seperate file
    };



    // TODO: May want to make this a child of dng::io::Pedigree::Member?
    struct Member {
      member_id mid; // interal usage
      std::string name; // name that shows up in pedigree and data file
      family_id fid;
      Gender sex;
      Member *mom; 
      Member *dad;
      bool hasDNA; // Used when generating the pedigree contigs

      // Used to keep track of differences between member's DNA and the reference. We need one for each chromatid.
      // TODO: Need a structure that better handles recombination. 
      //       An NxM matrix, where N = site location and M = Member, may be more compact and efficient.
      typedef std::unordered_map<size_t, char> reference_map;
      reference_map c1;
      reference_map c2;
      reference_map library;
      //std::vector<reference_map> library_mutations;
      //reference_map somatic1;
      //reference_map somatic2;

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

      void updateDNA(int cc, size_t site, char newNuc) {
	if(cc == 1)
	  c1.insert({site, newNuc});
	else
	  c2.insert({site, newNuc});
      }

      void inheritDNA(Member *mom, Member *dad) {
	int cn = rand() % 2;
	reference_map chrm1 = (cn == 0 ? mom->c1 : mom->c2);
	cn = rand() % 2;
	reference_map chrm2 = (cn == 0 ? dad->c1 : dad->c2);
	cn = rand() % 2;
	if(cn == 0) {
	  c1 = reference_map(chrm1);
	  c2 = reference_map(chrm2);
	}
	else{
	  c1 = reference_map(chrm2);
	  c2 = reference_map(chrm1);
	}
      }

      char getNuc(int cc, size_t site) {
	std::unordered_map<size_t, char> *chrom;
	if(cc == 1)
	  chrom = &c1;
	else
	  chrom = &c2;

	std::unordered_map<size_t, char>::const_iterator nuc = chrom->find(site);
	if(nuc == chrom->end())
	  return ' ';
	else
	  return nuc->second;

	//if(cc == 1)
	  //c1.insert({site, newNuc});
	  //else
	  //c2.insert({site, newNuc});
	//return 'N';
      }

      void addLibrary() {

      }

      void somatic_mutation(int chrom, int site, char nuc) {
	library[site*chrom] = nuc;
	//char nuc = getNuc(chrom, site);

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

	/*
	std::cout << "  " << ref << " :";
	for(int a = 0; a < 16; a++) {
	  std::cout << weights[a] << "\t";
	}
	std::cout << std::endl;

	double interval[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
	std::piecewise_constant_distribution<> dist(std::begin(interval), std::end(interval), weights.begin());// std::begin(weights));
	std::random_device rd;
	std::mt19937 gen(rd());
	std::cout << dist(gen) << std::endl;
	*/

	//std::cout << "  " << ref << " : " << weights << std::endl;


	/*
	for(int nuc1 = 0; nuc1 < 4; nuc1++) {
	  for(int nuc2 = 0; nuc2 < 4; nuc2++) {
	    double alp_ = alpha[nuc1]*alpha[nuc2];
	    if(nuc1 == nuc2)
	      
	  }
	}
	*/

	      /*
	weights[0] = alpha[0]*(1.0 + alpha[0]) / alpha_sum2; // AA
	weights[1]  2.0 * alpha[0]*(alpha[1]) / alpha_sum2, // AC
	  2.0 * alpha[0]*(alpha[2]) / alpha_sum2, // AG
	  2.0 * alpha[0]*(alpha[3]) / alpha_sum2, // AT
	  alpha[1]*(1.0 + alpha[1]) / alpha_sum2, // CC
	  2.0 * alpha[1]*(alpha[2]) / alpha_sum2, // CG
	  2.0 * alpha[1]*(alpha[3]) / alpha_sum2, // CT
	  alpha[2]*(1.0 + alpha[2]) / alpha_sum2, // GG
	  2.0 * alpha[2]*(alpha[3]) / alpha_sum2, // GT
	  alpha[3]*(1.0 + alpha[3]) / alpha_sum2; // TT
	//return ret;
	std::cout << "  " << ref << " : " << weights.transpose() << std::endl;

	double interval[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
	double weightsd[] = weights.transpose().data();
	//std::piecewise_constant_distribution<> dist(std::begin(interval), std::end(interval), std::begin(weights));
	*/
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

      

      std::pair<char, char> index2Genotype(int index) {
	std::vector<std::pair<char, char>> genotypes = {{'A','A'}, {'A','C'}, {'A','G'}, {'A','T'},
							{'C','A'}, {'C','C'}, {'C','G'}, {'C','T'},
							{'G','A'}, {'G','C'}, {'G','G'}, {'G','T'},
							{'T','A'}, {'T','C'}, {'T','G'}, {'T','T'}};
	return genotypes[index];
      }

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

	/*
	std::piecewise_constant_distribution<> distA(std::begin(interval), std::end(interval), ref_weights[0].begin());
	std::piecewise_constant_distribution<> distC(std::begin(interval), std::end(interval), ref_weights[1].begin());
	std::piecewise_constant_distribution<> distG(std::begin(interval), std::end(interval), ref_weights[2].begin());
	std::piecewise_constant_distribution<> distT(std::begin(interval), std::end(interval), ref_weights[3].begin());
	*/

	std::cout << "           ";
	for(int a = 0; a*10 < reference.size(); a++) {
	  std::cout << a << "        ";
	}
	std::cout << std::endl;
	std::cout << "reference: " << reference << std::endl;
	std::random_device rd;
	std::mt19937 gen(rd());
	for(int a = 0; a < members.size(); a++) {
	  Member *m = members[a];
	  if(m->mom == NULL && m->dad == NULL) {
	    for(int b = 0; b < reference.size(); b++) {
	      char ref = reference[b];
	      int ref_idx = Allele2Index(ref);
	      std::pair<char, char> genotype = index2Genotype(floor(dists[ref_idx](gen)));
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
	      std::cout << m->getNuc(1, a);
	    }
	    std::cout << std::endl;
	    
	    std::cout << "           ";
	    for(size_t a = 0; a < reference.size(); a++) {
	      std::cout << m->getNuc(2, a);
	    }
	    std::cout << std::endl;
	  }

	}

	/*
	genotype_dist weights_refA = genotype_weights(theta_, nuc_freqs_, {ref_weight_, 0, 0, 0});
	genotype_dist weights_refC = genotype_weights(theta_, nuc_freqs_, {0, ref_weight_, 0, 0});
	genotype_dist weights_refG = genotype_weights(theta_, nuc_freqs_, {0, 0, ref_weight_, 0});
	genotype_dist weights_refT = genotype_weights(theta_, nuc_freqs_, {0, 0, 0, ref_weight_});
	double interval[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
	std::piecewise_constant_distribution<> dist(std::begin(interval), std::end(interval), weights.begin());// std::begin(weights));
	std::random_device rd;
	std::mt19937 gen(rd());
	*/
      }

      inline char transition(char nuc){ 
	char ret;
	switch(nuc) {
	case 'A': ret = 'G'; break;
	case 'G': ret = 'A'; break;
	case 'C': ret = 'T'; break;
	case 'T': ret = 'C'; break;
	default: ret = 'N';
	}
	return ret;
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
		//char ref = reference[b];
		//int ref_idx = Allele2Index(ref);
		//std::pair<char, char> genotype = index2Genotype(floor(dists[ref_idx](gen)));
		int mut_type = floor(dist(gen));
		//std::cout << mut_type << std::endl;
		if(mut_type == 2)
		  continue;
		else if(mut_type == 0) {
		  char oldbase = m->getNuc(1, b);
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
	      std::cout << m->name << "   : ";
	      for(size_t a = 0; a < reference.size(); a++) {
		std::cout << m->getNuc(1, a);
	      }
	      std::cout << std::endl;
	      
	      std::cout << "           ";
	      for(size_t a = 0; a < reference.size(); a++) {
		std::cout << m->getNuc(2, a);
	      }
	      std::cout << std::endl;
	      
	    }
	  }
	} while(remaining_members > 0);

      }

      void createLibraryMutations() {
	std::uniform_real_distribution<double> unif(0, 1.0);
	std::random_device rd;


	double interval[] = {0, 1, 2, 3};
	//double probs[] =  {transitons_mut_, transversion_mut_, (1.0 - transitons_mut_ - transversion_mut_)};
	double probs[] = {0.33333, 0.33333, 1.0-0.66666};
	std::piecewise_constant_distribution<> dist(std::begin(interval), std::end(interval), std::begin(probs));
	std::random_device rd1;
	std::mt19937 gen(rd1());

	char options[4][3] = {{'C','G','T'}, {'A', 'G', 'T'}, {'A', 'C', 'T'}, {'A', 'C', 'G'}};
	
	for(int a = 0; a < members.size(); a++) {
	  for(int chromosome : {1, 2}) {
	    for(int b = 0; b < reference.size(); b++) {
	      Member *m = members[a];
	      double prob = unif(rd);
	      if(prob <= somatic_mutation_rate_) {
		char nuc = m->getNuc(chromosome, b);
		if(nuc == ' ') {
		  nuc = reference[b];
		}
		int nuc_indx = 0;
		if(nuc == 'C')
		  nuc_indx = 1;
		else if(nuc == 'G')
		  nuc_indx = 2;
		else if(nuc == 'T')
		  nuc_indx = 3;

		//std::piecewise_constant_distribution<> dist(std::begin(options[nuc_indx]), std::end(options[nuc_indx]), std::begin(probs));
		int opt = floor(dist(gen));
		char newnuc = options[nuc_indx][opt];
		m->somatic_mutation(chromosome, b, newnuc);

		//std::cout << nuc << " --> " << newnuc << std::endl;
	      } 
	    }
	  }
	}
      }


      void publishData() {

	createFounderGenotypes();
	createChildrenDNA();
	createLibraryMutations();
	//publishVCFData();
	
	

	
	/*
	std::cout.precision(5);
	std::cout << std::fixed;
	std::cout << "Ref :\tAA\tAC\tAG\tAT\tCA\tCC\tCG\tCT\tGA\tGC\tGG\tGT\tTA\tTC\tTG\tTT" << std::endl;
	*/
	/*
	genotype_probs('A', theta_, nuc_freqs_, {ref_weight_, 0, 0, 0});
	genotype_probs('C', theta_, nuc_freqs_, {0, ref_weight_, 0, 0});
	genotype_probs('G', theta_, nuc_freqs_, {0, 0, ref_weight_, 0});
	genotype_probs('T', theta_, nuc_freqs_, {0, 0, 0, ref_weight_});
	genotype_probs('N', theta_, nuc_freqs_, {0, 0, 0, 0});
	*/

	/*
	std::cout << std::endl << std::endl;



	for(int a = 0; a < members.size(); a++) {
	  Member *m = members[a];
	  if(m->mom == NULL && m->dad == NULL) {
	    //std::cout << m->name << " is a founder" << std::endl;
	    double exp_muts = reference.size()*pop_prior_mutation;
	    //std::cout << "Expected number of mutations: " << exp_muts << std::endl;
	    // randomly sample the reference for sites of mutations, can change prob. of mutation and type depending on model
	    for(int b = 0; b < exp_muts; b++) {
	      for(int chromatid : {1, 2}) {
		size_t site = rand() % reference.size();
		//int newallele = rand() % 4;
		char oldNuc = reference[site];
		char newNuc = oldNuc;
		switch(rand() % 4) {
		case 0: newNuc = 'A'; break;
		case 1: newNuc = 'C'; break;
		case 2: newNuc = 'G'; break;
		case 3: newNuc = 'T'; break;
		}
		if(oldNuc != newNuc) {
		  //std::cout << "\t" << "[" << site << "] " << oldNuc << " --> " << newNuc << std::endl;
		  //chromotid.insert({site, newNuc});
		  m->updateDNA(chromatid, site, newNuc);
		}
	      }	
	    }
	    m->hasDNA = true;
	    

	    std::cout << m->name << "   : ";
	    for(size_t a = 0; a < reference.size(); a++) {
	      std::cout << m->getNuc(1, a);
	    }
	    std::cout << std::endl;
	    
	    std::cout << "           ";
	    for(size_t a = 0; a < reference.size(); a++) {
	      std::cout << m->getNuc(2, a);
	    }
	    std::cout << std::endl;
	  }
	}

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
	      
	    }
	  }
	} while(remaining_members > 0);
	*/
	


	// Get Reference
	// for each root (founders) parents, 
	//    for each site in reference
	//       randomly mutate (SNP's, MNP, and In-Dels) sites and create heterozygousity
	// for each trio belonging to one of the founders
	//     randomly assign a parent's chromosome to the child
	//     randomly create a gametic mutation
	//     Keep track of differences between reference and child mutation
	// Repeat above using the new children.
	// Each Member should have a list of sites of mutations/heterozygousity that differ from the reference
	//
	// for each member
	//   
	
      }



      std::array<size_t, 4> depth_count(Member *m, size_t pos) {
	char ref = reference[pos];
	char allele1 = m->getNuc(1, pos);
	if(allele1 == ' ')
	  allele1 = ref;

	char allele2 = m->getNuc(2, pos);
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
	    //for(int a = 0; a < allele_order_map.size(); a++) {
	      int ad_indx = m_indx*4 + a;//allele_order_map[a];
	      //allele_depths[ad_indx] = gtcounts[allele_order_map[a]];
	      //allele_depths.push_back(gtcounts[allele_order_map[a]]);
	      allele_depths.push_back(gtcounts[ad_indx]);
	    }
	  }

	  rec.samples("AD", allele_depths);
	  out.WriteRecord(rec);
	  
	  /*
	  std::cout << "site " << pos << ": ";
	  std::cout << allele_order_str << "; ";
	  std::cout << allele_depths[0];
	  for(int a = 1; a < allele_depths.size(); a++) {
	    std::cout << ", " << allele_depths[a];	    
	  }
	  std::cout << std::endl;
	  */	  
	}
	
      }

      std::string getMContig(Member *m) {
	std::string contig;
	// randomly choose a contig to sample
	int chr = (rand() % 2) + 1;
	
	for(int pos = 0; pos < reference.size(); pos++) {
	  char ref = reference[pos];
	  char allele = m->getNuc(chr, pos);
	  if(allele == ' ') 
	    allele = ref;

	  // need to add random mutation
	  
	  contig += allele;  
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
	  std::string id = std::to_string(mem->mid);
	  //char name[8];
	  //sprintf(name, "NA%05d", id);
	  hdr_txt << "@RG\t"
		  << "ID:" << id << "\t"
		  << "LB:" << mem->name + "-" + mem->family_name() << "\t" //(std::string(name)+"-"+ring(id)) << "\t"
		  << "SM:" << mem->name << std::endl;
	  
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
	    sam_file << getMContig(mem) << "\t";
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
	reference = "TTAATAGGGCGTTGCTGGCGGGCGTTGGGTGTGGCCCGCAGTCCTGGTTGAGGATTGCCC";
      };
      int get_uid() { return id_list++; };
      family_id get_famid() { 
	family_id id = fam_list++;
	fam_rels[id] = std::vector<Member*>();
	return id;
      };

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
      std::string reference;
      double pop_prior_mutation = .1;
      double theta_ = 0.1;
      std::array<double, 4> nuc_freqs_ = {0.3, 0.2, 0.2, 0.3};
      double ref_weight_ = 1.0;
      double transitons_mut_ = 0.015;
      double transversion_mut_ = 0.005;

      double somatic_mutation_rate_ = 0.01;
      double homozygote_match_ = 0.98;
      double heterozygote_match_ = 0.99;
      
      std::string chrom = "20";
      size_t depth = 15;
    };

  
  } // namespace sim
} // namespace dng


#endif
