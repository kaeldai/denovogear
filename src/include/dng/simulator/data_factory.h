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


    // Format of output data file. 'None' indicates only PED file is generated.
    enum DataFormat { VCF, BCF, SAM, BAM, CRAM, None };

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

      void addLibrary(std::string &lname) {

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

      
      void publishData() {

	std::cout << "reference: " << reference << std::endl;
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


      void publishVCFData(std::string &file) {
	hts::bcf::File out(file.c_str(), "w");
	out.AddHeaderMetadata("##fileformat=VCFv4.2");
	out.AddContig("1", 243199373);
	for(int a = 0; a < members.size(); a++) {
	  Member *mem = members[a];
	  std::string sample = (boost::format("%s:LB-%s") % mem->name % mem->family_name()).str();
	  out.AddSample(sample.c_str());
	}
	out.WriteHeader();

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

    };

  
  } // namespace sim
} // namespace dng


#endif
