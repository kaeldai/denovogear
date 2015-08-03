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

#include <dng/simulator/ped_member.h>
#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>

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


class SimBuilder {

public:
	Member *AddTrio(Member *child, Member *mom, Member *dad) {
		return AddTrio(child, (child == nullptr ? Gender::Unknown : child->sex), mom, dad);
	}


	Member *AddTrio(Member *child, Gender sex, Member *mom, Member *dad) {
		// Make sure the trio includes a child
		if(!IS_EMPTY(child) && IS_UNKNOWN(child)) {
			// TODO: Attempting to add a childless mom/dad to pedigree. Should this be valid? Currently just skipping
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
			family = createFamilyID();
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
			mom = new Member(createMemberID(), family, Gender::Female);
			members.push_back(mom);
		}
		if(IS_EMPTY(dad)){
			dad = new Member(createMemberID(), family, Gender::Male);
			members.push_back(dad);
		}
		if(IS_EMPTY(child)) {
			child = new Member(createMemberID(), family, sex);
			members.push_back(child);
		}

		child->mom = (IS_UNKNOWN(mom) ? nullptr : mom);
		child->dad = (IS_UNKNOWN(dad) ? nullptr : dad);
		return child;
	}

	void AddLibrary(Member *m, std::string &libname, size_t depth) {


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



	template<typename Stream>
	void publishPed(Stream &output) {
		for(int a = 0; a < members.size(); a++) {
			Member *mem = members[a];

			// Get the names of the mom and dad,
			std::string mom = (mem->mom == nullptr ? PED_NOPARENT : mem->mom->name);
			std::string dad = (mem->dad == nullptr ? PED_NOPARENT : mem->dad->name);

			output << mem->get_family_name() << "\t"
				   << mem->name << "\t"
				   << dad << "\t"
				   << mom << "\t"
				   << getGender(mem->sex) << "\t"
				   << (mem->get_family_name() + mem->name)  << std::endl;
		}
	}

	void publishData(std::string &fname, SeqFormat format, DataScheme scheme = SINGLE_FILE) {
		createSeqData();

		if(scheme == SINGLE_FILE) {
			// Create a single file
			switch(format) {
			case BAM: publishDataSAM(fname.c_str(), "wb", members); break;
			case SAM: publishDataSAM(fname.c_str(), "w", members); break;
			case BCF: publishDataVCF(fname.c_str(), "wb", members); break;
			case VCF: publishDataVCF(fname.c_str(), "w", members); break;
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
			for(Member *m : members) {
				std::string file = base + "_" + std::to_string(m->mid) + postfix;
				std::vector<Member*> mems = {m};
				switch(format) {
				case BAM: publishDataSAM(fname.c_str(), "wb", members); break;
				case SAM: publishDataSAM(fname.c_str(), "w", members); break;
				case BCF: publishDataVCF(fname.c_str(), "wb", members); break;
				case VCF: publishDataVCF(fname.c_str(), "w", members); break;
				}
			}
		}
	}

	virtual void setReference(std::string &seqence, std::string &chrom, size_t start_pos) = 0;


	virtual void setReference(std::string &fasta, std::string &range) = 0;


	virtual void createSeqData() = 0;


	virtual void publishDataSAM(const char *file, const char *mode, std::vector<Member*> &mems) = 0;


	virtual void publishDataVCF(const char *file, const char *mode, std::vector<Member*> &mems) = 0;


	virtual ~SimBuilder() {

	}

private:

	// Returns a unique id that can be assigned to member of pedigree
	member_id createMemberID() {
		return mid_list++;
	}

	// Returns a unique id that can be assigned to family of pedigree
	family_id createFamilyID() {
		family_id id = fid_list++;

		// Add family ID to index, even though at the moment there are no members of said family
		family_index[id] = std::vector<Member*>();
		return id;
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


	// Put all members of pedigree family 'source' into pedigree family 'dest'.
	void mergeFamilies(family_id source, family_id dest) {
		// make sure source and destination IDs exist
		if(family_index.find(source) == family_index.end() || family_index.find(dest) == family_index.end())
			return;
		else {

			std::vector<Member*> old_fam = family_index.find(source)->second;
			std::vector<Member*> new_fam = family_index.find(dest)->second;
			for(Member *m : old_fam) {
				new_fam.push_back(m);
				m->fid = dest;
			}
			family_index.erase(source);
		}
	}

protected:
	// Just use an incrementing list of integers for assigning unique member and family ids. Pedigree building is not multithreaded.
	size_t mid_list = 1;
	size_t fid_list = 1;

	std::vector<Member*> members; // List of all members in the pedigree
	std::unordered_map<family_id, std::vector<Member*>> family_index; // Keep track of which families contain which members

	std::vector<Base> reference;
	std::string chrom_;
	size_t start_pos_;
	int sam_seq_len_;
};


} // namespace sim
} // namespace dng




#endif /* SIMULATOR_H_ */
