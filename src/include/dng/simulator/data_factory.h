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


class Factory1 {
public:

	/*
	static SimBuilder *newInstance(DataFormat df = DataFormat::None) {
		return new Factory(df);
	}
	*/

	/*
	Member *AddTrio(Member *child, Member *mom, Member *dad) {
		return AddTrio(child, (child == nullptr ? Gender::Unknown : child->sex), mom, dad);
	}
	*/

	/*
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
	*/

	/*
	void AddLibrary(Member *m, std::string &libname, size_t depth) {


	}
	*/

	/*
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
				   << get_gender(mem->sex) << "\t"
				   << (mem->get_family_name() + mem->name)  << std::endl;
		}
	}
	*/

	/*
    void publishData() {
		//createPriorsDist();
		//initGameteDNA();
		//createLibraryMutations();
		//publishVCFData();

		//////////////////////////
		std::cout << "ref:   ";
		for(Base b : reference) {
			std::cout << nt2char[b];
		}
		std::cout << std::endl;
		for(Member *m : members) {
			std::cout << m->name << " ";
			for(size_t chrom : {0, 1}) {
				for(int site = 0; site < reference.size(); site++) {
					Base base = m->get_gamete_nt(chrom, site);
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
    }
    */


      
      
	genotype_dist genotype_weights(double theta, std::array<double, 4> nuc_freq, std::array<double,4> prior) {
		std::array<double, 16> weights;
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

	/*
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
	*/
      
	/*
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
	 */
      
	std::array<std::pair<Base, Base>, 16> genotypes = {{{A, A}, {A, C}, {A, G}, {A, T},
												        {C, A}, {C, C}, {C, G}, {C, T},
												        {G, A}, {G, C}, {G, G}, {G, T},
												        {T, A}, {T, C}, {T, G}, {T, T}}};



	/*
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
	*/


	Base transitions[4] = {G, A, T, C};
	Base transversions[4][2] = {{C, T}, {C, T}, {A, G}, {A, G}};
      

	char mutation_options[4][3] = {{C, G, T}, {A, G, T}, {A, C, T}, {A, C, G}};

	/*
	void createLibraryMutations() {
		double interval[] = {0, 1, 2, 3};
		double probs[] =  {transitons_mut_, transversion_mut_, (1.0 - transitons_mut_ - transversion_mut_)};
		std::piecewise_constant_distribution<> dist(std::begin(interval), std::end(interval), std::begin(probs));
		std::random_device rd;
		std::mt19937 gen(rd());

		for(Member *m : members) {
			// If user hasn't specified any libraries create a default one
			if(m->libraries.size() == 0) {
				std::string libname = std::string(m->name) + "-" + m->get_family_name();
				m->add_library(libname, 10);
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
							Base oldbase = m->get_lib_dna(l_indx, chrom, site);
							if(oldbase == REF)
							oldbase = reference[site];
							m->update_lib_dna(l_indx, chrom, site, transitions[oldbase]);
						}
						else {
							// Transversion
							Base oldbase = m->get_lib_dna(l_indx, chrom, site);
							if(oldbase == REF)
							oldbase = reference[site];

							int r = rand() % 2; // 2 possible choices for each transverison
							m->update_lib_dna(l_indx, chrom, site, transversions[oldbase][r]);
						}
					}
				}
			}
		}
	}
	*/



	/*
	void createFounderDNA(Member *mem) {
		// Go through each site in the reference contig and randomly select a genotype based on the reference.
		for(size_t site = 0; site < reference.size(); site++) {
			Base ref = reference[site];
			std::pair<Base, Base> gt = genotypes[floor(pop_priors_dists[ref](ran_generator))];

			// Only need to keep track if site is different from the reference.
			if(gt.first != ref) {
				mem->update_gamete_dna(0, site, gt.first);
			}
			if(gt.second != ref) {
				mem->update_gamete_dna(1, site, gt.second);
			}
		}
	}
	*/

	/*
	void createChildDNA(Member *child, Member *parent) {
		// Temporarly create a parent
		Member *tmp_parent = new Member(-1, -1, Gender::Unknown);
		createFounderDNA(tmp_parent);

		// update child
		createChildDNA(child, tmp_parent, parent);

		// delete parent
		delete tmp_parent;
	}
	*/

	/*
	void createChildDNA(Member *child, Member *mom, Member *dad) {
		if(child == NULL || mom == NULL || dad == NULL)
			return; // TODO: throw error

		// Weighted probs of type of mutation; transtion, transversion, or none.
		double interval[] = {0, 1, 2, 3};
		double probs[] =  {transitons_mut_, transversion_mut_, (1.0 - transitons_mut_ - transversion_mut_)};
		std::piecewise_constant_distribution<> dist(std::begin(interval), std::end(interval), std::begin(probs));
		std::random_device rd;
		std::mt19937 gen(rd());

		child->inherit_dna(mom, dad);
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
					Base oldbase = child->get_gamete_nt(chrom, site);
					if(oldbase == REF)
						oldbase = reference[site];
					child->update_gamete_dna(chrom, site, transitions[oldbase]);
				}
				else {
					// Transversion
					Base oldbase = child->get_gamete_nt(chrom, site);
					if(oldbase == REF)
						oldbase = reference[site];

					int r = rand() % 2; // 2 possible choices for each transverison
					child->update_gamete_dna(chrom, site, transversions[oldbase][r]);

				}
			}
		}
	}
	*/

	/**
	* Iterate through all the members of the pedigree creating DNA based on the reference and parential DNA.
	*/
	/*
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
	*/

	/*
	std::array<size_t, 4> depth_count(Member *m, size_t pos) {
		char ref = reference[pos];
		char allele1 = m->get_gamete_nt(0, pos);
		if(allele1 == ' ')
			allele1 = ref;

		char allele2 = m->get_gamete_nt(1, pos);
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
	
		return ret;
	}
	*/


	/*
	void publishVCFData(std::string &file) {
		hts::bcf::File out(file.c_str() , "w");
		out.AddHeaderMetadata("##fileformat=VCFv4.2");
		out.AddHeaderMetadata("##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">");
		out.AddContig("1", 243199373);
		for(int a = 0; a < members.size(); a++) {
			Member *mem = members[a];
			std::string sample = (boost::format("%s:LB-%s") % mem->name % mem->get_family_name()).str();
			out.AddSample(sample.c_str());
		}
		out.WriteHeader();

		for(size_t pos  = 0; pos < reference.size(); pos++) {
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

			std::vector<int> allele_order_map;//    	  samFile *fp = sam_open(file.c_str(), "w");

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
	*/

	std::string getContig(Member *m) {
		std::string contig;
		// randomly choose a contig to sample
		int chr = (rand() % 2);

		for(int pos = 0; pos < reference.size(); pos++) {
			Base ref = reference[pos];
			Base allele = m->get_gamete_nt(chr, pos);
			if(allele == REF)
				allele = ref;

			// need to add random mutation

			contig += nt2char[allele];
		}
		return contig;
	}


	void publishSAMData(std::string &file) {
		// Build the bam header file from scratch
		std::stringstream hdr_txt;
		hdr_txt << "@HD\tVN:0.1\tSO:unknown\tGO:none" << std::endl;
		hdr_txt << "@SQ\tSN:1\tLN:249250621" << std::endl;
		int id = 0;
		for(int a = 0; a < members.size(); a++) {
			Member *mem = members[a];
			for(int lib_idx = 0; lib_idx < mem->libraries.size(); lib_idx++) {
				std::string id = std::to_string(mem->mid) + "-" + mem->libraries[lib_idx].name;
				hdr_txt << "@RG" << "\t"
				<< "ID:" << id << "\t"
				<< "LB:" << mem->libraries[lib_idx].name << "\t"
				<< "SM:" << mem->name << std::endl;
			}
		}

		dng::sim::BAMFile out(file.c_str(), "w");
		out.set_header(hdr_txt.str().c_str(), hdr_txt.str().size());

		for(Member *member : members) {
			for(int l = 0; l < member->libraries.size(); l++) {
				for(int d = 0; d < member->libraries[l].depth; d++) {
					dng::sim::BAMRec rec = out.init_rec();
					rec.set_qname(chrom);

					//std::string cigar = std::to_string(reference.size()) + "M";
					rec.add_cigar(BAM_CMATCH, reference.size());

					std::string seq;
					int cn = rand() % 2;
					for(int pos = 0; pos < reference.size(); pos++) {
						Base b = member->get_lib_dna(l, cn, pos);
						if(b == REF)
							b = reference[pos];

						seq += nt2char[b];
					}
					std::cout << seq << std::endl;

					rec.set_seq(seq);

					std::string qual = "*";
					rec.set_qual(qual);

					std::string aux = /*std::string("SMZ") + */ std::to_string(member->mid) + "-" + member->libraries[l].name;
					char sm[2] = {'S', 'M'};
					//std::string aux = "NA001";
					rec.add_aux(sm, 'Z', aux);
					//rec.add_aux(aux);

					out.write_record(rec);
				}
			}
		}

		out.save();
	}


private:
	/*
	Factory(DataFormat df) : dataformat_(df){
		std::string reference_str = "TTAATAGGGCGTTGCTGGCGGGCGTTGGGTGTGGCCCGCAGTCCTGGTTGAGGATTGCCCAA";
		setReference(reference_str);
	};
	*/

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

	//DataFormat dataformat_;
	// TODO: Either keep reference in file stream, or compact into 2 bits.
	//std::string reference;
	std::vector<Base> reference;

	double pop_prior_mutation = .1;
	double theta_ = 0.1;
	std::array<double, 4> nuc_freqs_ = {{0.3, 0.2, 0.2, 0.3}};
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
