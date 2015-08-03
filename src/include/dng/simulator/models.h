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

#ifndef DNG_SIM_MODELS_H_
#define DNG_SIM_MODELS_H_

#include "dng/simulator/simulator.h"
#include <iostream>
#include <vector>
#include <array>
#include <random>
#include <string>
#include <tuple>


namespace dng {
namespace sim {

typedef std::array<double, 16> genotype_dist;


class DNGModel : public SimBuilder {
public:

	void setReference(std::string &seqence, std::string &chrom, size_t start_pos) {

	}


	void setReference(std::string &fasta, std::string &range) {

	}

	virtual void setParameter(std::string &name, int val) {
		if(name == "gametic_mu") {
			gametic_mu_ = val;
		}
		else if(name == "somatic_mu") {
			somatic_mu_ = val;
		}
		else if(name == "read_mu") {
			read_mu_ = val;
		}
		else if(name == "ref_weight") {
			ref_weight_ = val;
		}
		else {
			SimBuilder::setParameter(name, val);
		}
	}

	virtual void setParameter(std::string &name, std::vector<double> val) {
		if(name == "pop_priors") {
			if(val.size() < 4)
				std::cerr << "Population priors requires 4 parameters. Skipping!" << std::endl;
			else {
				for(int i : {0, 1, 2, 3})
					pop_priors_[i] = val[i];
			}
		}
		else {
			SimBuilder::setParameter(name, val);
		}
	}


	void publishDataSAM(const char *file, const char *mode, std::vector<Member*> &mems) {
		createSeqData();
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

		dng::sim::BAMFile out(file, mode);
		out.set_header(hdr_txt.str().c_str(), hdr_txt.str().size());

		for(Member *member : mems) {
			for(int l = 0; l < member->libraries.size(); l++) {
				for(int d = 0; d < member->libraries[l].depth; d++) {
					dng::sim::BAMRec rec = out.init_rec();
					rec.set_qname(chrom_);

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


	void publishDataVCF(const char *file, const char *mode, std::vector<Member*> &mems) {

		hts::bcf::File out(file , mode);
		out.AddHeaderMetadata("##fileformat=VCFv4.2");
		out.AddHeaderMetadata("##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">");

		// Take a guess of the contig full length, should add as parameter?
		uint32_t l_chrom = static_cast<uint32_t>(start_pos_ + reference.size());
		out.AddContig(chrom_.c_str(), l_chrom);
		for(int a = 0; a < mems.size(); a++) {
			Member *mem = mems[a];
			std::string sample = (boost::format("%s:LB-%s") % mem->name % mem->get_family_name()).str();
			out.AddSample(sample.c_str());
		}
		out.WriteHeader();

		for(size_t pos  = start_pos_; pos < (start_pos_ + reference.size()); pos++) {
			char ref = reference[pos];
			int ref_indx = Allele2Index(ref);
			int gt_totals[4] = {0, 0, 0, 0};
			int total_reads = 0;
			std::vector<int> gtcounts(mems.size()*4);
			for(size_t m = 0; m < mems.size(); m++) {
				Member *mem = mems[m];
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
			rec.target(chrom_.c_str());
			rec.position(pos);

			std::string allele_order_str;
			allele_order_str = ref;
			for(int indx = 1; indx < allele_order_map.size(); indx++) {
				allele_order_str += std::string(",") + Index2Allele(allele_order_map[indx]);
			}

			rec.info("LL", static_cast<float>(1.0));

			rec.alleles(allele_order_str);



			std::vector<int32_t> allele_depths;
			for(int m_indx = 0; m_indx < mems.size(); m_indx++) {
				for(int a : allele_order_map) {
					int ad_indx = m_indx*4 + a;//allele_order_map[a];
					allele_depths.push_back(gtcounts[ad_indx]);
				}
			}

			rec.samples("AD", allele_depths);
			out.WriteRecord(rec);
		}
	}


	void createSeqData() {
		createPriorsDist();
		initGameteDNA();
		createLibraryMutations();

	}

	~DNGModel() {

	}

protected:

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

protected:
	std::vector<std::piecewise_constant_distribution<>> pop_priors_dists;
	double theta_ = 0.1;
	//double ref_weight_ = 1.0;
	double transitons_mut_ = 0.015;
	double transversion_mut_ = 0.005;
	double somatic_mutation_rate_ = 0.01;
	double homozygote_match_ = 0.98;
	double heterozygote_match_ = 0.99;
	std::array<double, 4> nuc_freqs_ = {{0.3, 0.2, 0.2, 0.3}};
	std::mt19937 ran_generator;

	Base transversions[4][2] = {{C, T}, {C, T}, {A, G}, {A, G}};
	Base transitions[4] = {G, A, T, C};
	int depth = 10;

	double gametic_mu_;
	double somatic_mu_;
	double read_mu_;
	double ref_weight_;
	std::array<double, 4> pop_priors_;
};


/**
 * Used for unit testing
 */
class Test1 : public SimBuilder {
public:

	void createSeqData() {

	}

	void publishDataSAM(const char *file, const char *mode, std::vector<Member*> &members) {

	}


	void publishDataVCF(const char *file, const char *mode, std::vector<Member*> &members) {

	}

	void setReference(std::string &seqence, std::string &chrom, size_t start_pos) {

	}


	void setReference(std::string &fasta, std::string &range) {

	}


	~Test1() {

	}

};


} // namespace sim
} // namespace dng


#endif
