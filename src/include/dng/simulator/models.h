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

	DNGModel() : ref_weight_(0.9), theta_(0.05), gametic_mu_(0.0001), somatic_mu_(0.0001), read_mu_(0.0001) {
		ran_generator = std::mt19937(time(0));
		nuc_freqs_ = {{0.3, 0.2, 0.2, 0.3}};
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
		// Build the bam header file from scratch
		std::stringstream hdr_txt;
		hdr_txt << "@HD\tVN:0.1\tSO:unknown\tGO:none" << std::endl;
		hdr_txt << "@SQ\tSN:1\tLN:249250621" << std::endl;
		int id = 0;
		for(Member *member : mems) {
			for(Sample &lib : member->libraries) {
				std::string id = member->name + "-" + lib.name;
				hdr_txt << "@RG" << "\t"
				        << "ID:" << id << "\t"
				        << "LB:" << lib.name << "\t"
				        << "SM:" << member->name << std::endl;
			}
		}

		dng::sim::BAMFile out(file, mode);
		out.set_header(hdr_txt.str().c_str(), hdr_txt.str().size());


		size_t contig_len_ = 30;
		for(Member *member : mems) {
			for(int l = 0; l < member->libraries.size(); l++) {
				Sample &library = member->libraries[l];

				for(int pos = 0; pos < reference.size(); pos += contig_len_) {
					size_t l_ccontig = (pos + contig_len_ < reference.size() ? contig_len_ : reference.size() - pos);
					//std::cout << "l_ccontig = " << l_ccontig << std::endl;

					// TODO: Figure out a better data structure for storing the results. Row base look-ups are not efficient.
					std::vector<std::vector<Base>> reads;
					for(int next = 0; next < l_ccontig; next++) {
						reads.push_back(baseCall(library, next+pos));
					}

					//std::cout << "HERE" << std::endl;

					//std::cout << "library depth = " << library.depth << std::endl;
					for(int d = 0; d < library.depth; d++) {
						//std::cout << "init_rec()" << std::endl;
						dng::sim::BAMRec rec = out.init_rec();
						rec.set_qname(chrom_);

						//std::cout << "Adding cigar " << std::endl;
						rec.add_cigar(BAM_CMATCH, l_ccontig);

						std::string seq;
						for(int site = 0; site < l_ccontig; site++) {
							seq += nt2char[reads[site][d]];
						}

						rec.set_seq(seq);

						std::string qual = "*";
						rec.set_qual(qual);

						std::string aux = std::to_string(member->mid) + "-" + member->libraries[l].name;
						char sm[2] = {'S', 'M'};
						rec.add_aux(sm, 'Z', aux);

						out.write_record(rec);
					}
				}
			}
		}
		out.save();
	}


	void publishDataVCF(const char *file, const char *mode, std::vector<Member*> &mem_subset) {


		// Build the VCF/BCF header
		hts::bcf::File out(file , mode);
		out.AddHeaderMetadata("##fileformat=VCFv4.2");
		out.AddHeaderMetadata("##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">");
		uint32_t l_chrom = static_cast<uint32_t>(2*start_pos_ + reference.size()); // Take a guess of the contig full length, should add as parameter?
		out.AddContig(chrom_.c_str(), l_chrom);

		// Add a format column for each library
		size_t n_libs = 0;
		for(int a = 0; a < mem_subset.size(); a++) {
			Member *mem = mem_subset[a];
			if(mem->libraries.size() == 0) {
				// If the user hasn't created a specific library for member then add a default one
				std::string sample = (boost::format("LB%s:%s") % mem->name % mem->get_family_name()).str();
				mem->add_library(sample, default_depth_);
			}

			for(int l = 0; l < mem->libraries.size(); l++) {
				out.AddSample(mem->libraries[l].name.c_str());
				n_libs++;
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
			for(size_t m = 0; m < mem_subset.size(); m++) {
				Member *mem = mem_subset[m];
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


			/*
			std::vector<int32_t> allele_depths;
			for(int m_indx = 0; m_indx < mem_subset.size(); m_indx++) {
				Member *mem = mem_subset[m];
				for(size_t l = 0; l < mem->libraries.size(); l++) {
					for(int a : allele_order_map) {
						int ad_indx = m_indx*4 + a;//allele_order_map[a];
						allele_depths.push_back(gtcounts[ad_indx]);
					}
				}
			}
			*/

			std::vector<int32_t> allele_depths;
			for(int lindx = 0; lindx+3 < gtcounts.size(); lindx += 4) {
				for(int a : allele_order_map) {
					int gt_loc = lindx + a;
					allele_depths.push_back(gtcounts[gt_loc]);
				}
			}

			/*
			//std::cout << allele_order_str << std::endl;
			std::cout << pos << ": [";
			for(int a = 0; a < gtcounts.size(); a++) {
				std::cout << gtcounts[a];
				if(a%4 == 3) {
					std::cout << "]";
					if(a+1 < gtcounts.size())
						std::cout << "[";
				}
				else {
					std::cout << ", ";
				}
			}
			std::cout << std::endl;
			*/

			rec.samples("AD", allele_depths);
			out.WriteRecord(rec);
		}
	}

	std::vector<Base> baseCall(Sample &lib, size_t site) {
		std::array<int, 5> interval = {0, 1, 2, 3, 4};
		Base ref = reference[site];
		Base genotype[2] = {lib.get_nt(0, site), lib.get_nt(1, site)};
		if(genotype[0] == REF)
			genotype[0] = ref;
		if(genotype[1] == REF)
			genotype[1] = ref;

		//if(genotype[0] != ref || genotype[1] != ref) {
		//	std::cout << "library " << lib.name << " has variance at site " << site << std::endl;
		//}

		std::array<unsigned int, 4> prev_reads = {0, 0, 0, 0};
		std::vector<Base> reads;
		for(int d = 0; d < lib.depth; d++) {

			// TODO: If homozygote don't call rand
			int chrom = rand() % 2;
			Base allele = genotype[chrom];

			std::array<double, 5> &weights = read_err_weights[ref][allele];
			std::array<double, 4> probs;
			probs[0] = (weights[0] + prev_reads[0])/(weights[4] + d);
			probs[1] = (weights[1] + prev_reads[1])/(weights[4] + d);
			probs[2] = (weights[2] + prev_reads[2])/(weights[4] + d);
			probs[3] = (weights[3] + prev_reads[3])/(weights[4] + d);

			std::piecewise_constant_distribution<> dist(std::begin(interval), std::end(interval), probs.begin());
			int read = floor(dist(ran_generator));
			//if(read != allele) {
			//	std::cout << lib.name << " read error : " << allele << " --> " << read << std::endl;
			//}

			prev_reads[read]++;
			reads.push_back(index2Base[read]);
		}

		return reads;
	}


	/*
	void readSeq(Sample &lib) {


	}
	*/


	void createSeqData() {
		createPriorsDist();
		initGameteDNA();
		createLibraryMutations();
		createReadWeights();

		/*
		/////////////////////////////////// For Testing Only //////////////////////////
		std::cout << "       ";
		for(int a = 0; a*10 < reference.size() && a < 10; a++) {
			std::cout << a << "         ";
		}
		std::cout << std::endl;
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
		///////////////////////////////////
		*/
	}

	~DNGModel() {

	}

protected:

	void createPriorsDist() {
		genotype_dist ref_weights[] = { genotype_weights(theta_, nuc_freqs_, {ref_weight_, 0, 0, 0}),
										genotype_weights(theta_, nuc_freqs_, {0, ref_weight_, 0, 0}),
										genotype_weights(theta_, nuc_freqs_, {0, 0, ref_weight_, 0}),
										genotype_weights(theta_, nuc_freqs_, {0, 0, 0, ref_weight_}),
										genotype_weights(theta_, nuc_freqs_, {0, 0, 0, 0})};
		double interval[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};

		pop_priors_dists.emplace_back(std::begin(interval), std::end(interval), ref_weights[0].begin());
		pop_priors_dists.emplace_back(std::begin(interval), std::end(interval), ref_weights[1].begin());
		pop_priors_dists.emplace_back(std::begin(interval), std::end(interval), ref_weights[2].begin());
		pop_priors_dists.emplace_back(std::begin(interval), std::end(interval), ref_weights[3].begin());
		pop_priors_dists.emplace_back(std::begin(interval), std::end(interval), ref_weights[4].begin());
		//genotype_prior_[4] = population_prior(p.theta, p.nuc_freq, {0, 0, 0, 0});
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

	void createReadWeights() {
		double ref_bias = 1.1; // A bias the sequencer gives towards the reference allele
		double err = read_mu_/3.0; // prob read doesn't match actual allele
		double match = 1.0 - read_mu_; // prob read is correct

		for(Base ref : {A, C, G, T}) {
			for(Base allele : {A, C, G, T}) {
				std::array<double, 5> &weights = read_err_weights[ref][allele];
				weights[allele] = match;
				weights[ref] *= ref_bias;
				weights[4] = weights[0] + weights[1] + weights[2] + weights[3];
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
			return;

		child->inherit_dna(mom, dad);


		// We expect there to be "ref_len * mu" sites of mutations. At every potential site of mutation there is a 1/4th chance
		// the randomly selected new NT matches the old. By solving for N in "N*3/4 = ref_len*mu", N tells us the number of potential
		// sites of mutations we should iterate through
		size_t l_reference = reference.size();
		size_t n_mutations = ceil((l_reference*gametic_mu_+1)*1.333333); // Added +1 term to get some variance when expected mut sites is less than 3/4


		// Use the geometric distribution to determine, based on gametic mutation prob, how many sites away from some start position
		// is the next mutation. To prevent many mutations bunching up and potentially overwriting each other we update the start
		// position every time a mutation site is selected
		double ln_no_mut = log(1.0 - gametic_mu_); // log(1-mu)
		srand(time(0));
		for(int chrom : {0, 1}) {
			size_t site = rand() % l_reference; // select a random place along contig to start
			for(int m = 0; m < n_mutations; m++) {
				double pr = ((double) rand() / RAND_MAX) * gametic_mu_; // choose a random prob normalized to the geometric prob
																		// that there's a mutation 0 sites away.
				double next_site = (log(pr/gametic_mu_)/ln_no_mut); // inverse of the gemetric dist to get an non-integer k
				site = (site + size_t(next_site)) % l_reference;
				Base old_nt = child->get_gamete_nt(0, site);
				if(old_nt == REF)
					old_nt = reference[site];

				// Randomly select a base, update child genotype if new base differs from the old.
				Base new_nt = index2Base[rand() % 4];
				//std::cout << "site " << site << ": " << old_nt << " --> " << new_nt << std::endl;

				if(new_nt != old_nt) {
					// TODO: If memory becomes an issue we should erase new mutations that matches the reference
					child->update_gamete_dna(chrom, site, new_nt);
				}
			}
		}
	}

	void createLibraryMutations() {
		// Creates random mutations along the library contigs using a geometric distribution similar to
		// createChildDNA(). If not going to change then consider merging code.
		size_t l_reference = reference.size();
		size_t n_mutations = ceil((l_reference*somatic_mu_+1)*1.333333);
		double ln_no_mut = log(1.0 - somatic_mu_); // log(1-mu)
		size_t site = rand() / l_reference; // select a random place along contig to start

		// Go through each member's libraries and create unique somatic cell mutations for each library.
		for(Member *m : members) {

			// If user hasn't specified any libraries create a default one
			if(m->libraries.size() == 0) {
				std::string libname = std::string(m->name) + "-" + m->get_family_name();
				m->add_library(libname, default_depth_);
			}

			srand(time(0));
			for(size_t l_indx = 0; l_indx < m->libraries.size(); l_indx++) { // each library belong to member
				for(int chrom : {0, 1}) { // each chromosome in the library
					size_t site = rand() % l_reference; // select a random place along contig to start
					for(int mut = 0; mut < n_mutations; mut++) {
						double pr = ((double) rand() / RAND_MAX);
						double next_site = (log(pr)/ln_no_mut);
						site = (site + size_t(next_site)) % l_reference;
						Base old_nt = m->get_gamete_nt(0, site);
						if(old_nt == REF)
							old_nt = reference[site];

						// Randomly select a base, update child genotype if new base differs from the old.
						Base new_nt = index2Base[rand() % 4];
						if(new_nt != old_nt) {
							m->update_lib_dna(l_indx, chrom, site, new_nt);
							//std::cout << "Adding " << m->libraries[l_indx].name << " mutation at site " << site << std::endl;
						}
					}
				}
			}
		}
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
	double ref_weight_; // weight given to reference
	double theta_; // population diversity
	std::array<double, 4> nuc_freqs_ = {{0.3, 0.2, 0.2, 0.3}}; // nucleotide freq's in A,C,G,T order
	double gametic_mu_ = 0.0001; // prob of gamete cell line mutation
	double somatic_mu_ = 0.0001; // prob of somatic cell mutation
	double read_mu_ = 0.0001; // prob of error in sequence read;

	double default_depth_ = 10; // Should be set uniquly for each library.


	std::vector<std::piecewise_constant_distribution<>> pop_priors_dists;
	std::array<std::array<std::array<double, 5>, 4>, 4> read_err_weights;

	double transitons_mut_ = 0.015;
	double transversion_mut_ = 0.005;
	double somatic_mutation_rate_ = 0.01;
	double homozygote_match_ = 0.98;
	double heterozygote_match_ = 0.99;

	std::mt19937 ran_generator;

	//Base transversions[4][2] = {{C, T}, {C, T}, {A, G}, {A, G}};
	//Base transitions[4] = {G, A, T, C};
	//int depth = 10;




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
