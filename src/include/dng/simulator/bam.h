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

#ifndef CXX_DNG_SIM_BAM_H_
#define CXX_DNG_SIM_BAM_H_

#include "htslib/sam.h"
#include <string>
#include <vector>
#include <iostream>


/**
 * TODO: I created a separate BAM class rather than modify dng/hts/bam.h because the later
 * would require altering some constructors and the way it's called in call.cc. For the future
 * it could be beneficial to move these methods into dng/hts/bam.h.
 */
namespace dng {
namespace sim {

class BAMRec : public bam1_t
{
public:

	/**
	 * Set the Query template Name (first col) in the sam file.
	 * qname should not have a '\n' or '\t' or the output file will not be properly formated
	 */
	void set_qname(std::string &qname) {
		std::copy(qname.begin(), qname.end(), std::back_inserter(bam_qname));

		// htslib requires a terminal stop to the QNAME
		if(qname[qname.size()-1] != '\0')
			bam_qname.push_back('\0');
	}

	/**
	 * Add an element of the cigar string [len][op]. Needs to be added in the order they will appear
	 * in the final cigar string.
	 *
	 * Use the CIGAR related macros found in htslib/sam.h (BAM_CMATCH, BAM_CINS, etc.) as the op.
	 * TODO: look for an htslib method that converts char into the macro.
	 */
	void add_cigar(int op, unsigned int len){

		// The CIGAR in htslib is stored as a 32-bit bit-set. To make the end results more managable
		// we'll split them up into their byte components.
		uint32_t cigar = bam_cigar_gen(len, op);//
		uint8_t byte1 = (cigar & 0xFF000000) >> 24;
		uint8_t byte2 = (cigar & 0xFF0000) >> 16;
		uint8_t byte3 = (cigar & 0xFF00) >> 8;
		uint8_t byte4 = (cigar & 0xFF);

		// Not sure why they need to be stored backwards, but htslib has problems if they aren't
		bam_cigars.push_back(byte4);
		bam_cigars.push_back(byte3);
		bam_cigars.push_back(byte2);
		bam_cigars.push_back(byte1);
	}

	/**
	 * Set the SEQ field in the sam file. This method will convert 'A','C','G','T' and 'N's into their
	 * proper htslib code.
	 */
	void set_seq(std::string &seq) {

		// Convert each character in seq to it's 4-bit htslib encoding, and save 2 nucleotide codes
		// into one byte.
		int seq_size = seq.size();
		int i = 0;
		for( ; (i+1) < seq_size; i += 2) {
			uint8_t shared_base = (seq_nt16_table[seq[i]] << 4) | seq_nt16_table[seq[i+1]];
			bam_seq.push_back(shared_base);
		}

		// If the seqence has an odd number of nucleotides save the last one along with a "="
		// TODO: Unable to get htslib to properly handle a half byte. Investigate further.
		if(i != seq_size) {
		  int8_t shared_base = seq_nt16_table[seq[i]] << 4;
		  bam_seq.push_back(shared_base);
		}

	}

	/**
	 * Set the quality, phred-scale, string. Needs to have the same length as the seq string.
	 * TODO: Figure out how to get htslib to put just an '*' when there is no qual (as per v1 specs).
	*/
	void set_qual(std::string &qual) {
		for(int i = 0; i < (bam_seq.size()*2); i++) {
			bam_qual.push_back(0);
		}
	}

	/**
	 * Add an auxilarly column tag:type:data.
	 */
	void add_aux(const char tag[2], char type, std::string &data) {
		size_t l_aux = 3 + data.size();
		bam_aux.push_back(tag[0]);
		bam_aux.push_back(tag[1]);
		bam_aux.push_back(type);
		for(char c : data) {
			bam_aux.push_back(c);
		}
		if(data[data.size()-1] != '\0')
			bam_aux.push_back('\0');
	}
public:
	// TODO: It's probably faster to store the data fields directly as uint8_t[] and just memcpy
	// directly into the data array.
	std::vector<uint8_t> bam_qname;
	std::vector<uint8_t> bam_cigars;
	std::vector<uint8_t> bam_seq;
	std::vector<uint8_t> bam_qual;
	std::vector<uint8_t> bam_aux;
};

class BAMFile
{
public:
	BAMFile(const char *file, const char *mode) {
		fp_ = sam_open(file, mode);
	}

	/**
	 * Need to build the header by hand as a string, then parse it into a bam_hdr_t object. Currently
	 * htslib doesn't appear to have functions for setting the header components individually.
	 *
	 */
	void set_header(const char *text, int l_text) {
		hdr_ = sam_hdr_parse(l_text, text); // this only copies the @SQ fields into the header

		// Copy text into bam_hdr_t
		hdr_->text = (char *)malloc(sizeof(char)*(l_text));
		strcpy(hdr_->text, text);
		sam_hdr_write(fp_, hdr_);
	}

	/**
	 * Returns a bam record/line for processings and writing.
	 */
	BAMRec init_rec() {
		return BAMRec();
	}

	void write_record(BAMRec &rec) {

		// Stores information about where in the data array each field is located
		rec.core.l_qname = rec.bam_qname.size();
		rec.core.l_qseq = rec.bam_seq.size()*2; // A nucleotide is 4 bits
		rec.core.n_cigar = rec.bam_cigars.size()/4; // A cigar string is 32 bits
		size_t data_size = rec.bam_qname.size() + rec.bam_cigars.size() + rec.bam_seq.size() + rec.bam_qual.size() + rec.bam_aux.size();
		rec.l_data = data_size; // actual size
		rec.m_data = data_size; // maximum size

		// Copy all of the 8 bit fields into the data array
		uint8_t *data = new uint8_t[data_size];
		uint8_t *data_ptr = data;
		for(uint8_t d : rec.bam_qname) {
			*(data_ptr++) = d;
		}

		for(uint32_t d : rec.bam_cigars) {
			*(data_ptr++) = d;
		}

		for(uint8_t d : rec.bam_seq) {
			*(data_ptr++) = d;
		}

		for(uint8_t d : rec.bam_qual) {
			*(data_ptr++) = d;
		}

		for(uint8_t d : rec.bam_aux) {
			*(data_ptr++) = d;
		}

		rec.data = data;

		sam_write1(handle(), hdr(), &rec);
	}

	bam_hdr_t *hdr() {
		return hdr_;
	}

	samFile *handle() {
		return fp_;
	}

	void save() {
	  sam_close(fp_);
	}

private:
	bam_hdr_t *hdr_;
    samFile *fp_;
};

} // namespace sim
} // namespace dng

#endif
