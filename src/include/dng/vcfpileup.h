/*
 * Copyright (c) 2014 Reed A. Cartwright
 * Authors:  Reed A. Cartwright <reed@cartwrig.ht>
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

#pragma once
#ifndef DNG_VCFPILEUP_H
#define DNG_VCFPILEUP_H

#include <vector>
#include <limits>
#include <unordered_map>

#include <dng/hts/bam.h>
#include <htslib/faidx.h>

#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>
#include <boost/range/size.hpp>
#include <boost/range/metafunctions.hpp>
#include <boost/range/algorithm/lower_bound.hpp>
#include <boost/range/algorithm/sort.hpp>

#include <dng/pool.h>
#include <dng/fileio.h>
#include <dng/cigar.h>
#include <dng/read_group.h>

#include <dng/hts/bcf.h>

namespace dng {
namespace pileup {
namespace vcf {


class VCFPileup {
public:
    //typedef void (callback_type)(bcf_hdr_t *, bcf1_t *);
    typedef void (callback_type)(hts::bcf::Variant&);

    VCFPileup() {
        // Not yet sure if ReadGroups should be passed in
    }

    template<typename Func>
      void operator()(hts::bcf::File &vcf, Func func);

private:

};

template<typename Func>
//void VCFPileup::operator()(const char *fname, Func func) {
void VCFPileup::operator()(hts::bcf::File &vcf, Func func) {
    // TODO? If using multiple input vcf files then we may require scanners to search each file for the same position
    // Note: synced_bcf_reader.h does exactly that

    // type erase callback function
    std::function<callback_type> call_back(func);

    //auto rec = vcf.InitVariant();
    vcf.InitReader();
    //hts::bcf::Variant rec = vcf.GetRecords();
    while(vcf.Next() > 0) {
      std::cout << "HERE" << std::endl;
      hts::bcf::Variant rec = vcf.GetRecord();
      call_back(rec);
      
    }


    /*
    // Open the VCF/BCF file
    htsFile *fp = hts_open(fname, "r");
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    bcf1_t *rec = bcf_init1();

    while(bcf_read1(fp, hdr, rec) >= 0) {
        // check that the current record is for an SNP and not an Indel, MNP, or something else
        if(bcf_get_variant_types(rec) != VCF_SNP) {
            continue;
        }
        // TODO? Check the QUAL field or PL,PP genotype fields
        // execute func
        call_back(hdr, rec);
    }
    */
}

}
}
} //namespace dng::pileup::vcf

#endif //DNG_VCFPILEUP_H

