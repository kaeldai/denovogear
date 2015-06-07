/*
 * Copyright (c) 2014-2015 Reed A. Cartwright <reed@cartwrig.ht>
 * Copyright (c) 2015 Kael Dai <kdai1@asu.edu>
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

#include <cstdlib>
#include <fstream>
#include <iterator>
#include <iosfwd>

#include <iostream>
#include <iomanip>
#include <ctime>
#include <chrono>
#include <sstream>
#include <string>

#include <boost/range/iterator_range.hpp>
#include <boost/range/algorithm/replace.hpp>

#include <boost/algorithm/string.hpp>

#include <dng/task/call.h>
#include <dng/pedigree.h>
#include <dng/fileio.h>
#include <dng/pileup.h>
#include <dng/read_group.h>
#include <dng/likelihood.h>
#include <dng/seq.h>
#include <dng/utilities.h>
#include <dng/hts/bcf.h>
#include <dng/hts/extra.h>
#include <dng/vcfpileup.h>

#include <htslib/faidx.h>
#include <htslib/khash.h>

#include "version.h"

using namespace dng::task;
using namespace dng;

// Helper function that mimics boost::istream_range
template<class Elem, class Traits> inline
boost::iterator_range<std::istreambuf_iterator<Elem, Traits> >
istreambuf_range(std::basic_istream<Elem, Traits> &in) {
    return boost::iterator_range<std::istreambuf_iterator<Elem, Traits>>(
               std::istreambuf_iterator<Elem, Traits>(in),
               std::istreambuf_iterator<Elem, Traits>());
}

// Helper function to determines if output should be bcf file, vcf file, or stdout. Also
// parses filename "bcf:<file>" --> "<file>"
std::pair<std::string, std::string> vcf_get_output_mode(
    Call::argument_type &arg) {
    using boost::algorithm::iequals;

    if(arg.output.empty() || arg.output == "-")
        return {"-", "w"};
    auto ret = hts::extra::extract_file_type(arg.output);
    if(iequals(ret.first, "bcf")) {
        return {ret.second, "wb"};
    } else if(iequals(ret.first, "vcf")) {
        return {ret.second, "w"};
    } else {
        throw std::runtime_error("Unknown file format '" + ret.second + "' for output '"
                                 + arg.output + "'.");
    }
    return {};
}

std::string vcf_timestamp() {
    using namespace std;
    using namespace std::chrono;
    std::string buffer(127, '\0');
    auto now = system_clock::now();
    auto now_t = system_clock::to_time_t(now);
    size_t sz = strftime(&buffer[0], 127, "Date=\"%FT%T%z\",Epoch=",
                         localtime(&now_t));
    buffer.resize(sz);
    auto epoch = std::chrono::duration_cast<std::chrono::milliseconds>(
                     now.time_since_epoch());
    buffer += to_string(epoch.count());
    return buffer;
}

template<typename VAL>
std::string vcf_command_line_text(const char *arg, VAL val) {
    return std::string("--") + arg + "=" + dng::util::to_pretty(val);
}

std::string vcf_command_line_text(const char *arg, std::string val) {
    return std::string("--") + arg + "=\'" + val + "\'";
}

// Helper function for writing the vcf header information
void vcf_add_header_text(hts::bcf::File &vcfout, Call::argument_type &arg) {
    using namespace std;
    string line{"##DeNovoGearCommandLine=<ID=dng-call,Version="
                PACKAGE_VERSION ","};
    line += vcf_timestamp();
    line += ",CommandLineOptions=\"";

// TODO: to_string loses precision.  Need to find/write a version
// that is as exact as possible.
#define XM(lname, sname, desc, type, def) \
	line += vcf_command_line_text(XS(lname),arg.XV(lname)) + " ";

#	include <dng/task/call.xmh>
#undef XM
    line.pop_back();
    line += "\">";
    vcfout.AddHeaderMetadata(line);

    // Add the available tags for INFO, FILTER, and FORMAT fields
    vcfout.AddHeaderMetadata("##INFO=<ID=LL,Number=1,Type=Float,Description=\"Log likelihood\">");
    vcfout.AddHeaderMetadata("##INFO=<ID=PMUT,Number=1,Type=Float,Description=\"Probability of mutation\">");
    // TODO: The commented lines are standard VCF fields that may be worth adding to dng output
    //vcfout.AddHeaderField("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">");
    //vcfout.AddHeaderField("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">");
    //vcfout.AddHeaderField("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">");
    //vcfout.AddHeaderField("##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">");
    //vcfout.AddHeaderField("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    //vcfout.AddHeaderField("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">");
    //vcfout.AddHeaderField("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">");
    //vcfout.AddHeaderField("##FORMAT=<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">");

    // AD defined http://gatkforums.broadinstitute.org/discussion/1268/how-should-i-interpret-vcf-files-produced-by-the-gatk
    vcfout.AddHeaderMetadata("##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">");
}

void vcf_add_record(hts::bcf::File &vcfout, const char *chrom, int pos,
                    const char ref,
                    double ll, double pmut, const std::vector<dng::depth5_t> &read_depths) {
    auto rec = vcfout.InitVariant();
    rec.target(chrom);
    rec.position(pos);
    //rec.filter("PASS");

    rec.info("LL", static_cast<float>(ll));
    rec.info("PMUT", static_cast<float>(pmut));

    // Based on all the samples determine what nucleotides show up and
    // the order they will appear in the REF and ALT field
    std::vector<uint16_t>
    allele_order; // List of nucleotides as they appear in the REF and ALT fields
    std::string allele_order_str; // alt_order in string format, used for SetAlleles

    std::size_t ref_index = seq::char_index(ref);
    allele_order.push_back(ref_index);
    allele_order_str = ref;
    for(std::size_t index = 1; index < 4; index++) {
        // iterate through the three remaining NTs, if the NT exists in one of the sample add it to alt_order
        int alt_allele_index = (ref_index + index) % 4;
        for(std::size_t sample = 0; sample < read_depths.size(); sample++) {
            if(read_depths[sample].counts[alt_allele_index] != 0) {
                allele_order.push_back(alt_allele_index);
                allele_order_str += std::string(",") + seq::indexed_char(alt_allele_index);
                break;
            }
        }
    }

    // Update REF, ALT fields
    rec.alleles(allele_order_str);

    // Turn allele frequencies into AD format; order will need to match REF+ALT ordering of nucleotides
    std::vector<int32_t> gtcounts;
    for(std::size_t sample = 0; sample < read_depths.size(); sample++) {
        for(std::size_t nt = 0; nt < allele_order.size(); nt++) {
            size_t allele_index = allele_order[nt];
            gtcounts.push_back(read_depths[sample].counts[allele_index]);
        }
    }
    rec.samples("AD", gtcounts);

    vcfout.WriteRecord(rec);
}

// Build a list of all of the possible contigs to add to the vcf header
std::vector<std::pair<std::string, uint32_t>> parse_contigs(
const bam_hdr_t *hdr) {
    if(hdr == nullptr)
        return {};
    std::vector<std::pair<std::string, uint32_t>> contigs;
    uint32_t n_targets = hdr->n_targets;
    for(size_t a = 0; a < n_targets; a++) {
        if(hdr->target_name[a] == nullptr) {
            continue;
        }
        contigs.emplace_back(hdr->target_name[a], hdr->target_len[a]);
    }
    return contigs;
}


// VCF header lacks a function to get sequence lengths
// So we will extract the contig lines from the input header
std::vector<std::string> extract_contigs(const bcf_hdr_t *hdr) {
    if(hdr == nullptr)
        return {};
    // Read text of header
    int len;
    std::unique_ptr<char[], void(*)(void *)> str{bcf_hdr_fmt_text(hdr, 0, &len), free};
    if(!str)
        return {};
    std::vector<std::string> contigs;

    // parse ##contig lines
    const char *text = str.get();
    if(strncmp(text, "##contig=", 9) != 0) {
        text = strstr(text, "\n##contig=");
    } else {
        text = text - 1;
    }
    const char *end;
    for(; text != nullptr; text = strstr(end, "\n##contig=")) {
        for(end = text + 10; *end != '\n' && *end != '\0'; ++end)
            /*noop*/;
        if(*end != '\n') {
            return contigs;    // bad header, return what we have.
        }
        contigs.emplace_back(text + 1, end);
    }

    return contigs;
}

// The main loop for dng-call application
// argument_type arg holds the processed command line arguments
int Call::operator()(Call::argument_type &arg) {
    using namespace std;

    // Parse pedigree from file
    dng::io::Pedigree ped;
    if(!arg.ped.empty()) {
        ifstream ped_file(arg.ped);
        if(!ped_file.is_open()) {
            throw std::runtime_error(
                "unable to open pedigree file '" + arg.ped + "'.");
        }
        ped.Parse(istreambuf_range(ped_file));
    } else {
        throw std::runtime_error("pedigree file was not specified.");
    }

    // Open Reference
    unique_ptr<char[], void(*)(void *)> ref{nullptr, free};
    int ref_sz = 0, ref_target_id = -1;
    unique_ptr<faidx_t, void(*)(faidx_t *)> fai{nullptr, fai_destroy};
    if(!arg.fasta.empty()) {
        fai.reset(fai_load(arg.fasta.c_str()));
        if(!fai)
            throw std::runtime_error("unable to open faidx-indexed reference file '"
                                     + arg.fasta + "'.");
    }

    // Parse Nucleotide Frequencies
    std::array<double, 4> freqs;
    // TODO: read directly into freqs????  This will need a wrapper that provides an "insert" function.
    // TODO: include the size into the pattern, but this makes it harder to catch the second error.
    // TODO: turn all of this into a template function that returns array<double,4>?
    {
        auto f = util::parse_double_list(arg.nuc_freqs, ',', 4);
        if(!f.second) {
            throw std::runtime_error("Unable to parse nuc-freq option. "
                                     "It must be a comma separated list of floating-point numbers.");
        }
        if(f.first.size() != 4) {
            throw std::runtime_error("Wrong number of values passed to nuc-freq. "
                                     "Expected 4; found " + std::to_string(f.first.size()) + ".");
        }
        std::copy(f.first.begin(), f.first.end(), &freqs[0]);
    }

    // quality thresholds
    int min_qual = arg.min_basequal;
    double min_prob = arg.min_prob;

    // Model genotype likelihoods as a mixture of two dirichlet multinomials
    // TODO: control these with parameters
    genotype::DirichletMultinomialMixture genotype_likelihood(
    {0.9, 0.001, 0.001, 1.05}, {0.1, 0.01, 0.01, 1.1});

    // Open input files
    dng::ReadGroups rgs;
    vector<hts::File> indata;
    vector<hts::bam::File> bamdata;
    vector<hts::bcf::File> bcfdata;
    for(auto && str : arg.input) {
        indata.emplace_back(str.c_str(), "r");
        if(indata.back().is_open()) {
            continue;
        }
        throw std::runtime_error("unable to open input file '" + str + "'.");
    }

    // Check to see if all inputs are of the same type
    const htsFormatCategory cat = indata[0].format().category;
    for(auto && f : indata) {
        if(f.format().category == cat) {
            continue;
        }
        throw std::runtime_error("mixing sequence data and variant data as input is not supported.");
    }

    // Begin writing VCF header
    auto out_file = vcf_get_output_mode(arg);
    hts::bcf::File vcfout(out_file.first.c_str(), out_file.second.c_str());
    vcf_add_header_text(vcfout, arg);

    if(cat == sequence_data) {
        // Wrap input in hts::bam::File
        for(auto && f : indata) {
            bamdata.emplace_back(std::move(f), arg.region.c_str(), arg.fasta.c_str(),
                                 arg.min_mapqual, arg.min_qlen);
        }
        // Read header from first file
        const bam_hdr_t *h = bamdata[0].header();

        // Add contigs to header
        for(auto && contig : parse_contigs(h)) {
            vcfout.AddContig(contig.first.c_str(), contig.second);
        }

        // Add each genotype/sample column
        rgs.ParseHeaderText(bamdata);
    } else if(cat == variant_data) {
        bcfdata.emplace_back(std::move(indata[0]));
        // Read header from first file
        const bcf_hdr_t *h = bcfdata[0].header();

        // Add contigs to header
        for(auto && contig : extract_contigs(h)) {
            vcfout.AddHeaderMetadata(contig.c_str());
        }
        // Add each genotype/sample column
        rgs.ParseSamples(bcfdata[0]);
    } else {
        throw runtime_error("unsupported file category.");
    }

    // Finish Header
    for(auto str : rgs.libraries()) {
        boost::replace(str, '\t', ':');
        vcfout.AddSample(str.c_str());
    }
    vcfout.WriteHeader();

    // Construct peeling algorithm from parameters and pedigree information
    dng::Pedigree peeler;
    peeler.Initialize({arg.theta, arg.mu, arg.mu_somatic, arg.mu_library,
                       arg.ref_weight, freqs
                      });
    if(!peeler.Construct(ped, rgs)) {
        throw std::runtime_error("Unable to construct peeler for pedigree; "
                                 "possible non-zero-loop pedigree.");
    }

    auto calculate = [&min_prob, &peeler, &vcfout,
                      &genotype_likelihood](const std::vector<depth5_t> &depths,
    const char *target_name, int position, char ref_base) {

        int ref_index = seq::char_index(ref_base);

        // calculate genotype likelihoods and store in the lower library vector
        double scale = 0.0, stemp;
        for(std::size_t u = 0; u < depths.size(); ++u) {
            std::tie(peeler.library_lower(u), stemp) =
                genotype_likelihood({depths[u].key}, ref_index);
            scale += stemp;
        }

        // Calculate probabilities
        double d = peeler.CalculateLogLikelihood(ref_index) + scale;
        double p = peeler.CalculateMutProbability(ref_index);

        // Skip this site if it does not meet lower probability threshold
        if(p < min_prob) {
            return;
        }

        vcf_add_record(vcfout, target_name, position, ref_base, d, p, depths);
        return;
    };

    // Pileup data
    std::vector<depth5_t> read_depths(rgs.libraries().size(), {0, 0});

    // Treat sequence_data and variant data separately
    if(cat == sequence_data) {
        const bam_hdr_t *h = bamdata[0].header();
        dng::BamPileup mpileup(rgs.groups());
        mpileup(bamdata, [&](const dng::BamPileup::data_type & data, uint64_t loc) {

            // Calculate target position and fetch sequence name
            int target_id = location_to_target(loc);
            int position = location_to_position(loc);

            if(target_id != ref_target_id && fai) {
                ref.reset(faidx_fetch_seq(fai.get(), h->target_name[target_id],
                                          0, 0x7fffffff, &ref_sz));
                ref_target_id = target_id;
            }

            // Calculate reference base
            char ref_base = (ref && 0 <= position && position < ref_sz) ?
                            ref[position] : 'N';

            // reset all depth counters
            read_depths.assign(read_depths.size(), {0, 0});

            // pileup on read counts
            // TODO: handle overflow?  Down sample?
            for(std::size_t u = 0; u < data.size(); ++u) {
                for(auto && r : data[u]) {
                    if(r.is_missing || r.qual.first[r.pos] < arg.min_basequal) {
                        continue;
                    }
                    read_depths[rgs.library_from_id(u)].counts[
                        seq::base_index(r.aln.seq_at(r.pos)) ] += 1;
                }
            }
            calculate(read_depths, h->target_name[target_id], position, ref_base);
        });
    } else if(cat == variant_data) {
        const char *fname = bcfdata[0].name();
        dng::pileup::vcf::VCFPileup vcfpileup;
        //vcfpileup(fname, [&](bcf_hdr_t *hdr, bcf1_t *rec) {
	vcfpileup(bcfdata[0], [&](hts::bcf::Variant &record) {
            // Won't be able to access ref->d unless we unpack the record first
            //bcf_unpack(rec, BCF_UN_STR);
	    record.unpack();

            // get chrom, position, ref from their fields
	    const char *chrom = record.chrom();
	    int32_t position = record.position();
	    int32_t n_alleles = record.n_alleles();
	    int32_t n_samples = record.n_samples();
	    std::cout << position << std::endl;

	    //char **alleles = record.alleles();
	    //const char ref_base = alleles[0][0];

	    //std::cout << chrom << "\t" << position << "\t" << ref_base << std::endl;

	    /*
            const char *chrom = bcf_hdr_id2name(hdr, rec->rid);
	    int32_t position = rec->pos;
            uint32_t n_alleles = rec->n_allele;
            uint32_t n_samples = bcf_hdr_nsamples(hdr);
            const char ref_base = *(rec->d.allele[0]);

            // Read all the Allele Depths for every sample into ad array
            int *ad = NULL;
            int n_ad = 0;
            int n_ad_array = 0;
            n_ad = bcf_get_format_int32(hdr, rec, "AD", &ad, &n_ad_array);

            // Create a map between the order of vcf alleles (REF+ALT) and their correct index in read_depths.counts[]
            vector<size_t> a2i;
            for(int a = 0; a < n_alleles; a++) {
                char base = *(rec->d.allele[a]);
                a2i.push_back(seq::char_index(base));
            }

            // Build the read_depths
            read_depths.assign(n_samples, {0, 0});
            for(size_t sample_ndx = 0; sample_ndx < n_samples; sample_ndx++) {
                for(size_t allele_ndx = 0; allele_ndx < n_alleles; allele_ndx++) {
                    int32_t depth = ad[n_alleles * sample_ndx + allele_ndx];
                    read_depths[sample_ndx].counts[a2i[allele_ndx]] = depth;
                }
            }
	    */

            //calculate(read_depths, chrom, position, ref_base);
        });
    } else {
        throw runtime_error("unsupported file category.");
    }

    return EXIT_SUCCESS;
}
