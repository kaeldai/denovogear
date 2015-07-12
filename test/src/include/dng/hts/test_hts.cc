#define BOOST_TEST_MODULE "hts"

#include <boost/test/unit_test.hpp>
#include "dng/hts/hts.h"
#include <iostream>
#include <htslib/vcf.h>
#include <stdio.h>


std::string tmpFileName() {
  char fname[] = "unit_test_XXXXXX";
  int fd = mkstemp(fname);
  assert(fd != -1);
  return std::string(fname);
}


std::string buildTmpHtsFile(bool bin) {
  std::string fname = tmpFileName();
  std::string mode = std::string("w") + (bin ? "b" : "");

  htsFile *fp = hts_open(fname.c_str(), mode.c_str());
  bcf_hdr_t *hdr = bcf_hdr_init("w");
  bcf_hdr_append(hdr, "##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species=\"Homo sapiens\",taxonomy=x>");
  bcf_hdr_append(hdr, "##phasing=partial");
  bcf_hdr_append(hdr, "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">");
  bcf_hdr_append(hdr, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">");
  bcf_hdr_append(hdr, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">");
  bcf_hdr_append(hdr, "##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">");
  bcf_hdr_append(hdr, "##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership, build 129\">");
  bcf_hdr_append(hdr, "##INFO=<ID=H2,Number=0,Type=Flag,Description=\"HapMap2 membership\">");
  bcf_hdr_append(hdr, "##FILTER=<ID=q10,Description=\"Quality below 10\">");
  bcf_hdr_append(hdr, "##FILTER=<ID=s50,Description=\"Less than 50% of samples have data\">");
  bcf_hdr_append(hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
  bcf_hdr_append(hdr, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">");
  bcf_hdr_append(hdr, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">");
  bcf_hdr_append(hdr, "##FORMAT=<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">");
  bcf_hdr_append(hdr, "##FORMAT=<ID=TS,Number=1,Type=String,Description=\"Test String\">");
  bcf_hdr_add_sample(hdr, "NA00001");
  bcf_hdr_add_sample(hdr, "NA00002");
  bcf_hdr_add_sample(hdr, "NA00003");
  bcf_hdr_add_sample(hdr, NULL);
  bcf_hdr_write(fp, hdr);

  bcf1_t *rec = bcf_init1();
  rec->rid = bcf_hdr_name2id(hdr, "20");
  rec->pos = 14369;
  bcf_update_id(hdr, rec, "rs6054257");
  bcf_update_alleles_str(hdr, rec, "G,A");
  rec->qual = 29;
  int32_t tmpi = bcf_hdr_id2int(hdr, BCF_DT_ID, "PASS");
  bcf_update_filter(hdr, rec, &tmpi, 1);
  tmpi = 3;
  bcf_update_info_int32(hdr, rec, "NS", &tmpi, 1);
  tmpi = 14;
  bcf_update_info_int32(hdr, rec, "DP", &tmpi, 1);
  float tmpf = 0.5;
  bcf_update_info_float(hdr, rec, "AF", &tmpf, 1);
  bcf_update_info_flag(hdr, rec, "DB", NULL, 1);
  bcf_update_info_flag(hdr, rec, "H2", NULL, 1);
  int32_t *tmpia = (int*)malloc(bcf_hdr_nsamples(hdr)*2*sizeof(int));
  tmpia[0] = bcf_gt_phased(0);
  tmpia[1] = bcf_gt_phased(0);
  tmpia[2] = bcf_gt_phased(1);
  tmpia[3] = bcf_gt_phased(0);
  tmpia[4] = bcf_gt_unphased(1);
  tmpia[5] = bcf_gt_unphased(1);
  bcf_update_genotypes(hdr, rec, tmpia, bcf_hdr_nsamples(hdr)*2);
  tmpia[0] = 48;
  tmpia[1] = 48;
  tmpia[2] = 43;
  bcf_update_format_int32(hdr, rec, "GQ", tmpia, bcf_hdr_nsamples(hdr));
  tmpia[0] = 1;
  tmpia[1] = 8;
  tmpia[2] = 5;
  bcf_update_format_int32(hdr, rec, "DP", tmpia, bcf_hdr_nsamples(hdr));
  tmpia[0] = 51;
  tmpia[1] = 51;
  tmpia[2] = 51;
  tmpia[3] = 51;
  tmpia[4] = bcf_int32_missing;
  tmpia[5] = bcf_int32_missing;
  bcf_update_format_int32(hdr, rec, "HQ", tmpia, bcf_hdr_nsamples(hdr)*2);
  char *tmp_str[] = {"String1","SomeOtherString2","YetAnotherString3"};
  bcf_update_format_string(hdr, rec, "TS", (const char**)tmp_str, 3);
  bcf_write1(fp, hdr, rec);
  
  //std::cout << "fp->is_bin = " << fp->is_bin << std::endl;
  

  hts_close(fp);

  return fname;
}



BOOST_AUTO_TEST_CASE(Test_BCF)
{
  // Open new bcf file for writing
  std::string fname = tmpFileName();
  hts::File bcfw(fname.c_str(), "wb");
  BOOST_CHECK(bcfw.is_open() == 1);
  BOOST_CHECK(bcfw.is_bin() == 1);
  BOOST_CHECK(bcfw.is_write() == 1);
  BOOST_CHECK(bcfw.is_cram() == 0);
  BOOST_CHECK(strcmp(bcfw.name(), fname.c_str()) == 0);
  remove(fname.c_str());

  // Test Reading BCF
  fname = buildTmpHtsFile(true);
  hts::File bcfr(fname.c_str(), "rb");
  BOOST_CHECK(bcfr.is_open() == 1);
  BOOST_CHECK(bcfr.is_bin() == 1);
  BOOST_CHECK(bcfr.is_write() == 0);
  BOOST_CHECK(bcfr.is_cram() == 0);
  BOOST_CHECK(strcmp(bcfr.name(), fname.c_str()) == 0);
  remove(fname.c_str());
}



BOOST_AUTO_TEST_CASE(Test_VCF)
{
  // Open new VCF file for writing
  std::string fname = tmpFileName();
  hts::File vcfw(fname.c_str(), "w");
  BOOST_CHECK(vcfw.is_open() == 1);
  BOOST_CHECK(vcfw.is_bin() == 0);
  BOOST_CHECK(vcfw.is_write() == 1);
  BOOST_CHECK(vcfw.is_cram() == 0);
  BOOST_CHECK(strcmp(vcfw.name(), fname.c_str()) == 0);
  remove(fname.c_str());

  // Open existing VCF file for reading
  fname = buildTmpHtsFile(false);
  hts::File vcfr(fname.c_str(), "r");
  BOOST_CHECK(vcfr.is_open() == 1);
  BOOST_CHECK(vcfr.is_bin() == 0);
  BOOST_CHECK(vcfr.is_write() == 0);
  BOOST_CHECK(vcfr.is_cram() == 0);
  BOOST_CHECK(strcmp(vcfr.name(), fname.c_str()) == 0);
  remove(fname.c_str());

}

// TODO: Test CRAM format
