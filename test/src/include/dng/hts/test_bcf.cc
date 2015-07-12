#define BOOST_TEST_MODULE "hts::bcf.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/output_test_stream.hpp>
#include "dng/hts/bcf.h"
#include "dng/hts/hts.h"
#include <htslib/vcf.h>
#include <iostream>


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








BOOST_AUTO_TEST_CASE(Test_VCF_Output)
{
  // Test all the functions, including on bad input (should return -1)
  std::string fname = tmpFileName();
  {
    hts::bcf::File vcf(fname.c_str(), "w");
    BOOST_CHECK(vcf.AddHeaderMetadata("##fileformat=VCFv4.2") == 0);
    BOOST_CHECK(vcf.AddHeaderMetadata(std::string("##FILTER=<ID=PASS,Description=\"All filters passed\">")) == 0);
    BOOST_CHECK(vcf.AddHeaderMetadata("##INFO=<ID=I1,Number=1,Type=String>") == 0);
    BOOST_CHECK(vcf.AddHeaderMetadata("##INFO=<ID=I2,Number=1,Type=Integer>") == 0);
    BOOST_CHECK(vcf.AddHeaderMetadata("##INFO=<ID=I3,Number=1,Type=Float>") == 0);
    BOOST_CHECK(vcf.AddHeaderMetadata("##INFO=<ID=I4,Number=1,Type=String>") == 0);
    BOOST_CHECK(vcf.AddHeaderMetadata("##FORMAT=<ID=S1,Number=1,Type=Float>") == 0);
    BOOST_CHECK(vcf.AddHeaderMetadata("##FORMAT=<ID=S2,Number=3,Type=Integer>") == 0);
    BOOST_CHECK(vcf.AddHeaderMetadata("##FORMAT=<ID=S3,Number=.,Type=String>") == 0);
    BOOST_CHECK(vcf.AddHeaderMetadata("##FORMAT=<ID=S4,Number=.,Type=String>") == 0);
    BOOST_CHECK(vcf.AddHeaderMetadata("string", "str") == 0);
    BOOST_CHECK(vcf.AddHeaderMetadata("int", 1) == 0);
    BOOST_CHECK(vcf.AddHeaderMetadata("float", 1.1) == 0);
    BOOST_CHECK(vcf.AddHeaderMetadata("cstring", std::string("cpp_string")) == 0);
    BOOST_CHECK(vcf.AddHeaderMetadata(NULL) == -1); // BAD INPUT
    BOOST_CHECK(vcf.AddHeaderMetadata(NULL, NULL) == -1); // BAD INPUT
    BOOST_CHECK(vcf.AddContig("1", 249250621) == 0);
    BOOST_CHECK(vcf.AddSample("sample1") == 0);
    BOOST_CHECK(vcf.AddSample("sample2") == 0);
    BOOST_CHECK(vcf.AddSample("sample3") == 0);
    BOOST_CHECK(vcf.WriteHeader() == 0);
    
    hts::bcf::Variant rec = vcf.InitVariant();
    rec.target("1");
    rec.position(1);
    rec.id("rs999");
    BOOST_CHECK(rec.alleles("G,T,C") == 0);
    rec.quality(0.0);
    rec.filter("PASS");
    BOOST_CHECK(rec.info("I1", "1") == 0);
    BOOST_CHECK(rec.info("I2", 2) == 0);
    BOOST_CHECK(rec.info("I3", 3.0f) == 0);
    BOOST_CHECK(rec.info("I4", std::string("str")) == 0);
    BOOST_CHECK(rec.samples("S1", std::vector<float>{1.2, 2.3, 4.5}) == 0);
    BOOST_CHECK(rec.samples("S2", std::vector<int32_t>{1, 2, 3, 4, 5, 6, 7, 8, 9}) == 0);
    BOOST_CHECK(rec.samples("S3", std::vector<const char*>{"one,two", "three,four", "five,six"}) == 0);
    BOOST_CHECK(rec.samples("S4", std::vector<std::string>{std::string("eins"), std::string("zwei"), std::string("drei")}) == 0);
    vcf.WriteRecord(rec);
    // TODO: Need to add a close() to bcf::File
  }


  // Read the file back and verify all the values above
  // TODO: Find a different way to verify the output. Using both htslib for reading the file and verifying it is problematic.
  //       Tried using output_test_stream, but is not portable and htslib does not guarentee to match order.
  htsFile *fp = hts_open(fname.c_str(), "r");
  bcf_hdr_t *hdr = bcf_hdr_read(fp);
  bcf1_t *rec = bcf_init1();
  BOOST_CHECK(bcf_read1(fp, hdr, rec) >= 0);
  BOOST_CHECK(strcmp(bcf_hdr_id2name(hdr, rec->rid), "1") == 0);
  BOOST_CHECK(rec->pos = 1);
  bcf_unpack(rec, BCF_UN_STR);
  BOOST_CHECK(strcmp(rec->d.id, "rs999") == 0);
  BOOST_CHECK(strcmp(rec->d.allele[0], "G") == 0);
  BOOST_CHECK(strcmp(rec->d.allele[1], "T") == 0);
  BOOST_CHECK(strcmp(rec->d.allele[2], "C") == 0);

  char *dst_s = NULL;
  int ndst;
  bcf_get_info_string(hdr, rec, "I1", &dst_s, &ndst);
  BOOST_CHECK(strcmp(dst_s, "1") == 0);
  int32_t *dst_i = NULL;
  ndst = 0;
  bcf_get_info_int32(hdr, rec, "I2", &dst_i, &ndst);
  BOOST_CHECK(ndst == 1 && dst_i[0] == 2);
  float *dst_f = NULL;
  ndst = 0;
  bcf_get_info_float(hdr, rec, "I3", &dst_f, &ndst);
  BOOST_CHECK(ndst == 1 && dst_f[0] == 3.0f);
  // TODO: Make sure memory allocated is freed, especially in call.cc

  dst_f = NULL;
  ndst = 0;
  bcf_get_format_float(hdr, rec, "S1", &dst_f, &ndst);
  BOOST_CHECK(ndst == 3 && dst_f[0] == 1.2f && dst_f[1] == 2.3f && dst_f[2] == 4.5f);

  dst_i = NULL;
  ndst = 0;
  bcf_get_format_int32(hdr, rec, "S2", &dst_i, &ndst);
  for(int a = 0; a < 9; a++)
    BOOST_CHECK(dst_i[a] == (a+1));

  char **dst_str = NULL;
  ndst = 0;
  bcf_get_format_string(hdr, rec, "S3", &dst_str, &ndst);
  BOOST_CHECK(strcmp(dst_str[0], "one,two") == 0);
  BOOST_CHECK(strcmp(dst_str[1], "three,four") == 0);
  BOOST_CHECK(strcmp(dst_str[2], "five,six") == 0);
  
}



BOOST_AUTO_TEST_CASE(Test_VCF_Input)
{
  // Read in a vcf/bcf file and test all the getter methods
  std::string fname = buildTmpHtsFile(false);
  hts::bcf::File vcfr(fname.c_str(), "r");
  std::pair<char**, int> samples = vcfr.samples();
  BOOST_CHECK(samples.second == 3);
  BOOST_CHECK(strcmp(samples.first[0], "NA00001") == 0);
  BOOST_CHECK(strcmp(samples.first[1], "NA00002") == 0);
  BOOST_CHECK(strcmp(samples.first[2], "NA00003") == 0);
  

}
