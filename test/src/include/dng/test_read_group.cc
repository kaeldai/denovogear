#define BOOST_TEST_MODULE "dng::io::ped"

#include <boost/test/unit_test.hpp>

#include "dng/simulator/models.h"
#include "dng/read_group.h"
#include "dng/hts/bcf.h"
#include "dng/hts/bam.h"

#include <stdio.h>
#include <iostream>
#include <iosfwd>
#include <string>
#include <fstream>

// Helper function for test cases that use VCF header as input
void testVCF(dng::sim::SimBuilder &model) {

  // Create temporary VCF file
  char fname[] = "vcf_XXXXXX";
  int fd = mkstemp(fname);
  std::string vcffile = std::string(fname);
  model.publishData(vcffile, dng::sim::SeqFormat::VCF);

  // Parse in VCF file to ReadGroups
  dng::ReadGroups rgs;
  hts::bcf::File vcf(vcffile.c_str(), "r");
  rgs.ParseSamples(vcf);

  // Compare simulated data members/libraries to what was parsed by ReadGroups.
  size_t lib_count = 0;
  for(dng::sim::Member *m : model.getMembers()) {
    for(dng::sim::Sample lib : m->libraries) {
      lib_count++;
      dng::detail::DataBase::index<dng::rg::id>::type::iterator db = rgs.data().get<dng::rg::id>().find(lib.id_vcf());
      BOOST_CHECK(db != rgs.data().get<dng::rg::id>().end());
      BOOST_CHECK(db->sample == lib.sm);
      BOOST_CHECK(db->library == lib.id_vcf()); // ParseSamples uses the "id" for "lb" tag in VCF files.
    }
  }

  BOOST_CHECK(rgs.samples().size() == model.getMembers().size());
  BOOST_CHECK(rgs.groups().size() == lib_count);
  BOOST_CHECK(rgs.libraries().size() == lib_count);

  // Remove temp vcf file
  remove(vcffile.c_str());
}

// Helper function for SAM/BAM input headers
void testSAM(dng::sim::SimBuilder &model)
{

  // Create temporary file
  char fname[] = "sam_XXXXXX";
  int fd = mkstemp(fname);
  std::string samfile = std::string(fname);
  model.publishData(samfile, dng::sim::SeqFormat::SAM);

  // Parse SAM file into ReadGroups
  dng::ReadGroups rgs;
  std::vector<hts::bam::File> bamdata;
  bamdata.emplace_back(samfile.c_str(), "r");
  rgs.ParseHeaderText(bamdata);

  // Check data
  size_t lib_count = 0;
  for(dng::sim::Member *m : model.getMembers()) {
    for(dng::sim::Sample lib : m->libraries) {
      lib_count++;
      
      dng::detail::DataBase::index<dng::rg::id>::type::iterator db = rgs.data().get<dng::rg::id>().find(lib.id_sam());
      BOOST_CHECK(db != rgs.data().get<dng::rg::id>().end());
      BOOST_CHECK(db->sample == lib.sm);
      BOOST_CHECK(db->library == lib.id_sam());//lib.sm + "-" + lib.lb);
    }
  }

  BOOST_CHECK(rgs.samples().size() == model.getMembers().size());
  BOOST_CHECK(rgs.groups().size() == lib_count);
  BOOST_CHECK(rgs.libraries().size() == lib_count);

  // Delete SAM file 
  remove(samfile.c_str());
}

  
BOOST_AUTO_TEST_CASE(Test_Trio_VCF)
{
  dng::sim::Trio trio;
  testVCF(trio);
}

BOOST_AUTO_TEST_CASE(Test_Trio_SAM)
{
  dng::sim::Trio trio;
  testSAM(trio);
}

BOOST_AUTO_TEST_CASE(Test_FamilyTree_VCF)
{
  dng::sim::FamilyTree ft;
  testVCF(ft);
}

BOOST_AUTO_TEST_CASE(Test_FamilyTree_SAM)
{
  dng::sim::FamilyTree ft;
  testSAM(ft);
}


BOOST_AUTO_TEST_CASE(Test_ExtendedFamily_VCF)
{
  dng::sim::ExtendedTree et;
  testVCF(et);
}

BOOST_AUTO_TEST_CASE(Test_ExtendedFamily_SAM)
{
  dng::sim::ExtendedTree et;
  testSAM(et);
}


