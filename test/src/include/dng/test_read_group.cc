#define BOOST_TEST_MODULE "dng::io::ped"

#include <boost/test/unit_test.hpp>

#include "dng/simulator/data_factory.h"
#include "dng/read_group.h"
#include "dng/hts/bcf.h"
#include "dng/hts/bam.h"

#include <stdio.h>
#include <iostream>
#include <iosfwd>
#include <string>
#include <fstream>

using namespace dng;

#define GROUP_STR(m) (boost::format("%s:LB-%s") % (m)->name % (m)->family_name()).str()
#define LIBRARY_STR(m) (boost::format("%s\tLB-%s") % (m)->name % (m)->family_name()).str()
#define SAM_GROUP_STR(m) std::to_string((m)->mid)
#define SAM_LIBRARY_STR(m) ((m)->name + "\t" + (m)->name + "-" + (m)->family_name())

std::string filename() {
  char fname[] = "pedigree_XXXXXX";
  int fd = mkstemp(fname);
  assert(fd != -1);
  return std::string(fname);
}

BOOST_AUTO_TEST_CASE(Test_VCF)
{
  sim::Factory *fac = sim::Factory::newInstance(sim::VCF);
  sim::Member *child = fac->AddTrio(nullptr, sim::Gender::Female, nullptr, nullptr);
  std::string file = filename();
  fac->publishVCFData(file);

  //vector<hts::File> indata;
  //vector<hts::bcf::File> bcfdata;
  //bcf.emplace_back(file, "r");


  ReadGroups rgs;
  hts::bcf::File vcf(file.c_str(), "r");
  rgs.ParseSamples(vcf);

  ReadGroups::StrSet groups = rgs.groups();
  BOOST_CHECK(groups.size() == 3);
  BOOST_CHECK(groups.find(GROUP_STR(child)) != groups.end());
  BOOST_CHECK(groups.find(GROUP_STR(child->mom)) != groups.end());
  BOOST_CHECK(groups.find(GROUP_STR(child->dad)) != groups.end());

  ReadGroups::StrSet libs = rgs.libraries();
  BOOST_CHECK(libs.size() == 3);
  BOOST_CHECK(libs.find(LIBRARY_STR(child)) != libs.end());
  BOOST_CHECK(libs.find(LIBRARY_STR(child->mom)) != libs.end());
  BOOST_CHECK(libs.find(LIBRARY_STR(child->dad)) != libs.end());

  ReadGroups::StrSet sm = rgs.samples();
  BOOST_CHECK(sm.size() == 3);
  BOOST_CHECK(sm.find(child->name) != sm.end());
  BOOST_CHECK(sm.find(child->mom->name) != sm.end());
  BOOST_CHECK(sm.find(child->dad->name) != sm.end());


  /*
  for(std::string str : rgs.groups()) {
    std::cout << str << std::endl;
  }
  
  for(std::string str : rgs.libraries()) {
    std::cout << str << std::endl;

  }
  for(std::string str : rgs.samples()) {
    std::cout << str << std::endl;
  }
  */

  //remove(file.c_str());
}

BOOST_AUTO_TEST_CASE(Test_BCF)
{


}

BOOST_AUTO_TEST_CASE(Test_SAM)
{
  sim::Factory *fac = sim::Factory::newInstance(sim::SAM);
  sim::Member *child = fac->AddTrio(nullptr, sim::Gender::Female, nullptr, nullptr);
  std::string file = filename();
  fac->publishSAMData(file);

  std::vector<hts::bam::File> bamdata;
  bamdata.emplace_back(file.c_str(), "r");
  ReadGroups rgs;
  rgs.ParseHeaderText(bamdata);

  ReadGroups::StrSet groups = rgs.groups();
  BOOST_CHECK(groups.size() == 3);
  BOOST_CHECK(groups.find(SAM_GROUP_STR(child)) != groups.end());
  BOOST_CHECK(groups.find(SAM_GROUP_STR(child->mom)) != groups.end());
  BOOST_CHECK(groups.find(SAM_GROUP_STR(child->dad)) != groups.end());

  ReadGroups::StrSet libs = rgs.libraries();
  BOOST_CHECK(libs.size() == 3);
  BOOST_CHECK(libs.find(SAM_LIBRARY_STR(child)) != libs.end());
  BOOST_CHECK(libs.find(SAM_LIBRARY_STR(child->mom)) != libs.end());
  BOOST_CHECK(libs.find(SAM_LIBRARY_STR(child->dad)) != libs.end());

  ReadGroups::StrSet sm = rgs.samples();
  BOOST_CHECK(sm.size() == 3);
  BOOST_CHECK(sm.find(child->name) != sm.end());
  BOOST_CHECK(sm.find(child->mom->name) != sm.end());
  BOOST_CHECK(sm.find(child->dad->name) != sm.end());

  /*

  for(std::string str : rgs.groups()) {
    std::cout << str << std::endl;
  }
  
  for(std::string str : rgs.libraries()) {
    std::cout << str << std::endl;

  }
  for(std::string str : rgs.samples()) {
    std::cout << str << std::endl;
  }
  */

}

BOOST_AUTO_TEST_CASE(Test_BAM)
{


}


