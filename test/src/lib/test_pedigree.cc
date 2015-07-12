#define BOOST_TEST_MODULE "pedigree.cc"

#include <boost/test/unit_test.hpp>
#include "dng/pedigree.h"
#include "dng/io/ped.h"
#include "dng/read_group.h"
#include "dng/hts/bcf.h"

#include <iostream>
#include <fstream>
#include <iosfwd>
#include <iterator>
#include <stdio.h>

#include "htslib/sam.h"

#include <boost/range/iterator_range.hpp>
#include "dng/simulator/data_factory.h"

std::string filename() {
  char fname[] = "pedigree_XXXXXX";
  int fd = mkstemp(fname);
  assert(fd != -1);
  return std::string(fname);
}

template<class Elem, class Traits> inline
boost::iterator_range<std::istreambuf_iterator<Elem, Traits> >
istreambuf_range(std::basic_istream<Elem, Traits> &in) {
    return boost::iterator_range<std::istreambuf_iterator<Elem, Traits>>(
               std::istreambuf_iterator<Elem, Traits>(in),
               std::istreambuf_iterator<Elem, Traits>());
}




BOOST_AUTO_TEST_CASE(Test_Trio_VCF)
{
  dng::sim::Factory *fac = dng::sim::Factory::newInstance(dng::sim::VCF);
  dng::sim::Member *child = fac->AddTrio(nullptr, dng::sim::Gender::Female, nullptr, nullptr);
  std::string file = filename();
  std::string ped_file = file + ".ped";
  std::string vcf_file = file + ".vcf";

  {
    std::ofstream fp;
    fp.open(ped_file);
    fac->publishPed(fp);
    fp.flush();
    fp.close();
    
    fac->publishVCFData(file);
    fac->publishData();
  }

  // theta, mu, mu-somatic, mu-library, ref-weigth, nuc-freq
  dng::Pedigree::params_t params = {0.001, 1e-9, 0.0, 0.0, 1.0, std::array<double, 4>{0.3, 0.2, 0.2, 0.3}};

  dng::io::Pedigree ped;
  std::ifstream ped_stream(ped_file);
  ped.Parse(istreambuf_range(ped_stream));
	    
  dng::ReadGroups rgs;
  hts::bcf::File vcf(file.c_str(), "r");
  rgs.ParseSamples(vcf);

  

  dng::Pedigree peeler;
  //peeler.Initialize(params);
  //peeler.Construct(ped, rgs);

  



  //peeler.library_lower(u);
  //double d = peeler.CalculateLogLikelihood(ref_index) + scale;
  //double p = peeler.CalculateMutProbability(ref_index);
  //double s = peeler.LogPeelAll();


}


/*
BOOST_AUTO_TEST_CASE(Test_Pedigree)
{
  dng::io::Pedigree ped;
  std::string fname = getPedName();
  buildPedigree(fname, 3);
  std::ifstream ped_in(fname);
  ped.Parse(istreambuf_range(ped_in));
  //std::cout << ped.member_count() << std::endl;
  const dng::io::Pedigree::MemberTable family = ped.table();
  for(dng::io::Pedigree::Member member : family) {
    std::cout << "<" << member.sample_tree << ">" << std::endl;
    //std::cout << "family = " << member.fam << std::endl;
    std::cout << "\tchild = " <<member.child << std::endl;
    std::cout << "\tmom = " << member.mom << std::endl;
    std::cout << "\tdad = " << member.dad << std::endl;
    std::cout << "\tfam = " << member.fam << std::endl;
    //std::cout << "\tsex = " << member.sex << std::endl;
    //std::cout << "id = " << ped.id(member.sample_tree) << std::endl;
 
  }
  //boost::iterator_range<std::istreambuf_iterator< ped_range = boost::iterator_range<std::istreambuf_iterator<Elem, Traits>>(std::istreambuf_iterator<Elem, Traits>(ped_in), std::istreambuf_iterator<Elem, Traits>());
  //ped.Parse(ped_range);

}
*/



