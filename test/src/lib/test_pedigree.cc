#define BOOST_TEST_MODULE "pedigree.cc"

#include "dng/pedigree.h"
#include "dng/io/ped.h"
#include "dng/read_group.h"
#include "dng/hts/bcf.h"
#include "dng/hts/bam.h"
#include "dng/matrix.h"
#include "dng/mutation.h"

#include "dng/simulator/simulator.h"
#include "dng/simulator/models.h"
#include <boost/test/unit_test.hpp>


#include <iostream>
#include <fstream>
#include <iosfwd>
#include <iterator>
#include <stdio.h>


#include "htslib/sam.h"

#include <boost/range/iterator_range.hpp>

// Copied pop priors based on an A reference w/ 20-30-30-20 percentages.
#define PRIORS 0.9989511188,0.0001997602597,0.0001997602597,0.0002996403896,9.987014485e-05,3.994006993e-08,5.991010489e-08,9.987014485e-05,5.991010489e-08,0.0001498201948

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

std::string create_ped(dng::sim::SimBuilder &sim) {
  char fname[] = "pedigree_XXXXXX";
  int fd = mkstemp(fname);
  std::ofstream ped_io;
  ped_io.open(fname);
  sim.publishPed(ped_io);
  ped_io.flush();
  ped_io.close();
  return std::string(fname);
}

std::string create_sam(dng::sim::SimBuilder &sim) {
  char fname[] = "sam_XXXXXX";
  int fd = mkstemp(fname);
  std::string samfile = std::string(fname);
  sim.publishData(samfile, dng::sim::SeqFormat::SAM);
  return samfile;
}

// Attempts to create non-mutational transition matricies for the tree. Using non-mutational matricies to 
// test the the peeling functions since they are the simplest - may change in the future.
dng::TransitionVector create_transition_matricies(dng::Pedigree &ped) {
  dng::TransitionVector tv;
  for(dng::Pedigree::transition_t t : ped.transitions()) {
    if(t.type == dng::Pedigree::TransitionType::Germline) {
      dng::MutationMatrix m1 = dng::f81::matrix(t.length1, {0.3, 0.2, 0.2, 0.3});
      dng::MutationMatrix m2 = dng::f81::matrix(t.length2, {0.3, 0.2, 0.2, 0.3});
      tv.push_back(dng::meiosis_diploid_matrix(m1, m2, 0));
    } else if(t.type == dng::Pedigree::TransitionType::Library || (t.type == dng::Pedigree::TransitionType::Somatic)) {
      dng::MutationMatrix m1 = dng::f81::matrix(t.length1, {0.3, 0.2, 0.2, 0.3});
      tv.emplace_back(dng::mitosis_diploid_matrix(m1, 0));
    } else {
      tv.push_back({});
    }
  }
  return tv;
}



// Check the appropiate labels exists
void check_labels(dng::Pedigree &ped, dng::sim::SimBuilder &sim) { 
  std::vector<std::string> labels = ped.labels();
  std::vector<std::shared_ptr<dng::sim::Member>> sim_fam = sim.GetPedigree();
  for(std::shared_ptr<dng::sim::Member> m : sim_fam) {
    for(dng::sim::Library lib : m->libraries) {
      // A label should exists for each unique library
      BOOST_CHECK(std::find(labels.begin(), labels.end(), "LB-"+lib.id) != labels.end());
    }
    
    // If individual has mom and/or dad then check that GL node exits
    if(m->mom_ptr != nullptr) {
      BOOST_CHECK(std::find(labels.begin(), labels.end(), "GL-"+m->mom_ptr->id) != labels.end());
    }
    if(m->dad_ptr != nullptr) {
      BOOST_CHECK(std::find(labels.begin(), labels.end(), "GL-"+m->dad_ptr->id) != labels.end());
    }      
  }
}





BOOST_AUTO_TEST_CASE(Test_Trio)
{
  dng::sim::Trio trio;
  std::string pfile = create_ped(trio);
  std::string sfile = create_sam(trio);

  dng::io::Pedigree ped;
  std::ifstream ped_file(pfile);
  ped.Parse(istreambuf_range(ped_file));
    
  dng::ReadGroups rgs;
  std::vector<hts::bam::File> bamdata;
  bamdata.emplace_back(sfile.c_str(), "r");
  rgs.ParseHeaderText(bamdata);

  dng::Pedigree pedigree;
  pedigree.Construct(ped, rgs, 1e-8, 1e-8, 1e-8);

  BOOST_CHECK(pedigree.num_nodes() == 5);
  BOOST_CHECK(pedigree.library_nodes().second == 5);
  check_labels(pedigree, trio);

  dng::peel::workspace_t ws = pedigree.CreateWorkspace();
  dng::TransitionVector tv = create_transition_matricies(pedigree);
  dng::GenotypeArray pp{10};
  pp << PRIORS;
  ws.SetFounders(pp);
  BOOST_CHECK_CLOSE(pedigree.PeelForwards(ws, tv), -1.32439989e-07, 1);
  BOOST_CHECK_CLOSE(pedigree.PeelBackwards(ws, tv), -1.32439989e-07, 1);

  remove(pfile.c_str());
  remove(sfile.c_str());
}

BOOST_AUTO_TEST_CASE(Test_FamilyTree)
{
  dng::sim::FamilyTree ft;
  std::string pfile = create_ped(ft);
  std::string sfile = create_sam(ft);

  dng::io::Pedigree ped;
  std::ifstream ped_file(pfile);
  ped.Parse(istreambuf_range(ped_file));
    
  dng::ReadGroups rgs;
  std::vector<hts::bam::File> bamdata;
  bamdata.emplace_back(sfile.c_str(), "r");
  rgs.ParseHeaderText(bamdata);

  dng::Pedigree pedigree;
  pedigree.Construct(ped, rgs, 1e-8, 1e-8, 1e-8);

  // TODO: double-check the math. Maybe reduce the number of libraries.
  BOOST_CHECK(pedigree.num_nodes() == 27);
  BOOST_CHECK(pedigree.library_nodes().second == 27);
  check_labels(pedigree, ft);

  dng::peel::workspace_t ws = pedigree.CreateWorkspace();
  dng::TransitionVector tv = create_transition_matricies(pedigree);
  dng::GenotypeArray pp{10};
  pp << PRIORS;
  ws.SetFounders(pp);

  BOOST_CHECK_CLOSE(pedigree.PeelForwards(ws, tv), -4.5426485e-07, 1);
  BOOST_CHECK_CLOSE(pedigree.PeelBackwards(ws, tv), -4.5426485e-07, 1);

  remove(pfile.c_str());
  remove(sfile.c_str());
}


BOOST_AUTO_TEST_CASE(Test_ExtendedTree)
{
  dng::sim::ExtendedTree et;
  std::string pfile = create_ped(et);
  std::string sfile = create_sam(et);

  dng::io::Pedigree ped;
  std::ifstream ped_file(pfile);
  ped.Parse(istreambuf_range(ped_file));
    
  dng::ReadGroups rgs;
  std::vector<hts::bam::File> bamdata;
  bamdata.emplace_back(sfile.c_str(), "r");
  rgs.ParseHeaderText(bamdata);

  dng::Pedigree pedigree;
  pedigree.Construct(ped, rgs, 1e-8, 1e-8, 1e-8);

  BOOST_CHECK(pedigree.num_nodes() == 17);
  BOOST_CHECK(pedigree.library_nodes().second == 17);
  check_labels(pedigree, et);

  dng::peel::workspace_t ws = pedigree.CreateWorkspace();
  dng::TransitionVector tv = create_transition_matricies(pedigree);
  dng::GenotypeArray pp{10};
  pp << PRIORS;
  ws.SetFounders(pp);
  BOOST_CHECK_CLOSE(pedigree.PeelForwards(ws, tv), -3.5961863e-07, 1);
  BOOST_CHECK_CLOSE(pedigree.PeelBackwards(ws, tv),-3.5961863e-07, 1);

  remove(pfile.c_str());
  remove(sfile.c_str());
}


