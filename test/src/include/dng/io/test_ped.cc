#define BOOST_TEST_MODULE "dng::io::ped"

#include <boost/test/unit_test.hpp>
#include "dng/io/ped.h"

#include "dng/simulator/models.h"
#include <fstream>
#include <iterator>
#include <iosfwd>
#include <iostream>


// Copied from call.cc
template<class Elem, class Traits> inline
boost::iterator_range<std::istreambuf_iterator<Elem, Traits> >
istreambuf_range(std::basic_istream<Elem, Traits> &in) {
    return boost::iterator_range<std::istreambuf_iterator<Elem, Traits>>(
               std::istreambuf_iterator<Elem, Traits>(in),
               std::istreambuf_iterator<Elem, Traits>());
}


/*
 * Compare the input, dng::sim::Member, with the dng::io::Pedigree::Member created from reading in the
 * the pedigree file.
 */
void confirm(dng::sim::Member *orig, dng::io::Pedigree &ped) {
  size_t found = 0;
  for(dng::io::Pedigree::Member m : ped.table()) {
    if(ped.name(m.child) == orig->id) {
      BOOST_CHECK(orig->mom_ptr == nullptr ? ped.name(m.mom) == "" : orig->mom_ptr->id == ped.name(m.mom)); 
      BOOST_CHECK(orig->dad_ptr == nullptr ? ped.name(m.dad) == "" : orig->dad_ptr->id == ped.name(m.dad)); 
      // Family names are not yet implements in ped.h
      // BOOST_CHECK(orig->family_name() == ped.name(m.fam));
      BOOST_CHECK(orig->sex == m.sex);
      found++;
    }
  }

  // Make sure exactly one instance of indiviual was found
  BOOST_CHECK(found == 1);

}


void run_test(dng::sim::SimBuilder &sim)
{
  // Use the simulator to create a pedigree and write to a temporary file
  //dng::sim::ExtendedTree trio;
  char fname[] = "pedigree_XXXXXX";
  int fd = mkstemp(fname);
  std::ofstream ped_io;
  ped_io.open(fname);
  sim.publishPed(ped_io);
  ped_io.flush();
  ped_io.close();

  // Mimic the way call.cc opens the pedigree file
  dng::io::Pedigree ped;
  std::ifstream ped_file(fname);
  ped.Parse(istreambuf_range(ped_file));

  // Check that pedigree is correct
  const dng::io::Pedigree::MemberTable family = ped.table();
  std::vector<std::shared_ptr<dng::sim::Member>> sim_fam = sim.GetPedigree();
  BOOST_CHECK(family.size() == (sim_fam.size()+1)); // Check ped size
  for(std::shared_ptr<dng::sim::Member> m : sim_fam) {
    confirm(m.get(), ped); // compares each member with ped with the original data from the simulator
  }

  // delete the temporary ped file
  remove(fname);
}

// Parents and child
BOOST_AUTO_TEST_CASE(Test_Trio)
{
  dng::sim::Trio trio;
  run_test(trio);
}

// Parents, child and grandparents
BOOST_AUTO_TEST_CASE(Test_FamilyTree)
{
  dng::sim::FamilyTree ft;
  run_test(ft);
}


// Parents, two child, and the female child has a daughter
BOOST_AUTO_TEST_CASE(Test_ExtendedFamily)
{
  dng::sim::ExtendedTree et;
  run_test(et);
}

