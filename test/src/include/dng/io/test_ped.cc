#define BOOST_TEST_MODULE "dng::io::ped"

#include <boost/test/unit_test.hpp>
#include "dng/io/ped.h"

#include "dng/simulator/data_factory.h"
#include <fstream>
#include <iterator>
#include <iosfwd>
#include <iostream>
#include <boost/tokenizer.hpp>

using namespace dng::sim;

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
    if(ped.name(m.child) == orig->name) {
      BOOST_CHECK(orig->mom == nullptr ? ped.name(m.mom) == "" : orig->mom->name == ped.name(m.mom)); 
      BOOST_CHECK(orig->dad == nullptr ? ped.name(m.dad) == "" : orig->dad->name == ped.name(m.dad)); 
      // Family names are not yet implements in ped.h
      // BOOST_CHECK(orig->family_name() == ped.name(m.fam));
      BOOST_CHECK(orig->sex == m.sex);
      found++;
    }
  }

  // Make sure exactly one instance of indiviual was found
  BOOST_CHECK(found == 1);

}


/**
 * Creates a temporary pedigree file for ped.h input
 */
std::string create_p(dng::sim::Factory *fac) {
  char fname[] = "pedigree_XXXXXX";
  int fd = mkstemp(fname);
  assert(fd != -1);

  std::ofstream ped_io;
  ped_io.open(fname);
  fac->publishPed(ped_io);
  ped_io.flush();
  ped_io.close();

  return std::string(fname);
}

/**
 * Cleans up the pedigree file
 */
void clean_p(std::string &filename) {
  remove(filename.c_str());
}


BOOST_AUTO_TEST_CASE(Test_Trio)
{
  // Build a single trio
  dng::sim::Factory *sb = dng::sim::Factory::newInstance();
  dng::sim::Member *child = sb->AddTrio(nullptr, dng::sim::Gender::Male, nullptr, nullptr);
 
  // Create a pedigree file and pass it into dng::io::Pedigree
  // TODO: Find a way to pass in the output stream directly instead of writing to a file?
  std::string file_name = create_p(sb);
  dng::io::Pedigree ped;
  std::ifstream ped_file(file_name);
  ped.Parse(istreambuf_range(ped_file));

  // Pedigree should have 3 members plus the extra 0 root
  const dng::io::Pedigree::MemberTable family = ped.table();
  BOOST_CHECK(family.size() == (3+1));
  
  confirm(child, ped); // Check child
  confirm(child->mom, ped); // Check dad
  confirm(child->dad, ped); // Check mom
  
  clean_p(file_name);
}


BOOST_AUTO_TEST_CASE(Test_Tree)
{
  
  // Build an extend family; child, parents, and all grandparents
  dng::sim::Factory *sb = dng::sim::Factory::newInstance();
  dng::sim::Member *child = sb->AddTrio(nullptr, dng::sim::Gender::Male, nullptr, nullptr);
  sb->AddTrio(child->mom, nullptr, nullptr); // materal grandparents
  sb->AddTrio(child->dad, nullptr, nullptr); // paternal grandparents

  std::string file_name = create_p(sb);
  dng::io::Pedigree ped;
  std::ifstream ped_file(file_name);
  ped.Parse(istreambuf_range(ped_file));

  const dng::io::Pedigree::MemberTable family = ped.table();
  BOOST_CHECK(family.size() == (7+1));

  confirm(child, ped); // child
  confirm(child->mom, ped); // mom
  confirm(child->dad, ped); // dad

  // Grandparents
  confirm(child->mom->mom, ped); 
  confirm(child->mom->dad, ped); 
  confirm(child->dad->mom, ped); 
  confirm(child->dad->dad, ped); 

  clean_p(file_name);
}


BOOST_AUTO_TEST_CASE(Test_Branched_Ped)
{
  // Dad has two children with two mothers. The other mother's grandparents are not in pedigree
  dng::sim::Factory *sb = dng::sim::Factory::newInstance();
  dng::sim::Member *child1 = sb->AddTrio(nullptr, dng::sim::Gender::Male, nullptr, nullptr);
  sb->AddTrio(child1->mom, nullptr, nullptr); // materal grandparents
  sb->AddTrio(child1->dad, nullptr, nullptr); // paternal grandparents
  dng::sim::Member *child2 = sb->AddTrio(nullptr, dng::sim::Gender::Female, nullptr, child1->dad);

  std::string file_name = create_p(sb);
  dng::io::Pedigree ped;
  std::ifstream ped_file(file_name);
  ped.Parse(istreambuf_range(ped_file));

  const dng::io::Pedigree::MemberTable family = ped.table();
  BOOST_CHECK(family.size() == (9+1));

  confirm(child1, ped); // child
  confirm(child1->mom, ped); // mom
  confirm(child1->dad, ped); // dad
  confirm(child1->mom->mom, ped); 
  confirm(child1->mom->dad, ped); 
  confirm(child1->dad->mom, ped); 
  confirm(child1->dad->dad, ped); 

  confirm(child2, ped);
  confirm(child2->dad, ped);
  confirm(child2->mom, ped);

  clean_p(file_name);
}

BOOST_AUTO_TEST_CASE(Test_Multiple_Family)
{
  // At the moment ped.h does not handle multiple families in the same pedigree
}
