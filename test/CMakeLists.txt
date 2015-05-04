find_package(Boost COMPONENTS unit_test_framework REQUIRED)
find_package(HTSLIB)

include_directories(../src/include)
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${EIGEN3_INCLUDE_DIR}/unsupported)
add_definitions(-DBOOST_TEST_DYN_LINK)

file(GLOB_RECURSE TESTS test_*[cc|cpp])
foreach(test ${TESTS})
  get_filename_component(testName ${test} NAME_WE)
  add_executable(${testName} ${test})
  set_target_properties(${testName} PROPERTIES CXX_STANDARD 11)
  set_target_properties(${testName} PROPERTIES CXX_STANDARD_REQUIRED ON)
  target_link_libraries(${testName} ${Boost_LIBRARIES} ${HTSLIB_LIBRARIES})
  add_test(${testName} ${testName})
endforeach(test)


