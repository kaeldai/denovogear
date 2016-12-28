/*
 * Copyright (c) 2016 Reed A. Cartwright
 * Authors:  Reed A. Cartwright <reed@cartwrig.ht>
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

#pragma once
#ifndef DNG_APP_SIM_H
#define DNG_APP_SIM_H

#include <dng/task.h>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

namespace dng {
namespace task {
namespace sim {

// use X-Macros to specify argument variables
struct arg_t : public task::arg_t {
#define XM(lname, sname, desc, type, def) type JS(_,lname) ;
#	include "sim.xmh"
#undef XM
};

inline void add_app_args(po::options_description &desc, arg_t &arg) {
    desc.add_options()
#define XM(lname, sname, desc, type, def) ( \
	XS(lname) IFD(sname, "," BOOST_PP_STRINGIZE sname), \
	po::value< type >(&arg.XV(lname))->default_value(def), \
	desc )
#	include "sim.xmh"
#undef XM
    ;
}

} // namespace sim




class Sim : public Task<sim::arg_t> {
public:
    int operator()(argument_type &arg);
};

}
} // namespace dng::task


#endif
