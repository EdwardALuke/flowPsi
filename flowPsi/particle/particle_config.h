//#############################################################################
//#
//# Copyright 2015, Mississippi State University
//#
//# This file is part of the flowPsi computational fluid dynamics solver.
//#
//# The flowPsi solver is free software: you can redistribute it and/or modify
//# it under the terms of the GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The flowPsi solver is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# GNU General Public License for more details.
//#
//# You should have received a copy of the GNU General Public License
//# along with the flowPsi solver.  If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################

#ifndef PARTICLE_CONFIG_H
#define PARTICLE_CONFIG_H

// this file is mainly used to define global option in compiling
// the particle tracking program

// 1. flag to perform runtime profiling
//#define PROFILING

#ifdef PROFILING
// profiling items
#define TIMING
#define ALGORITHM_STAT
#define FACE_HISTORY_STAT
#endif

// this flag defines whether or not to compute the
// average face cross value each iteration
//#define SHOW_AFC

// this activates some of the runtime consistency
// checking code, should only be used for diagnostic
// purpose as the self checking is slow
//#define SELF_CHECKING

// this flag will activate the code that performs
// a checking on the particle walking routine to
// ensure that the particles are located relatively correctly
//#define CHECK_PARTICLE_WALK

// this flag activates the PAPI instrumentation code
//#define USE_PAPI

// this flag will active the "no split" mode, i.e., the code
// uses a single buffer for the communicated cell data regardless
// of local or remote data
//#define NO_SPLIT_MODE

//this flag will activate the particle mixture mode
#define PARTICLE_MIXTURE

// this flag enables checking for particles crossing a mesh hole
// #define PARTICLE_HOLE_CHECK
// this flag enables some on-screen reports on particles entering a hole
// #define PARTICLE_HOLE_VERBOSE

#endif
