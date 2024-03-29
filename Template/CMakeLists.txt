# This file was automatically generated by MoMEMta-MaGMEE,
# A MadGraph Matrix Element Exporter plugin for MoMEMta.
#
# It is subject to MoMEMta-MaGMEE's license and copyright:
#
# Copyright (C) 2016  Universite catholique de Louvain (UCL), Belgium
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

cmake_minimum_required(VERSION 2.8...3.20)

if(${CMAKE_VERSION} VERSION_LESS 3.12)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

project(%(dir_name)s CXX)

# Find MoMEMta
find_package(MoMEMta CONFIG REQUIRED)

# Flags necessary to ensure complex arithmetic performances on-par with
# Fortran code:
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-variable -std=c++11 -O3 -fcx-fortran-rules -fcx-limited-range")

file(GLOB_RECURSE SOURCES "*.cc")

include_directories("include")
%(include_subprocs_list)s

add_library(me_%(dir_name)s SHARED ${SOURCES})
target_link_libraries(me_%(dir_name)s momemta::momemta)
