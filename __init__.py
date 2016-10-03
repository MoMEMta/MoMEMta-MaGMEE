# MoMEMta-MaGMEE: a MadGraph Matrix Element Exporter plugin for MoMEMta
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

## import the required files

import ProcessExporter

## Define a typical error for the plugin
class MoMEMta_Error(Exception): pass

## We do not modify the MadGraph interface
new_interface = False
 
## We add a new output class 
new_output = { 'MoMEMta_standalone': ProcessExporter.ProcessExporterMoMEMta }

## Plugin version
__version__ = (0,1,0)

## MG5 version requirements
minimal_mg5amcnlo_version = (2,3,4) 
maximal_mg5amcnlo_version = (1000,1000,1000)
latest_validated_version = (2,4,0)
