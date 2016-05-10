## import the required files

import ProcessExporter

## Define a typical error for the plugin
class MoMEMta_Error(Exception): pass

## Does it define a new interface (will be avaible with ./bin/mg5_aMC --mode=maddm
## Put None if no dedicated command are required
new_interface = False
 
## Does it define a new output mode. Need to be define
new_output = { 'MoMEMta': ProcessExporter.ProcessExporterMoMEMta }

## The test/code have been validated up to this version
latest_validated_version = '2.3.4'

