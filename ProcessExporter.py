## For the present version of this code,
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
#
## For the original version of this code,
# Copyright (c) 2009 The MadGraph5_aMC@NLO Development team and Contributors

"""Methods and classes to export models and matrix elements to MoMEMta-dedicated C++ format."""

import fractions
import itertools
import logging
import os
import re

import madgraph.core.base_objects as base_objects
import madgraph.core.color_algebra as color
import madgraph.core.helas_objects as helas_objects
import madgraph.iolibs.helas_call_writers as helas_call_writers
import madgraph.iolibs.file_writers as writers
import madgraph.iolibs.group_subprocs as group_subprocs
from madgraph import MadGraph5Error, MG5DIR
from madgraph.iolibs.export_cpp import UFOModelConverterCPP
from madgraph.iolibs.export_v4 import VirtualExporter

_file_path = os.path.join(MG5DIR, "PLUGIN", "MoMEMta-MaGMEE")
_template_dir = os.path.join(_file_path, "Template")
logger = logging.getLogger('madgraph.export_pythia8')


#===============================================================================
# ProcessExporterCPP
#===============================================================================
class OneProcessExporterMoMEMta(object):
    """Class to take care of exporting a set of matrix elements to MoMEMta-dedicated C++ format."""

    # Static variables (for inheritance)
    process_template_h                 = 'process.h'
    process_template_cc                = 'process.cc'
    process_class_template             = 'process_class.inc'
    process_definition_template        = 'function_definitions.inc'
    process_wavefunction_template      = 'wavefunctions.inc'
    process_matrix_averaging_template  = 'matrix_averaging.inc'
    single_process_template            = 'matrix.inc'

    class OneProcessExporterMoMEMtaError(Exception):
        pass
    
    def __init__(self, subproc_group, helicity_model, parent_folder, namespace):
    
        # Path to parent folder where the whole process is exported
        self.parent_folder = parent_folder

        # We want the leptons to be split, but still be inside the same class,
        # so we do the splitting here.
        # What we get is an instance of 'HelasMatrixElementList', which is just like a 
        # python list of 'HelasMatrixElement' objects
        
        if isinstance(subproc_group, group_subprocs.SubProcessGroupList):
            self.matrix_elements = subproc_group.split_lepton_grouping().get_matrix_elements()
        
        elif isinstance(subproc_group, group_subprocs.SubProcessGroup):
            temp_group = group_subprocs.SubProcessGroupList([subproc_group])
            self.matrix_elements = temp_group.split_lepton_grouping().get_matrix_elements()
        
        else:
            raise base_objects.PhysicsObject.PhysicsObjectError, "Wrong object type for subproc_group"

        self.processes = sum([me.get('processes') for \
                              me in self.matrix_elements], [])
        self.processes.extend(sum([me.get_mirror_processes() for \
                              me in self.matrix_elements], []))

        self.nprocesses = len(self.matrix_elements)
        self.nprocesses += sum([1 for me in self.matrix_elements if me.get('has_mirror_process')])

        self.process_string = self.processes[0].base_string()
        self.process_number = self.processes[0].get('id')

        # Retrieve model, set model name
        self.model = self.matrix_elements[0].get('processes')[0].get('model')
        self.model_name = OneProcessExporterMoMEMta.get_model_name(self.model.get('name'))

        # Namespace for the C++ code, ie the basename of the parent folder _ the model name
        self.namespace = namespace + "_" + self.model_name

        # Process name: string of type "pp_ttx_..."
        self.process_name = self.get_process_name()

        # C++ class for the whole process
        self.process_class = "P%d_%s" % (self.process_number, self.process_name)

        # Directory where all the files will be written
        self.path = os.path.join(parent_folder, 'SubProcesses', self.process_class)

        self.helas_call_writer = helicity_model

        if not isinstance(self.helas_call_writer, helas_call_writers.CPPUFOHelasCallWriter):
            raise self.OneProcessExporterMoMEMtaError, \
                "helas_call_writer not CPPUFOHelasCallWriter"

        self.nexternal, self.ninitial = \
                        self.matrix_elements[0].get_nexternal_ninitial()
        self.nfinal = self.nexternal - self.ninitial

        # Check if we can use the same helicities for all matrix elements
        hel_matrix = self.get_helicity_matrix(self.matrix_elements[0])
        for me in self.matrix_elements[1:]:
            if self.get_helicity_matrix(me) != hel_matrix:
                raise Exception('Multiple helicities are not supported in mode standalone_cpp_mem.')
            
        # Since all processes have the same helicity structure, this
        # allows us to reuse the same wavefunctions for the
        # different processes
        
        self.wavefunctions = []
        wf_number = 0

        for me in self.matrix_elements:
            for iwf, wf in enumerate(me.get_all_wavefunctions()):
                try:
                    old_wf = \
                           self.wavefunctions[self.wavefunctions.index(wf)]
                    wf.set('number', old_wf.get('number'))
                except ValueError:
                    wf_number += 1
                    wf.set('number', wf_number)
                    self.wavefunctions.append(wf)

        # Also combine amplitudes
        self.amplitudes = helas_objects.HelasAmplitudeList()
        amp_number = 0
        for me in self.matrix_elements:
            for iamp, amp in enumerate(me.get_all_amplitudes()):
                try:
                    old_amp = \
                           self.amplitudes[self.amplitudes.index(amp)]
                    amp.set('number', old_amp.get('number'))
                except ValueError:
                    amp_number += 1
                    amp.set('number', amp_number)
                    self.amplitudes.append(amp)
        diagram = helas_objects.HelasDiagram({'amplitudes': self.amplitudes})
        self.amplitudes = helas_objects.HelasMatrixElement({\
            'diagrams': helas_objects.HelasDiagramList([diagram])})
    
    #===============================================================================
    # Global helper methods
    #===============================================================================
    @classmethod
    def read_template_file(cls, filename, classpath=False):
        """Open a template file and return the contents."""
         
        if isinstance(filename, tuple):
            file_path = filename[0]
            filename = filename[1]
        elif isinstance(filename, str):
            if classpath:
                file_path = cls.__template_path
            else:
                file_path = cls.template_path
        else:
            raise MadGraph5Error('Argument should be string or tuple.')
        
        return open(os.path.join(file_path, filename)).read()

    # Methods for generation of process files for C++

    def generate_process_files(self):
        """Generate the .h and .cc files containing the matrix elements"""

        # Create the files
        filename = os.path.join(self.path, '%s.h' % self.process_class)
        
        self.write_process_h_file(writers.CPPWriter(filename))

        filename = os.path.join(self.path, '%s.cc' % self.process_class)

        self.write_process_cc_file(writers.CPPWriter(filename))

        logger.info('Created files %(process)s.h and %(process)s.cc in %(dir)s' % \
                    {'process': self.process_class, 'dir': self.path})



    #===========================================================================
    # write_process_h_file
    #===========================================================================
    def write_process_h_file(self, writer):
        """Write the class definition (.h) file for the process"""
        
        if not isinstance(writer, writers.CPPWriter):
            raise writers.CPPWriter.CPPWriterError(\
                "writer not CPPWriter")

        replace_dict = {}

        replace_dict['namespace'] = self.namespace
        replace_dict['model_name'] = self.model_name
        replace_dict['process_class_definition'] = self.get_process_class_definition()
        
        # Extract process info lines for all processes
        process_lines = "\n".join([self.get_process_info_lines(me) for me in \
                                   self.matrix_elements])
        replace_dict['process_lines'] = process_lines

        file_content = self.read_template_file((_template_dir, self.process_template_h)) % replace_dict

        # Write the file
        writer.writelines(file_content)

    #===========================================================================
    # write_process_cc_file
    #===========================================================================
    def write_process_cc_file(self, writer):
        """Write the class member definition (.cc) file for the process
        described by matrix_element"""

        if not isinstance(writer, writers.CPPWriter):
            raise writers.CPPWriter.CPPWriterError(\
                "writer not CPPWriter")

        replace_dict = {}

        replace_dict['namespace'] = self.namespace
        replace_dict['process_class'] = self.process_class
        replace_dict['model_name'] = self.model_name
        replace_dict['process_function_definitions'] = self.get_process_function_definitions()

        # Extract process info lines
        replace_dict['process_lines'] = \
                             "\n".join([self.get_process_info_lines(me) for \
                                        me in self.matrix_elements])

        file_content = self.read_template_file((_template_dir, self.process_template_cc)) % replace_dict

        # Write the file
        writer.writelines(file_content)

    #===========================================================================
    # Process export helper functions
    #===========================================================================
    def get_process_class_definition(self):
        """Template values for the class definition in the header file of the process"""

        replace_dict = {}
        
        replace_dict['model_name'] = self.model_name
        replace_dict['process_class'] = self.process_class
        replace_dict['nfinal'] = self.nfinal
        replace_dict['ninitial'] = self.ninitial
        replace_dict['nexternal'] = self.nexternal
        replace_dict['helicity_matrix'] = self.get_helicity_matrix(self.matrix_elements[0])

        replace_dict['all_wavefunction_definitions'] = \
                      """// Wavefunctions
                      void calculate_wavefunctions(const int perm[], const int hel[]);
                      std::complex<double> amp[%d];\n""" % (len(self.amplitudes.get_all_amplitudes()))

        replace_dict['all_matrix_definitions'] = "// Matrix elements\n" + \
                       "\n".join(["double matrix_%s();" % \
                                  me.get('processes')[0].shell_string().\
                                  replace("0_", "") \
                                  for me in self.matrix_elements])

        return self.read_template_file((_template_dir, self.process_class_template)) % replace_dict


    def get_process_function_definitions(self):
        """Template values for the class definition in the source file of the process"""

        replace_dict = {}

        replace_dict['model_name'] = self.model_name
        replace_dict['process_class'] = self.process_class

        color_amplitudes = [ me.get_color_amplitudes() for me in self.matrix_elements ]
        replace_dict['constructor_lines'] = self.get_constructor_lines(self.matrix_elements[0], color_amplitudes)

        replace_dict['nexternal'] = self.nexternal
        replace_dict['finalstates_map'] = self.get_finalstates_map()
        replace_dict['matrix_averaging'] = self.get_matrix_averaging(color_amplitudes)
        replace_dict['matrix_evaluations'] = self.get_matrix_evaluations(color_amplitudes)

        return self.read_template_file((_template_dir, self.process_definition_template)) % replace_dict


    def get_process_name(self):
        """Return process file name for the process in matrix_element"""

        process_string = self.process_string

        # Extract process number
        proc_number_pattern = re.compile("^(.+)@\s*(\d+)\s*(.*)$")
        proc_number_re = proc_number_pattern.match(process_string)
        proc_number = 0
        if proc_number_re:
            proc_number = int(proc_number_re.group(2))
            process_string = proc_number_re.group(1) + \
                             proc_number_re.group(3)

        # Remove order information
        order_pattern = re.compile("^(.+)\s+(\w+)\s*=\s*(\d+)\s*$")
        order_re = order_pattern.match(process_string)
        while order_re:
            process_string = order_re.group(1)
            order_re = order_pattern.match(process_string)
        
        process_string = process_string.replace(' ', '')
        process_string = process_string.replace('>', '_')
        process_string = process_string.replace('+', 'p')
        process_string = process_string.replace('-', 'm')
        process_string = process_string.replace('~', 'x')
        process_string = process_string.replace('/', '_no_')
        process_string = process_string.replace('$', '_nos_')
        process_string = process_string.replace('|', '_or_')
        if proc_number != 0:
            process_string = "%d_%s" % (proc_number, process_string)

        process_string = "Sigma_%s_%s" % (self.model_name,
                                          process_string)
        return process_string


    def get_process_info_lines(self, matrix_element):
        """Return info lines describing the processes for this matrix element"""

        return"\n".join([ "# " + process.nice_string().replace('\n', '\n# * ') \
                         for process in matrix_element.get('processes')])


    def get_constructor_lines(self, matrix_element, color_amplitudes):
        """Get constructor lines for function definition for process source file"""

        constructor_lines = []

        constructor_lines.append("// Set external particle masses for this matrix element")

        for part in matrix_element.get_external_wavefunctions():
            constructor_lines.append("mME.push_back(std::ref(params->%s));" % part.get('mass'))

        return "\n".join(constructor_lines)


    def get_finalstates_map(self):
        """Build map of final states (instanciates SubProcess class)  """
        
        final_states = {}
       
        # First retrieve all the final state's SubProcess definitions

        for me in self.matrix_elements:
            
            proc = me.get('processes')[0]
            
            final_ids = "{" + ",".join( [ str(i) for i in proc.get_final_ids_after_decay() ] ) + "}"
            
            iproc = {}
            
            iproc["function"] = "&%s::matrix_%s" % (self.process_class, proc.shell_string().replace("0_", ""))
            iproc["mirror"] = ( me.get('has_mirror_process') and "true" ) or "false"
            iproc["istates"] = "{" + ",".join( [ "std::make_pair(%i, %i)" % (proc.get('legs')[0].get('id'), proc.get('legs')[1].get('id')) for proc in me.get('processes') ] ) + "}"
            iproc["ncomb"] = me.get_helicity_combinations() 
            iproc["denom"] = me.get_denominator_factor() 
            
            final_states[final_ids] = final_states.get(final_ids, []) + [iproc]

        # Then define the actual final states map using these

        out  = ""
        
        for final, data in final_states.items():
            
            out += "mapFinalStates[%s] =\n" % (final)
            out += "{\n"
            out += ",\n".join( [ "{%(function)s,\n %(mirror)s,\n %(istates)s,\n %(ncomb)i,\n %(denom)i\n}\n" % dati for dati in data ] )
            out += "};\n"

        return out
            

    def get_matrix_averaging(self, color_amplitudes):
        """Get matrix call and averaging loop for process source file"""
        
        replace_dict = {}

        replace_dict['ncomb'] = self.matrix_elements[0].get_helicity_combinations()
        replace_dict['nexternal'] = self.nexternal
        
        return self.read_template_file((_template_dir, self.process_matrix_averaging_template)) % replace_dict


    def get_matrix_evaluations(self, color_amplitudes):
        """Get matrix evaluation functions for process source file""" 

        ret_lines = []
        
        ret_lines.append("void %s::calculate_wavefunctions(const int perm[], const int hel[]) {" % self.process_class)
        ret_lines.append("// Calculate wavefunctions for all processes")
        ret_lines.append(self.get_calculate_wavefunctions(self.wavefunctions, self.amplitudes))
        ret_lines.append("}")
        
        ret_lines.extend( 
                [ self.get_matrix_single_process(i, me, color_amplitudes[i]) for i, me in enumerate(self.matrix_elements) ]
                )
        return "\n".join(ret_lines)


    def get_calculate_wavefunctions(self, wavefunctions, amplitudes):
        """Return the lines for optimized calculation of the wavefunctions for all subprocesses"""

        replace_dict = {}

        replace_dict['nwavefuncs'] = len(wavefunctions)
        
        # Ensure no recycling of wavefunction ! incompatible with some output
        for me in self.matrix_elements:
            me.restore_original_wavefunctions()

        replace_dict['wavefunction_calls'] = "\n".join( self.helas_call_writer.get_wavefunction_calls(helas_objects.HelasWavefunctionList(wavefunctions)) )

        # Change vector format for 4-vectors
        replace_dict['wavefunction_calls'] = re.sub(r"p\[perm\[(\d+)]\]", r'&momenta[perm[\1]][0]', replace_dict['wavefunction_calls'])

        replace_dict['amplitude_calls'] = "\n".join( self.helas_call_writer.get_amplitude_calls(amplitudes) )

        # Change way parameters are called from Parameters_X class
        replace_dict['wavefunction_calls'] = replace_dict['wavefunction_calls'].replace('pars->', 'params->') 
        replace_dict['amplitude_calls'] = replace_dict['amplitude_calls'].replace('pars->', 'params->') 
        
        return self.read_template_file((_template_dir, self.process_wavefunction_template)) % replace_dict


    def get_matrix_single_process(self, i, matrix_element, color_amplitudes):
        """Write matrix() for each process"""

        replace_dict = {}

        # Process name
        replace_dict['proc_name'] = matrix_element.get('processes')[0].shell_string().replace("0_", "")

        # Process class
        replace_dict['process_class'] = self.process_class 
        
        # Process number
        replace_dict['proc_number'] = i

        # Number of color flows
        replace_dict['ncolor'] = len(color_amplitudes)

        # Get color matrix
        replace_dict['color_matrix_lines'] = self.get_color_matrix_lines(matrix_element)
        
        # Get color flow coefficients
        replace_dict['jamp_lines'] = self.get_jamp_lines(color_amplitudes)
        
        return self.read_template_file((_template_dir, self.single_process_template)) % replace_dict


    def get_helicity_matrix(self, matrix_element):
        """Return the Helicity matrix definition lines for this matrix element"""

        helicity_line = "const int helicities[%s][%s] = {" % (self.matrix_elements[0].get_helicity_combinations(), self.nexternal);
        helicity_line_list = []

        for helicities in matrix_element.get_helicity_matrix(allow_reverse=False):
            helicity_line_list.append("{" + ",".join(['%d'] * len(helicities)) % tuple(helicities) + "}")

        return helicity_line + ",".join(helicity_line_list) + "};"

    
    def get_den_factor_line(self, matrix_element):
        """Return the denominator factor line for this matrix element"""

        return "const int denominator = %d;" % matrix_element.get_denominator_factor()


    def get_color_matrix_lines(self, matrix_element):
        """Return the color matrix definition lines for this matrix element. 
        Split rows in chunks of size n."""

        ncolor = str(len(matrix_element.get_color_amplitudes()))

        if not matrix_element.get('color_matrix'):
            
            return "static const double denom[1] = {1.};\nstatic const double cf[1][1] = {1.};"
        
        else:
            
            # First define denominator array
            color_denominators = matrix_element.get('color_matrix').get_line_denominators()
            denom_string = "static const double denom[" + ncolor + "] = {%s};" % \
                           ",".join( [ "%i" % denom for denom in color_denominators ] )

            matrix_strings = []
            my_cs = color.ColorString()
            
            for index, denominator in enumerate(color_denominators):
                # Then write the numerators for the matrix elements
                num_list = matrix_element.get('color_matrix').get_line_numerators(index, denominator)

                matrix_strings.append( "{%s}" % ",".join( [ "%d" % i for i in num_list ] ) )
            
            matrix_string = "static const double cf[" + ncolor + "][" + ncolor + "] = {" + ",".join(matrix_strings) + "};"
            
            return "\n".join([denom_string, matrix_string])


    def get_jamp_lines(self, color_amplitudes):
        """Return the jamp = sum(fermionfactor * amp[i]) lines"""

        declare_ci = False
        res_list = []

        for i, coeff_list in enumerate(color_amplitudes):

            res = "jamp[%i]=" % i

            # Optimization: if all contributions to that color basis element have
            # the same coefficient (up to a sign), put it in front
            list_fracs = [ abs(coefficient[0][1]) for coefficient in coeff_list ]
            common_factor = False
            diff_fracs = list(set(list_fracs))
            if len(diff_fracs) == 1 and abs(diff_fracs[0]) != 1:
                common_factor = True
                global_factor = diff_fracs[0]
                res = res + '%s(' % coeff(1, global_factor, False, 0)

            for (coefficient, amp_number) in coeff_list:
                if not declare_ci:
                    declare_ci = coefficient[2]
                if common_factor:
                    res = res + "%samp[%d]" % (coeff(coefficient[0],
                                               coefficient[1] / abs(coefficient[1]),
                                               coefficient[2],
                                               coefficient[3]),
                                               amp_number - 1)
                else:
                    res = res + "%samp[%d]" % (coeff(coefficient[0],
                                               coefficient[1],
                                               coefficient[2],
                                               coefficient[3]),
                                               amp_number - 1)

            if common_factor:
                res = res + ')'

            res += ';'

            res_list.append(res)

        if declare_ci:
            res_list.insert(0, "static const std::complex<double> cI(0., 1.);")

        return "\n".join(res_list)


    @staticmethod
    def get_model_name(name):
        """Replace - with _, + with _plus_ in a model name."""

        name = name.replace('-', '_')
        name = name.replace('+', '_plus_')
        return name



#===============================================================================
# Global helper methods
#===============================================================================

def expand_initial_state(i_id):
    if( -4 <= i_id <= 4  and i_id != 0):
        return [i_id < 0 and -i or i  for i in [1,2,3,4]]
    return [i_id]

def coeff(ff_number, frac, is_imaginary, Nc_power, Nc_value=3):
    """Returns a nicely formatted string for the coefficients in JAMP lines"""

    total_coeff = ff_number * frac * fractions.Fraction(Nc_value) ** Nc_power

    if total_coeff == 1:
        if is_imaginary:
            return '+cI*'
        else:
            return '+'
    elif total_coeff == -1:
        if is_imaginary:
            return '-cI*'
        else:
            return '-'

    res_str = '%+i.' % total_coeff.numerator

    if total_coeff.denominator != 1:
        # Check if total_coeff is an integer
        res_str = res_str + '/%i.' % total_coeff.denominator

    if is_imaginary:
        res_str = res_str + '*cI'

    return res_str + '*'


class ProcessExporterMoMEMta(VirtualExporter):
    """Plugin class handling the export of processes for MoMEMta in C++"""

    # Check status of the directory (ask to remove it if already exists)
    check = True 
    # Language type: 'v4' for f77/ 'cpp' for C++ output
    exporter = 'cpp'
    # Output type:
    #[Template/dir/None] copy the Template, just create dir or do nothing 
    output = 'Template'
    # Decide which type of merging if used [madevent/madweight]
    grouped_mode = 'madweight'
    # If no grouping on can decide to merge uu~ and u~u anyway:
    sa_symmetry = True
    
    def __init__(self, dir_path="", opt=None):
       
        # Output directory
        self.dir_path = dir_path

        # Not used by us
        self.opt = dict()
        if opt:
            self.opt.update(opt)

        # The class actually handling the exporting of the ME
        self.Exporter = OneProcessExporterMoMEMta

        # Will contain all sub-process directories created hereafter
        self.sub_dirs = []

    #===============================================================================
    # copy_template
    #===============================================================================
    def copy_template(self, model):
        """Prepare dir_path as output directory, including:
        include (for model and ALOHA header files)
        src (for model and ALOHA source files)
        lib (with compiled libraries from src)
        SubProcesses (with makefile and Pxxxxx directories)
        """
   
        self.model = model
        self.model_name = self.Exporter.get_model_name(self.model.get('name'))

        self.dir_name = os.path.basename(os.path.normpath(self.dir_path))
        # Check that the name used will compile in c++
        name_check = re.compile('\\W') # match any non alphanumeric (excluding '_') character
        if self.dir_name[0].isdigit() or name_check.search(self.dir_name):
            raise Exception('Exported directory name is used as C++ namespace for the process and must therefore be a legal C++ variable name.')
    
        cwd = os.getcwd()
    
        try:
            os.mkdir(self.dir_path)
        except os.error as error:
            logger.warning(error.strerror + " " + self.dir_path)
        
        try:
            os.chdir(self.dir_path)
        except os.error:
            logger.error('Could not cd to directory %s' % self.dir_path)
            return 0
    
        logger.info('Creating subdirectories in directory %s' % self.dir_path)
    
        try:
            os.mkdir('include')
        except os.error as error:
            logger.warning(error.strerror + " " + self.dir_path)
        
        try:
            os.mkdir('src')
        except os.error as error:
            logger.warning(error.strerror + " " + self.dir_path)
        
        try:
            os.mkdir('lib')
        except os.error as error:
            logger.warning(error.strerror + " " + self.dir_path)
        
        try:
            os.mkdir('Cards')
        except os.error as error:
            logger.warning(error.strerror + " " + self.dir_path)
        
        try:
            os.mkdir('SubProcesses')
        except os.error as error:
            logger.warning(error.strerror + " " + self.dir_path)
    
        # Write param_card
        with open(os.path.join("Cards","param_card.dat"), 'w') as m_file:
            m_file.write(model.write_param_card())
    
        # Copy the SubProcess base class file into 'include' directory
        subprocess = self.Exporter.read_template_file((_template_dir, 'SubProcess.h')) % \
                               {'namespace': self.dir_name + "_" + self.model_name }
        with open(os.path.join('include', 'SubProcess.h'), 'w') as m_file:
            m_file.write(subprocess)
    
        # Return to original PWD
        os.chdir(cwd)
        self.opt = dict()
    
    

    #===============================================================================
    # generate_subprocess_directory
    #===============================================================================
    def generate_subprocess_directory(self, subproc_group, helicity_model, proc_number=None):
        """Generate the Pxxxxx directory for a subprocess in C++ standalone,
        including the necessary .h and .cc files"""
    
        cwd = os.getcwd()
        
        # Create the process_exporter
        process_exporter = self.Exporter(subproc_group, helicity_model, self.dir_path, self.dir_name)
    
        # Create the directory PN_xx_xxxxx in the specified path
        sub_dir_path = process_exporter.path
        try:
            os.mkdir(sub_dir_path)
        except os.error as error:
            logger.warning(error.strerror + " " + sub_dir_path)
    
        try:
            os.chdir(sub_dir_path)
        except os.error:
            logger.error('Could not cd to directory %s' % sub_dir_path)
            return 0
    
        logger.info('Creating files in directory %s' % sub_dir_path)
    
        # Create the process .h and .cc files
        process_exporter.generate_process_files()

        # Log created dir
        self.sub_dirs.append(sub_dir_path)
    
        # Return to original PWD
        os.chdir(cwd)

        return 0
    
    #===============================================================================
    # Routines to export/output UFO models in C++ format
    #===============================================================================
    
    def convert_model(self, model, wanted_lorentz = [], wanted_couplings = []):
        """Create a full valid C++ model from an MG5 model (coming from UFO)"""
    
        # Create the files for model parameter and amplitude calls
        model_builder = UFOModelConverterMoMEMta(
                                            self.dir_name,
                                            self.model,
                                            self.dir_path,
                                            wanted_lorentz,
                                            wanted_couplings)
    
        model_builder.write_files()

    def finalize(self, matrix_element, cmdhistory, MG5options, outputflag):
        # Copy CMakeLists.txt and fill template
        include_commands = "\n".join( [ 'include_directories("SubProcesses/{}")'.format(os.path.basename(dir)) for dir in self.sub_dirs ] )
        makefile = OneProcessExporterMoMEMta.read_template_file((_template_dir, 'CMakeLists.txt')) % {
                'dir_name': self.dir_name,
                'include_subprocs_list': include_commands
                }
        with open(os.path.join(self.dir_path, 'CMakeLists.txt'), 'w') as m_file:
            m_file.write(makefile)

    def modify_grouping(self, matrix_element):
        return False, matrix_element


#===============================================================================
# UFOModelConverterMoMEMta
#===============================================================================

class UFOModelConverterMoMEMta(UFOModelConverterCPP):
    """Subclass of UFOModelConverterCPP, since we modify certain things, such as the parameters class"""

    include_dir = "include"
    cc_file_dir = "src"
    param_template_h = (_template_dir, 'model_parameters.h')
    param_template_cc = (_template_dir, 'model_parameters.cc')

    def __init__(self, namespace, *args, **kwargs):
        self.namespace = namespace
        return UFOModelConverterCPP.__init__(self, *args, **kwargs)
    
    def generate_parameters_class_files(self):
        """Create the content of the Parameters_model.h and .cc files"""

        replace_dict = {}

        replace_dict['model_name'] = self.model_name
        replace_dict['namespace'] = self.namespace

        # Unchanged compared to base class
        replace_dict['independent_parameters'] = \
                                   "// Model parameters independent of aS\n" + \
                                   self.write_parameters(self.params_indep)
        replace_dict['independent_couplings'] = \
                                   "// Model parameters dependent on aS\n" + \
                                   self.write_parameters(self.params_dep)
        replace_dict['dependent_parameters'] = \
                                   "// Model couplings independent of aS\n" + \
                                   self.write_parameters(self.coups_indep)
        replace_dict['dependent_couplings'] = \
                                   "// Model couplings dependent on aS\n" + \
                                   self.write_parameters(self.coups_dep.values())

        replace_dict['set_independent_couplings'] = \
                               self.write_set_parameters(self.coups_indep)
        replace_dict['set_dependent_parameters'] = \
                               self.write_set_parameters(self.params_dep)
        replace_dict['set_dependent_couplings'] = \
                               self.write_set_parameters(self.coups_dep.values())

        # This part is modified by us
        
        # First retrieve list of params read from the card, or not:
        params_indep_card = []
        params_indep_nocard = []
        for param in self.params_indep:
            if 'slha' in param.expr:
                params_indep_card.append(param)
            else:
                params_indep_nocard.append(param)
       
        # This goes in the constructor: initialise map of parameters
        replace_dict['parameter_map_lines'] = self.write_set_parameters(params_indep_card)
        replace_dict['parameter_map_lines'] = re.sub(r'(.*) = slha', r'm_card_parameters["\1"] = card', replace_dict['parameter_map_lines'])
        
        # In the method: retrieve parameters from map, or from expression for other parameters
        replace_dict['set_independent_parameters'] = \
                               self.write_parameters_from_map(params_indep_card)
        replace_dict['set_independent_parameters'] += \
                               self.write_set_parameters(params_indep_nocard)
        
        file_h = self.read_template_file(self.param_template_h) % replace_dict
        file_cc = self.read_template_file(self.param_template_cc) % replace_dict
        
        return file_h, file_cc

    
    def write_parameters_from_map(self, params):
        """Write out the lines of independent parameters"""

        res_strings = []
        for param in params:
            res_strings.append('%(param)s = m_card_parameters["%(param)s"];' % {'param': param.name})

        return "\n".join(res_strings)
