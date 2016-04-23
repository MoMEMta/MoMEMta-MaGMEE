################################################################################
#
# Copyright (c) 2009 The MadGraph5_aMC@NLO Development team and Contributors
#
# This file is a part of the MadGraph5_aMC@NLO project, an application which 
# automatically generates Feynman diagrams and matrix elements for arbitrary
# high-energy processes in the Standard Model and beyond.
#
# It is subject to the MadGraph5_aMC@NLO license which should accompany this 
# distribution.
#
# For more information, visit madgraph.phys.ucl.ac.be and amcatnlo.web.cern.ch
#
################################################################################

"""Methods and classes to export models and matrix elements to MEM-dedicated C++ Standalone format."""

import fractions
import glob
import itertools
import logging
from math import fmod
import os
import re
import shutil
import subprocess

import madgraph.core.base_objects as base_objects
import madgraph.core.color_algebra as color
import madgraph.core.helas_objects as helas_objects
import madgraph.iolibs.drawing_eps as draw
import madgraph.iolibs.files as files
import madgraph.iolibs.helas_call_writers as helas_call_writers
import madgraph.iolibs.file_writers as writers
import madgraph.iolibs.template_files as template_files
import madgraph.iolibs.ufo_expression_parsers as parsers
import madgraph.iolibs.group_subprocs as group_subprocs
from madgraph import MadGraph5Error, InvalidCmd, MG5DIR
from madgraph.iolibs.files import cp, ln, mv
from madgraph.iolibs.export_cpp import UFOModelConverterCPP

import madgraph.various.misc as misc

import aloha.create_aloha as create_aloha
import aloha.aloha_writers as aloha_writers

_file_path = os.path.join(MG5DIR, "PLUGIN", "MoMEMta-MaGMEE")
_template_dir = os.path.join(_file_path, "Template")
logger = logging.getLogger('madgraph.export_pythia8')


#===============================================================================
# ProcessExporterCPP
#===============================================================================
class OneProcessExporterMoMEMta(object):
    """Class to take care of exporting a set of matrix elements to
    C++ format."""

    # Static variables (for inheritance)
    process_dir = '.'
    include_dir = '.'
    process_template_h                 = 'process_h.inc'
    process_template_cc                = 'process_cc.inc'
    process_class_template             = 'process_class.inc'
    process_definition_template        = 'function_definitions.inc'
    process_wavefunction_template      = 'wavefunctions.inc'
    process_sigmaKin_function_template = 'me_function.inc'
    single_process_template            = 'matrix.inc'

    class OneProcessExporterMoMEMtaError(Exception):
        pass
    
    def __init__(self, matrix_elements, cpp_helas_call_writer, process_string = "",
                 process_number = 0, path = os.getcwd()):
        """Initiate with matrix elements, helas call writer, process
        string, path. Generate the process .h and .cc files."""

        self.ufolder = ''

        if isinstance(matrix_elements, helas_objects.HelasMultiProcess) or isinstance(matrix_elements, group_subprocs.SubProcessGroup):
            self.matrix_elements = matrix_elements.get('matrix_elements')
        elif isinstance(matrix_elements, helas_objects.HelasMatrixElement):
            self.matrix_elements = \
                         helas_objects.HelasMatrixElementList([matrix_elements])
        elif isinstance(matrix_elements, helas_objects.HelasMatrixElementList):
            self.matrix_elements = matrix_elements
        else:
            raise base_objects.PhysicsObject.PhysicsObjectError,\
                  "Wrong object type for matrix_elements"

        if not self.matrix_elements:
            raise MadGraph5Error("No matrix elements to export")

        self.model = self.matrix_elements[0].get('processes')[0].get('model')
        self.model_name = OneProcessExporterMoMEMta.get_model_name(self.model.get('name'))

        self.processes = sum([me.get('processes') for \
                              me in self.matrix_elements], [])
        self.processes.extend(sum([me.get_mirror_processes() for \
                              me in self.matrix_elements], []))

        self.nprocesses = len(self.matrix_elements)
        self.nprocesses += sum([1 for me in self.matrix_elements if me.get('has_mirror_process')])

        if process_string:
            self.process_string = process_string
        else:
            self.process_string = self.processes[0].base_string()

        if process_number:
            self.process_number = process_number
        else:
            self.process_number = self.processes[0].get('id')

        self.process_name = self.get_process_name()
        self.process_class = "CPPProcess"

        self.path = path
        self.helas_call_writer = cpp_helas_call_writer

        if not isinstance(self.helas_call_writer, helas_call_writers.CPPUFOHelasCallWriter):
            raise self.OneProcessExporterMoMEMtaError, \
                "helas_call_writer not CPPUFOHelasCallWriter"

        self.nexternal, self.ninitial = \
                        self.matrix_elements[0].get_nexternal_ninitial()
        self.nfinal = self.nexternal - self.ninitial

        # Check if we can use the same helicities for all matrix
        # elements
        
        hel_matrix = self.get_helicity_matrix(self.matrix_elements[0])

        # Test if all subprocesses have the same helicity structure, abort export otherwise
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
        """Generate the .h and .cc files needed for C++, for the
        processes described by multi_matrix_element"""

        # Create the files
        if not os.path.isdir(os.path.join(self.path, self.include_dir)):
            os.makedirs(os.path.join(self.path, self.include_dir))
        filename = os.path.join(self.path, self.include_dir,
                                '%s.h' % self.process_class)
        
        self.write_process_h_file(writers.CPPWriter(filename))

        if not os.path.isdir(os.path.join(self.path, self.process_dir)):
            os.makedirs(os.path.join(self.path, self.process_dir))
        filename = os.path.join(self.path, self.process_dir,
                                '%s.cc' % self.process_class)

        self.write_process_cc_file(writers.CPPWriter(filename))

        logger.info('Created files %(process)s.h and %(process)s.cc in' % \
                    {'process': self.process_class} + \
                    ' directory %(dir)s' % {'dir': os.path.split(filename)[0]})



    #===========================================================================
    # write_process_h_file
    #===========================================================================
    def write_process_h_file(self, writer):
        """Write the class definition (.h) file for the process"""
        
        if not isinstance(writer, writers.CPPWriter):
            raise writers.CPPWriter.CPPWriterError(\
                "writer not CPPWriter")

        replace_dict = {}

        # Extract version number and date from VERSION file
        info_lines = get_mg5_info_lines()
        replace_dict['info_lines'] = info_lines

        # Extract model name
        replace_dict['model_name'] = \
                         self.model_name

        # Extract process file name
        replace_dict['process_file_name'] = self.process_name

        # Extract class definitions
        process_class_definitions = self.get_process_class_definitions()
        replace_dict['process_class_definitions'] = process_class_definitions

        file = self.read_template_file((_template_dir, self.process_template_h)) % replace_dict

        # Write the file
        writer.writelines(file)

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

        # Extract version number and date from VERSION file
        info_lines = get_mg5_info_lines()
        replace_dict['info_lines'] = info_lines

        # Extract process file name
        replace_dict['process_file_name'] = self.process_name

        replace_dict['process_class'] = self.process_class

        # Extract model name
        replace_dict['model_name'] = self.model_name
                         

        # Extract class function definitions
        process_function_definitions = \
                         self.get_process_function_definitions()
        replace_dict['process_function_definitions'] = \
                                                   process_function_definitions

        file = self.read_template_file((_template_dir, self.process_template_cc)) % replace_dict

        # Write the file
        writer.writelines(file)

    #===========================================================================
    # Process export helper functions
    #===========================================================================
    def get_process_class_definitions(self):
        """The complete class definition for the process"""

        replace_dict = {}

        # Extract model name
        replace_dict['model_name'] = self.model_name

        # Extract process info lines for all processes
        process_lines = "\n".join([self.get_process_info_lines(me) for me in \
                                   self.matrix_elements])
        
        replace_dict['process_lines'] = process_lines

        # Extract number of external particles
        replace_dict['nfinal'] = self.nfinal

        # Extract number of external particles
        replace_dict['ninitial'] = self.ninitial

        replace_dict['nexternal'] = self.nexternal

        # Extract process class name 
        replace_dict['process_class'] = self.process_class

        # Extract process definition
        process_definition = "%s (%s)" % (self.process_string,
                                          self.model_name)
        replace_dict['process_definition'] = process_definition

        process = self.processes[0]

        # Extract helicity matrix
        replace_dict['helicity_matrix'] = self.get_helicity_matrix(self.matrix_elements[0])

        replace_dict['all_sigma_kin_definitions'] = \
                      """// Calculate wavefunctions
                      void calculate_wavefunctions(const int perm[], const int hel[]);
                      std::complex<double> amp[%d];""" % (len(self.amplitudes.get_all_amplitudes()))

        replace_dict['all_matrix_definitions'] = \
                       "\n".join(["double matrix_%s();" % \
                                  me.get('processes')[0].shell_string().\
                                  replace("0_", "") \
                                  for me in self.matrix_elements])

        file = self.read_template_file((_template_dir, self.process_class_template)) % replace_dict

        return file

    def get_process_function_definitions(self):
        """The complete Pythia 8 class definition for the process"""

        replace_dict = {}

        # Extract model name
        replace_dict['model_name'] = self.model_name

        # Extract process info lines
        replace_dict['process_lines'] = \
                             "\n".join([self.get_process_info_lines(me) for \
                                        me in self.matrix_elements])

        # Extract process class name 
        replace_dict['process_class'] = self.process_class

        color_amplitudes = [me.get_color_amplitudes() for me in \
                            self.matrix_elements]

        replace_dict['nexternal'] = self.nexternal

        replace_dict['initProc_lines'] = \
                                self.get_initProc_lines(self.matrix_elements[0],
                                                        color_amplitudes)

        replace_dict['finalstates_map'] = self.get_finalstates_map()

        replace_dict['sigmaKin_lines'] = \
                                     self.get_sigmaKin_lines(color_amplitudes)

        replace_dict['all_sigmaKin'] = \
                                  self.get_all_sigmaKin_lines(color_amplitudes,
                                                              self.process_class)

        file = self.read_template_file((_template_dir, self.process_definition_template)) %\
               replace_dict

        return file

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


    def get_initProc_lines(self, matrix_element, color_amplitudes):
        """Get initProc_lines for function definition for Pythia 8 .cc file"""

        initProc_lines = []

        initProc_lines.append("// Set external particle masses for this matrix element")

        for part in matrix_element.get_external_wavefunctions():
            initProc_lines.append("mME.push_back(params.%s);" % part.get('mass'))

        return "\n".join(initProc_lines)

    def get_calculate_wavefunctions(self, wavefunctions, amplitudes):
        """Return the lines for optimized calculation of the
        wavefunctions for all subprocesses"""

        replace_dict = {}

        replace_dict['nwavefuncs'] = len(wavefunctions)
        
        #ensure no recycling of wavefunction ! incompatible with some output
        for me in self.matrix_elements:
            me.restore_original_wavefunctions()

        replace_dict['wavefunction_calls'] = "\n".join(\
            self.helas_call_writer.get_wavefunction_calls(\
            helas_objects.HelasWavefunctionList(wavefunctions)))

        # change vector format for 4-vectors
        replace_dict['wavefunction_calls'] = re.sub(r"p\[perm\[(\d+)]\]",r'&momenta[perm[\1]][0]',replace_dict['wavefunction_calls'])

        replace_dict['amplitude_calls'] = "\n".join(\
            self.helas_call_writer.get_amplitude_calls(amplitudes))

        # Change way parameters are called from Parameters_X class
        replace_dict['wavefunction_calls'] = replace_dict['wavefunction_calls'].replace('pars->', 'params.') 
        replace_dict['amplitude_calls'] = replace_dict['amplitude_calls'].replace('pars->', 'params.') 
        
        file = self.read_template_file((_template_dir, self.process_wavefunction_template)) % \
                replace_dict

        return file

    def get_finalstates_map(self):
        """Build map of final states with __matrixElements  """
        final_states = {}
        for me in self.matrix_elements:
            proc = me.get('processes')[0]
            #f_s = "{"+",".join([str(part.get('id')) for part in proc.get('legs')[2:]])+"}"
            f_s = "{"+",".join([str(i) for i in proc.get_final_ids_after_decay()])+"}"
            iproc = {}
            iproc["function"] = "&%s::matrix_%s" % (self.process_class, proc.shell_string().replace("0_", ""))
            iproc["mirror"] = me.get('has_mirror_process') and "true" or "false"
            iproc["istates"] = "{"+",".join(["std::make_pair(%i,%i)" % (proc.get('legs')[0].get('id'), proc.get('legs')[1].get('id')) for proc in me.get('processes') ])+"}"
            iproc["ncomb"] = me.get_helicity_combinations() 
            iproc["denom"] = me.get_denominator_factor() 
            final_states[f_s] = final_states.get(f_s,[])+[iproc]

        out  = ""
        for final, data in final_states.items():
            out += "mapFinalStates[%s] =\n" % (final)
            out += "{\n"
            out += ",\n".join(["{ %(function)s,\n %(mirror)s,\n %(istates)s ,\n %(ncomb)i,\n %(denom)i\n }\n" % dati for dati in data])
            out += "};\n"

        return out
            

    def get_sigmaKin_lines(self, color_amplitudes):
        """Get sigmaKin_lines for function definition for .cc file"""

        
        replace_dict = {}

        # Number of helicity combinations
        replace_dict['ncomb'] = self.matrix_elements[0].get_helicity_combinations()
        replace_dict['nexternal'] = self.nexternal

        # Process name
        replace_dict['process_class_name'] = self.process_name
        
        file = self.read_template_file((_template_dir, self.process_sigmaKin_function_template)) % replace_dict

        return file

    def get_all_sigmaKin_lines(self, color_amplitudes, class_name):
        """Get sigmaKin_process for all subprocesses for .cc file"""

        ret_lines = []
        
        ret_lines.append(\
            "void %s::calculate_wavefunctions(const int perm[], const int hel[]){" % \
            class_name)
        ret_lines.append("// Calculate wavefunctions for all processes")
        ret_lines.append(self.get_calculate_wavefunctions(\
            self.wavefunctions, self.amplitudes))
        ret_lines.append("}")
        
        ret_lines.extend([self.get_matrix_single_process(i, me,
                                                         color_amplitudes[i],
                                                         class_name) \
                                for i, me in enumerate(self.matrix_elements)])
        return "\n".join(ret_lines)


    def get_matrix_single_process(self, i, matrix_element, color_amplitudes,
                                  class_name):
        """Write matrix() for each process"""

        # Write matrix() for the process

        replace_dict = {}

        # Process name
        replace_dict['proc_name'] = \
          matrix_element.get('processes')[0].shell_string().replace("0_", "")
        

        # Wavefunction and amplitude calls
        replace_dict['matrix_args'] = ""

        # Process name
        replace_dict['process_class'] = class_name
        
        # Process number
        replace_dict['proc_number'] = i

        # Number of color flows
        replace_dict['ncolor'] = len(color_amplitudes)

        # Extract color matrix
        replace_dict['color_matrix_lines'] = \
                                     self.get_color_matrix_lines(matrix_element)

                                     
        replace_dict['jamp_lines'] = self.get_jamp_lines(color_amplitudes)


        #specific exporter hack
        replace_dict =  self.get_class_specific_definition_matrix(replace_dict, matrix_element)
        
        file = self.read_template_file((_template_dir, self.single_process_template)) % \
                replace_dict

        return file

    def get_class_specific_definition_matrix(self, converter, matrix_element):
        """place to add some specific hack to a given exporter.
        Please always use Super in that case"""

        return converter

    def get_helicity_matrix(self, matrix_element):
        """Return the Helicity matrix definition lines for this matrix element"""

        helicity_line = "const int helicities[%s][%s] = {" % (self.matrix_elements[0].get_helicity_combinations(), self.nexternal);
        helicity_line_list = []

        for helicities in matrix_element.get_helicity_matrix(allow_reverse=False):
            helicity_line_list.append("{"+",".join(['%d'] * len(helicities)) % \
                                       tuple(helicities) + "}")

        return helicity_line + ",".join(helicity_line_list) + "};"

    def get_den_factor_line(self, matrix_element):
        """Return the denominator factor line for this matrix element"""

        return "const int denominator = %d;" % \
               matrix_element.get_denominator_factor()

    def get_color_matrix_lines(self, matrix_element):
        """Return the color matrix definition lines for this matrix element. Split
        rows in chunks of size n."""

        ncolor = str(len(matrix_element.get_color_amplitudes()))

        if not matrix_element.get('color_matrix'):
            return "\n".join(["static const double denom[1] = {1.};",
                              "static const double cf[1][1] = {1.};"])
        else:
            color_denominators = matrix_element.get('color_matrix').\
                                                 get_line_denominators()
            denom_string = "static const double denom["+ncolor+"] = {%s};" % \
                           ",".join(["%i" % denom for denom in color_denominators])

            matrix_strings = []
            my_cs = color.ColorString()
            for index, denominator in enumerate(color_denominators):
                # Then write the numerators for the matrix elements
                num_list = matrix_element.get('color_matrix').\
                                            get_line_numerators(index, denominator)

                matrix_strings.append("{%s}" % \
                                     ",".join(["%d" % i for i in num_list]))
            matrix_string = "static const double cf["+ncolor+"]["+ncolor+"] = {" + \
                            ",".join(matrix_strings) + "};"
            return "\n".join([denom_string, matrix_string])

    def get_jamp_lines(self, color_amplitudes):
        """Return the jamp = sum(fermionfactor * amp[i]) lines"""

        declare_ci = False
        res_list = []

        for i, coeff_list in enumerate(color_amplitudes):

            res = "jamp[%i]=" % i

            # Optimization: if all contributions to that color basis element have
            # the same coefficient (up to a sign), put it in front
            list_fracs = [abs(coefficient[0][1]) for coefficient in coeff_list]
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

#def read_template_file(filename):
#    """Open a template file and return the contents."""
#
#    return open(os.path.join(_template_dir, filename)).read()

def get_mg5_info_lines():
    """Return info lines for MG5, suitable to place at beginning of
    Fortran files"""

    info = misc.get_pkg_info()
    info_lines = ""
    if info and info.has_key('version') and  info.has_key('date'):
        info_lines = "#  MadGraph5_aMC@NLO v. %s, %s\n" % \
                     (info['version'], info['date'])
        info_lines = info_lines + \
                     "#  By the MadGraph5_aMC@NLO Development Team\n" + \
                     "#  Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch"
    else:
        info_lines = "#  MadGraph5_aMC@NLO\n" + \
                     "#  By the MadGraph5_aMC@NLO Development Team\n" + \
                     "#  Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch"        

    return info_lines

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


class ProcessExporterMoMEMta(object):

    default_opt = {'clean': False, 'complex_mass':False,
            'export_format': 'madevent', 'mp': False,
            'v5_model': True
            }
    grouped_mode = True 

    def __init__(self, dir_path="", opt=None):
        self.dir_path = dir_path

        self.opt = dict(self.default_opt)
        if opt:
            self.opt.update(opt)

    #===============================================================================
    # setup_cpp_standalone_dir
    #===============================================================================
    def setup_cpp_standalone_dir(self, model):
        """Prepare export_dir as standalone_cpp directory, including:
        src (for model and ALOHA files + makefile)
        lib (with compiled libraries from src)
        SubProcesses (with makefile and Pxxxxx directories)
        """
   
        self.model = model

        process_class = os.path.basename(os.path.normpath(self.dir_path))
        # Check that the name used will compile in c++
        name_check = re.compile('\\W') # match any non alphanumeric (excluding '_') character
        if process_class[0].isdigit() or name_check.search(process_class):
            raise Exception('Exported directory name is used as C++ class name for the process and must therefore be a legal C++ variable name.')
    
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
        open(os.path.join("Cards","param_card.dat"), 'w').write(\
            model.write_param_card())
    
        src_files = ['read_slha.h', 'read_slha.cc']
        
        # Copy the needed src files
        for f in src_files:
            cp(_template_dir + '/' + f, 'src')
    
        # Copy the base class file into 'src' directory
        cp(_template_dir + '/baseclasses.h', 'src/process_base_classes.h')
    
        # Copy src Makefile
        makefile = OneProcessExporterMoMEMta.read_template_file((_template_dir, 'Makefile_model.inc')) % \
                               {'model': OneProcessExporterMoMEMta.get_model_name(model.get('name'))}
        open(os.path.join('src', 'Makefile'), 'w').write(makefile)
    
        # Copy SubProcesses Makefile
        makefile = OneProcessExporterMoMEMta.read_template_file((_template_dir, 'Makefile_process.inc')) % \
                                        {'model': OneProcessExporterMoMEMta.get_model_name(model.get('name')),
                                         'process_class': process_class
                                        }
        open(os.path.join('SubProcesses', 'Makefile'), 'w').write(makefile)
    
        # Copy CMakeLists.txt
        makefile = OneProcessExporterMoMEMta.read_template_file((_template_dir, 'CMakeLists.inc')) % {'process_class': process_class}
        open('CMakeLists.txt', 'w').write(makefile)
    
        # Return to original PWD
        os.chdir(cwd)
        self.opt = dict()
    
    

    #===============================================================================
    # generate_subprocess_directory_standalone_cpp
    #===============================================================================
    def generate_subprocess_directory(self, matrix_element, cpp_helas_call_writer, proc_number=None):
    
        """Generate the Pxxxxx directory for a subprocess in C++ standalone,
        including the necessary .h and .cc files"""
    
        cwd = os.getcwd()
        # Create the process_exporter
        process_exporter_cpp = OneProcessExporterMoMEMta(matrix_element, cpp_helas_call_writer)
    
        # extract user defined folder name (dirty way)
        #f_search = re.search(".*/(.*)/SubProcess",path)
        #ufolder = f_search and f_search.groups()[0] or "user_folder_not_found"
        #process_exporter_cpp.ufolder = ufolder
        #process_exporter_cpp.process_class = ufolder

        process_exporter_cpp.ufolder = self.dir_path
        process_exporter_cpp.process_class = os.path.basename(self.dir_path)
    
        # Create the directory PN_xx_xxxxx in the specified path
        sub_dirpath = os.path.join(self.dir_path, "SubProcesses",
                       "P%d_%s" % (process_exporter_cpp.process_number,
                                   process_exporter_cpp.process_name))
        try:
            os.mkdir(sub_dirpath)
        except os.error as error:
            logger.warning(error.strerror + " " + sub_dirpath)
    
        try:
            os.chdir(sub_dirpath)
        except os.error:
            logger.error('Could not cd to directory %s' % sub_dirpath)
            return 0
    
        logger.info('Creating files in directory %s' % sub_dirpath)
    
        process_exporter_cpp.path = sub_dirpath
        # Create the process .h and .cc files
        process_exporter_cpp.generate_process_files()
    
        linkfiles = ['Makefile']
        for file in linkfiles:
            ln('../%s' % file)
    
        # Return to original PWD
        os.chdir(cwd)

        return 0
    
    #def make_model_cpp(self, dir_path):
    #    """Make the model library in a C++ standalone directory"""
    #
    #    source_dir = os.path.join(dir_path, "src")
    #    # Run standalone
    #    logger.info("Running make for src")
    #    misc.compile(cwd=source_dir)
    
    #===============================================================================
    # Routines to export/output UFO models in C++ format
    #===============================================================================
    
    def convert_model(self, model, wanted_lorentz = [],
                             wanted_couplings = []):
        """Create a full valid C++ model from an MG5 model (coming from UFO)"""
    
        # create the model parameter files
        model_builder = UFOModelConverterCPP(self.model,
                                             os.path.join(self.dir_path, 'src'),
                                             wanted_lorentz,
                                             wanted_couplings)
    
        # We just want to change the Parameters_model class templates
        model_builder.param_template_h = (_template_dir, 'model_parameters_h.inc')
        model_builder.param_template_cc = (_template_dir, 'model_parameters_cc.inc')
    
        model_builder.write_files()

    def finalize(self, *args, **kwargs):
        pass


ProcessExporterMoMEMta_cfg = {
    # check status of the directory. Remove it if already exists
    'check': True, 
    # Language type: 'v4' for f77/ 'cpp' for C++ output
    'exporter': 'cpp',
    # Output type:
    #[Template/dir/None] copy the Template, just create dir  or do nothing 
    'output': 'Template', 
    # Grouping SubProcesses into quark/lepton (even if not identical matrix element)
    'group_subprocesses': True,
    # Decide which type of merging if used [madevent/madweight]
    'group_mode': 'madweight', 
    # if no grouping on can decide to merge uu~ and u~u anyway:
    'sa_symmetry': True, 
    # The most important part where the exporter is defined:
    # a plugin typically defines this file in another file (here tak one of MG5)
    'exporter_class': ProcessExporterMoMEMta
    }


