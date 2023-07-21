# -*- coding: utf-8 -*-
"""Plugin for the pynball calculation."""
import json

from pathlib import Path

from aiida import orm
from aiida.common import datastructures
from aiida.engine import CalcJob

from aiida_pseudo.data.pseudo.upf import UpfData
from aiida_quantumespresso.calculations.pw import PwCalculation

class PynballCalculation(CalcJob):
    """``CalcJob`` implementation for the Python pynball module."""

    DEFAULT_OUTPUT_FILE = 'speriamobene.txt'
    DEFAULT_INPUT_FILE = 'pinball.json'

    @classmethod
    def define(cls, spec):
        """Define the process."""
        super().define(spec)
        spec.input('metadata.options.parser_name', valid_type=str, 
                   default='hydrogen_restorer.pynball')
        spec.input('parent_folder', valid_type=orm.RemoteData,
                   help='the folder of a completed SCF `PwCalculation`')
        spec.input('all_peaks', valid_type=orm.ArrayData,
                   help='All the possible hydrogen positions.')
        spec.input('number_hydrogen', valid_type=orm.Int,
                   help='Number of hydrogen atoms to place.')
        spec.input('hydrogen_pseudo', valid_type=UpfData,
                   help='The pseudopotential to use for hydrogen.')
        spec.output('final_structure', valid_type=orm.StructureData)
        spec.output('output_parameters', valid_type=orm.Dict)
        spec.exit_code(321, 'ERROR_READING_CIF',
                       'Failed to parse the resulting CIF file.')
        spec.exit_code(322, 'ERROR_READING_JSON',
                    'Failed to parse the resulting JSON file.')
     
    def prepare_for_submission(self, folder):
        """Prepare the calculation."""

        # Prepare the contents of the JSON input file
        peak_positions = self.inputs.all_peaks.get_array('peak_positions')        
        peak_positions = peak_positions.tolist()

        pinball_input = {
            'pwin': PwCalculation._DEFAULT_INPUT_FILE,
            "tot_pinballs": self.inputs.number_hydrogen.value, 
            "pseudo": self.inputs.hydrogen_pseudo.filename, 
            "positions": peak_positions
        }

        with folder.open(self.DEFAULT_INPUT_FILE, 'w') as handle:
            handle.write(json.dumps(pinball_input))


        remote_copy_list = []
        local_copy_list = []

        scf_path = self.inputs.parent_folder.get_remote_path()

        remote_copy_list.append((
            self.inputs.parent_folder.computer.uuid,
            Path(scf_path, PwCalculation._DEFAULT_INPUT_FILE).as_posix(),
            PwCalculation._DEFAULT_INPUT_FILE
        ))
        folder.get_subfolder(PwCalculation._PSEUDO_SUBFOLDER, create=True)
        remote_copy_list.append((
            self.inputs.parent_folder.computer.uuid,
            Path(scf_path, PwCalculation._PSEUDO_SUBFOLDER, '*').as_posix(),
            PwCalculation._PSEUDO_SUBFOLDER
        ))
        folder.get_subfolder(PwCalculation._OUTPUT_SUBFOLDER, create=True)
        remote_copy_list.append((
            self.inputs.parent_folder.computer.uuid,
            Path(scf_path, PwCalculation._OUTPUT_SUBFOLDER, '*').as_posix(),
            PwCalculation._OUTPUT_SUBFOLDER
        ))
        local_copy_list.append((
            self.inputs.hydrogen_pseudo.uuid,
            self.inputs.hydrogen_pseudo.filename,
            Path(PwCalculation._PSEUDO_SUBFOLDER, self.inputs.hydrogen_pseudo.filename).as_posix()
        ))

        codeinfo = datastructures.CodeInfo()
        codeinfo.cmdline_params = (['-m', 'pynball', self.DEFAULT_INPUT_FILE])
        codeinfo.stdout_name = self.DEFAULT_OUTPUT_FILE
        codeinfo.code_uuid = self.inputs.code.uuid

        calcinfo = datastructures.CalcInfo()
        calcinfo.codes_info = [codeinfo]
        calcinfo.remote_copy_list = remote_copy_list
        calcinfo.local_copy_list = local_copy_list

        calcinfo.retrieve_list = [
            self.DEFAULT_OUTPUT_FILE,
            'output.cif',
            'output.json'
        ]
        return calcinfo
