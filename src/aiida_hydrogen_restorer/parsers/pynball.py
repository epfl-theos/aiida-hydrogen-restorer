# -*- coding: utf-8 -*-
"""`Parser` implementation for the `PynballCalculation` calculation job class."""

import json
from pymatgen.core import Structure

from aiida import orm

from aiida.parsers import Parser
from aiida.common.exceptions import NotExistent


class PynballParser(Parser):
    """`Parser` for the `PynballCalculation` calculation job class."""

    def parse(self, **kwargs):
        """Parse the final structure and energy from the retrieved files."""
        try:
            out_folder = self.retrieved
        except NotExistent:
            return self.exit(self.exit_codes.ERROR_NO_RETRIEVED_FOLDER)

        try:
            with out_folder.open('output.cif', 'r') as handle:
                structure_cif = handle.read()
                structure = Structure.from_str(structure_cif, fmt='cif')
        except OSError:
            return self.exit(self.exit_codes.ERROR_READING_CIF)

        try:
            with out_folder.open('output.json', 'r') as handle:
                output_dict = json.load(handle)
        except OSError:
            return self.exit(self.exit_codes.ERROR_READING_JSON)

        self.out('final_structure', orm.StructureData(pymatgen=structure))
        self.out('output_parameters', orm.Dict(output_dict))
