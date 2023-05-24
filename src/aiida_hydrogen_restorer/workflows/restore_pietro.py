# -*- coding: utf-8 -*-
"""Work chain to restore hydrogens to an inputs structure."""

from aiida.engine import ToContext, WorkChain, while_, calcfunction
from aiida import orm
from aiida.common import AttributeDict
from aiida_pseudo.data.pseudo.upf import UpfData
from aiida_quantumespresso.workflows.pw.base import PwBaseWorkChain
from aiida_quantumespresso.calculations.pp import PpCalculation

from aiida_hydrogen_restorer.calculations.add_hydrogens_to_structure import add_hydrogens_to_structure

@calcfunction
def subtract_potentials(array_1, array_2):
    difference = array_1.get_array('data') - array_2.get_array('data') 
    potential_difference = orm.ArrayData()
    potential_difference.set_array('data', difference)
    return {'potential_difference': potential_difference}


class RestorePietroWorkChain(WorkChain):

    @classmethod
    def define(cls, spec):
        """Define the process specification"""
        super().define(spec)

        spec.expose_inputs(PwBaseWorkChain, namespace='scf',
            exclude=('clean_workdir', 'pw.structure', 'pw.parent_folder'),
            namespace_options={'help': 'Inputs for the `PwBaseWorkChain` for the initial scf calculation.'})
        spec.expose_inputs(PpCalculation, namespace='pp',
            exclude=('parent_folder'),
            namespace_options={'help': 'Inputs for the `pp.x` process to find electrostatic potential.'})

        spec.input('structure', valid_type=orm.StructureData, help='The input structure.')
        spec.input('number_hydrogen', valid_type=orm.Int, help='Number of expected hydrogen in the structure.')
        spec.input('do_supercell', valid_type=orm.Bool, default=lambda: orm.Bool(True), help='If True a supercell 3x3x3 is created.')
        spec.input('equiv_peak_threshold', valid_type=orm.Float, default=lambda: orm.Float(0.995), help='Threshold for selecting maxima peaks.')
        spec.input('hydrogen_pseudo', valid_type=UpfData)
        spec.input('clean_workdir', valid_type=orm.Bool, default=lambda: orm.Bool(False))
        spec.output('all_peaks', valid_type=orm.ArrayData, help='List of the maxima peaks')
        spec.output('final_structure', valid_type=orm.StructureData, help='The final structure.')

        spec.outline(
            cls.setup,
            while_(cls.should_add_hydrogens)(
                cls.run_scf,
                cls.inspect_scf,
                cls.run_pp,
                cls.inspect_pp,
                cls.add_hydrogen,
            ),
            cls.results
        )

        spec.exit_code(401, 'ERROR_SUB_PROCESS_FAILED_SCF',
            message='the `scf` PwBaseWorkChain sub process failed')
        spec.exit_code(402, 'ERROR_SUB_PROCESS_FAILED_PP',
            message='the `pp` PwBaseWorkChain sub process failed')
        spec.exit_code(403, 'ERROR_SUB_PROCESS_FAILED_RELAX',
            message='the `relax` PwBaseWorkChain sub process failed')
        spec.exit_code(501, 'WARNING_FINAL_STRUCTURE_NOT_COMPLETE',
            message='the final obtained structure does not have the required number of hydrogen.')

    @classmethod
    def get_builder_from_protocol(
        cls,
        pw_code,
        pp_code,
        structure,
        number_hydrogen,
        protocol=None,
        overrides=None,
        **kwargs
    ):
        overrides = {} if overrides is None else overrides
        base_inputs = PwBaseWorkChain.get_protocol_inputs(protocol, overrides.get('scf', None))
        pseudo_family = orm.load_group(base_inputs.pop('pseudo_family'))

        base_scf = PwBaseWorkChain.get_builder_from_protocol(
            code=pw_code, structure=structure, protocol=protocol, overrides=overrides.get('scf', None)
        )
        base_scf['pw'].pop('structure', None)
        base_scf.pop('clean_workdir', None)

        builder = cls.get_builder()

        builder.scf = base_scf

        pp_builder = PpCalculation.get_builder()
        pp_builder.code = pp_code
        parameters = {
            'INPUTPP': {
                'plot_num': 11,
            },
            'PLOT': {
                'iflag': 3,
            },
        }
        pp_builder.metadata.options.resources = {
            'num_machines': 1,
            'num_mpiprocs_per_machine': 1,
            'num_cores_per_machine': 1,
        }
        pp_builder.metadata.options.max_wallclock_seconds = 1800
        pp_builder.parameters = orm.Dict(parameters)
        builder.pp = pp_builder
        builder.structure = structure
        builder.number_hydrogen = orm.Int(number_hydrogen)
        builder.hydrogen_pseudo = pseudo_family.get_pseudo('H')

        return builder

    def setup(self):
        """Set up the initial context variables."""
        self.ctx.current_structure = self.inputs.structure
        self.ctx.current_folder = None
        self.ctx.failed_to_add_hydrogen = False
        self.ctx.all_peaks = None
        self.ctx.num_peaks = None

    def run_scf(self):
        """Run the `PwBaseWorkChain` that calculations the initial potential."""
        structure = self.ctx.current_structure

        # Full SCF with electronic charge equal to number of missing hydrogen
        full_inputs = AttributeDict(self.exposed_inputs(PwBaseWorkChain, namespace='scf'))
        full_inputs.pw.structure = structure
        parameters = full_inputs.pw.parameters.get_dict()
        parameters['SYSTEM']['tot_charge'] = - (
            self.inputs.number_hydrogen.value - self.ctx.current_structure.get_pymatgen().composition['H']
        )
        full_inputs.pw.parameters = orm.Dict(parameters)
        base_full = self.submit(PwBaseWorkChain, **full_inputs)

        # Partial SCF with electronic charge equal to number of missing hydrogen minus one
        partial_inputs = AttributeDict(self.exposed_inputs(PwBaseWorkChain, namespace='scf'))
        partial_inputs.pw.structure = structure
        partial_params = partial_inputs.pw.parameters.get_dict()
        partial_params['SYSTEM']['tot_charge'] = - (
            self.inputs.number_hydrogen.value - self.ctx.current_structure.get_pymatgen().composition['H']
        ) + 1
        partial_inputs.pw.parameters = orm.Dict(partial_params)
        base_partial = self.submit(PwBaseWorkChain, **partial_inputs)

        self.report(f'launched two PwBaseWorkChain for initial scf: {base_full.pk} & {base_partial.pk}')

        return ToContext(base_full=base_full, base_partial=base_partial)

    def inspect_scf(self):
        """Inspect the results of both SCF runs."""
        full_scf_workchain = self.ctx.base_full
        partial_scf_workchain = self.ctx.base_partial

        if any((
            not full_scf_workchain.is_finished_ok,
            not partial_scf_workchain.is_finished_ok,
        )):
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_SCF

        self.ctx.current_full_folder = full_scf_workchain.outputs.remote_folder
        self.ctx.current_partial_folder = partial_scf_workchain.outputs.remote_folder

    def run_pp(self):
        """Run the `PwBaseWorkChain` that calculations the initial potential."""
        full_inputs = AttributeDict(self.exposed_inputs(PpCalculation, namespace='pp'))
        full_inputs.parent_folder = self.ctx.current_full_folder
        pp_calc_full = self.submit(PpCalculation, **full_inputs)

        partial_inputs = AttributeDict(self.exposed_inputs(PpCalculation, namespace='pp'))
        partial_inputs.parent_folder = self.ctx.current_partial_folder
        pp_calc_partial = self.submit(PpCalculation, **partial_inputs)

        self.report(
            'launched two `pp.x` calculations to find electrostatic potential with PKs:'
            f'<{pp_calc_full.pk}> and <{pp_calc_partial.pk}>'
        )
        return ToContext(pp_calc_full=pp_calc_full, pp_calc_partial=pp_calc_partial)

    def inspect_pp(self):
        """Inspect the results of the BLABLA"""
        pp_calc_full = self.ctx.pp_calc_full
        pp_calc_partial = self.ctx.pp_calc_partial
    
        if any((
            not pp_calc_full.is_finished_ok,
            not pp_calc_partial.is_finished_ok,
        )):
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_PP

    def add_hydrogen(self):
        """Add hydrogen to the current structure based on the Pietro method."""
        structure = self.ctx.current_structure
        potential_array_full = self.ctx.pp_calc_full.outputs.output_data
        potential_array_partial = self.ctx.pp_calc_partial.outputs.output_data

        potential_difference = subtract_potentials(
            array_1 = potential_array_full,
            array_2 = potential_array_partial
        )['potential_difference']

        results = add_hydrogens_to_structure(
            structure, 
            potential_difference,
            self.inputs.do_supercell,
            self.inputs.equiv_peak_threshold,
            self.inputs.number_hydrogen
        )
        if structure.get_pymatgen().composition['H'] == results['new_structure'].get_pymatgen().composition['H']:
            self.ctx.failed_to_add_hydrogen = True
        else:
            self.ctx.current_structure = results['new_structure']
            self.ctx.all_peaks = results['all_peaks']
            self.ctx.num_peaks = len(results['all_peaks'].get_array('peak_values'))
            current_H = self.ctx.current_structure.get_pymatgen().composition['H']
            self.report(
                f'Now there are {current_H} out of {self.inputs.number_hydrogen.value} hydrogens '
                f'(I found {self.ctx.num_peaks} maxima).'
            )

    def should_add_hydrogens(self):
        """Check if more hydrogens should be added to the structure."""
        not_enough_hydrogen = (
            self.ctx.current_structure.get_pymatgen().composition['H'] != self.inputs.number_hydrogen.value
        )
        return not_enough_hydrogen and not self.ctx.failed_to_add_hydrogen == True

    def results(self):
        """Add the results to the outputs."""
        structure=self.ctx.current_structure
        all_peaks = self.ctx.all_peaks
        enough_hydrogen = (
            self.ctx.current_structure.get_pymatgen().composition['H'] == self.inputs.number_hydrogen.value
        )
        self.out('all_peaks', all_peaks)
        self.out('final_structure', structure)

        if enough_hydrogen:
            self.report('Good job!')
        else:
            self.report('You need to change method.')
            return self.exit_codes.WARNING_FINAL_STRUCTURE_NOT_COMPLETE
