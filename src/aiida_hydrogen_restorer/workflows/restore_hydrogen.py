# -*- coding: utf-8 -*-
"""Work chain to restore hydrogens to an inputs structure."""

from aiida.engine import ToContext, WorkChain, while_

from aiida import orm
from aiida.common import AttributeDict
from aiida_pseudo.data.pseudo.upf import UpfData

from aiida_quantumespresso.workflows.pw.base import PwBaseWorkChain
from aiida_quantumespresso.calculations.pp import PpCalculation

from aiida_hydrogen_restorer.calculations.add_hydrogens_to_structure import add_hydrogens_to_structure

class RestoreHydrogenWorkChain(WorkChain):

    @classmethod
    def define(cls, spec):
        """Define the process specification"""
        super().define(spec)

        spec.expose_inputs(PwBaseWorkChain, namespace='scf',
            exclude=('clean_workdir', 'pw.structure', 'pw.parent_folder'),
            namespace_options={'help': 'Inputs for the `PwBaseWorkChain` for the initial scf calculation.'})
        spec.expose_inputs(PpCalculation, namespace='pp',
            exclude=('parent_folder'),
            namespace_options={'help': 'Inputs for the `PwBaseWorkChain` for the initial scf calculation.'})

        spec.input('structure', valid_type=orm.StructureData, help='The input structure.')
        spec.input('number_hydrogen', valid_type=orm.Int, help='Number of expected hydrogen in the structure.')
        spec.input('do_supercell', valid_type=orm.Bool, default=lambda: orm.Bool(True), help='If True a supercell 3x3x3 is created.')
        spec.input('equiv_peak_threshold', valid_type=orm.Float, default=lambda: orm.Float(0.995), help='Threshold for selecting maxima peaks.')
        spec.input('hydrogen_pseudo', valid_type=UpfData)
        spec.input('clean_workdir', valid_type=orm.Bool, default=lambda: orm.Bool(False))

        spec.output('final_structure', valid_type=orm.StructureData, help='The final structure.')

        spec.outline(
            cls.run_scf,
            cls.inspect_scf,
            cls.run_pp,
            cls.add_hydrogen,
            while_(cls.more_hydrogen_needed)(
                cls.run_relax,
                cls.inspect_relax,
                cls.run_pp,
                cls.add_hydrogen,
            ),
            #cls.run_final_relax,
            cls.results
        )

        spec.exit_code(401, 'ERROR_SUB_PROCESS_FAILED_SCF',
            message='the `scf` PwBaseWorkChain sub process failed')
        spec.exit_code(402, 'ERROR_SUB_PROCESS_FAILED_PP',
            message='the `pp` PwBaseWorkChain sub process failed')
        spec.exit_code(403, 'ERROR_SUB_PROCESS_FAILED_RELAX',
            message='the `relax` PwBaseWorkChain sub process failed')

    @classmethod
    def get_builder_from_protocol(
        cls,
        pw_code,
        pp_code,
        structure,
        number_hydrogen,
        protocol=None,
        overrides=None, #ok same
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

    def run_scf(self):
        """Run the `PwBaseWorkChain` that calculations the initial potential."""
        self.ctx.current_structure = self.inputs.structure

        inputs = AttributeDict(self.exposed_inputs(PwBaseWorkChain, namespace='scf'))
        inputs.pw.structure = self.inputs.structure

        parameters = inputs.pw.parameters.get_dict()
        parameters['SYSTEM']['tot_charge'] = - (
            self.inputs.number_hydrogen.value - self.ctx.current_structure.get_pymatgen().composition['H']
        )
        inputs.pw.parameters = orm.Dict(parameters)

        running = self.submit(PwBaseWorkChain, **inputs)

        self.report(f'launching PwBaseWorkChain<{running.pk}> for initial scf.')

        return ToContext(workchain_scf=running)

    def inspect_scf(self):
        """Inspect the results of the BLABLA"""
        scf_workchain = self.ctx.workchain_scf
    
        if not scf_workchain.is_finished_ok:
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_SCF
        
        self.ctx.current_folder = scf_workchain.outputs.remote_folder

    def run_pp(self):
        """Run the `PwBaseWorkChain` that calculations the initial potential."""

        inputs = AttributeDict(self.exposed_inputs(PpCalculation, namespace='pp'))
        inputs.parent_folder = self.ctx.current_folder

        pp_calcnode = self.submit(PpCalculation, **inputs)

        self.report(f'launching PwBaseWorkChain<{pp_calcnode.pk}> for initial scf.')

        return ToContext(pp_calculation=pp_calcnode)

    def inspect_pp(self):
        """Inspect the results of the BLABLA"""
        pp_calculation = self.ctx.pp_calculation
    
        if not pp_calculation.is_finished_ok:
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_PP

    def add_hydrogen(self):
        """TODO"""
        structure = self.ctx.current_structure
        potential_array = self.ctx.pp_calculation.outputs.output_data

        results = add_hydrogens_to_structure(
            structure, 
            potential_array,
            self.inputs.do_supercell,
            self.inputs.equiv_peak_threshold,
            self.inputs.number_hydrogen
        )

        self.ctx.current_structure = results['new_structure']

    def more_hydrogen_needed(self):
        """Check if more hydrogen needs to be added to the structure."""
        return self.ctx.current_structure.get_pymatgen().composition['H'] != self.inputs.number_hydrogen.value

    def run_relax(self):
        """"""
        inputs = AttributeDict(self.exposed_inputs(PwBaseWorkChain, namespace='scf'))
        inputs.pw.structure = self.ctx.current_structure
        
        parameters = inputs.pw.parameters.get_dict()
        parameters['CONTROL']['calculation'] = 'relax'
        parameters['SYSTEM']['tot_charge'] = - (
            self.inputs.number_hydrogen.value - self.ctx.current_structure.get_pymatgen().composition['H']
        )
        inputs.pw.parameters = orm.Dict(parameters)

        inputs.pw.pseudos['H'] = self.inputs.hydrogen_pseudo

        settings = inputs.pw.get('settings', {})
        settings['FIXED_COORDS'] = [
            [False, False, False] if site.kind_name == 'H' else [True, True, True]
            for site in self.ctx.current_structure.sites
        ]
        inputs.pw.settings = orm.Dict(settings)

        running = self.submit(PwBaseWorkChain, **inputs)

        self.report(f'launching PwBaseWorkChain<{running.pk}> for relaxation babay.')

        return ToContext(workchain_relax=running)

    def inspect_relax(self):
        """Inspect the results of the relax calc"""
        workchain_relax = self.ctx.workchain_relax
    
        if not workchain_relax.is_finished_ok:
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_RELAX
        
        self.ctx.current_structure = workchain_relax.outputs.output_structure
        self.ctx.current_folder = workchain_relax.outputs.remote_folder

    def results(self):
        structure=self.ctx.current_structure
        self.out('final_structure', structure)

        self.report('Sono alla fine, ma non credo')
