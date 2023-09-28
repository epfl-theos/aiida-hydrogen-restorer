# -*- coding: utf-8 -*-
"""Work chain to restore hydrogens to an inputs structure."""

from aiida.engine import ToContext, WorkChain, while_, if_, calcfunction
from aiida import orm
from aiida.common import AttributeDict
from aiida_pseudo.data.pseudo.upf import UpfData
from aiida_quantumespresso.workflows.pw.base import PwBaseWorkChain
from aiida_quantumespresso.calculations.pp import PpCalculation
from aiida_hydrogen_restorer.calculations.pynball import PynballCalculation
from qe_tools import CONSTANTS


from aiida_hydrogen_restorer.calculations.add_hydrogens_to_structure import add_hydrogens_to_structure

@calcfunction
def get_energy(energy):
    initial_energy = energy / CONSTANTS.ry_to_ev
    return orm.Float(initial_energy)


class RestoreHydrogenPWorkChain(WorkChain):

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
        spec.expose_inputs(PynballCalculation, namespace='pinball',
            exclude=('parent_folder', 'all_peaks', 'number_hydrogen'),
            namespace_options={'help': 'Inputs for the `pinball.x` process .'})
        


        spec.input('structure', valid_type=orm.StructureData, help='The input structure.')
        spec.input('number_hydrogen', valid_type=orm.Int, help='Number of expected hydrogen in the structure.')
        spec.input('do_supercell', valid_type=orm.Bool, default=lambda: orm.Bool(True), help='If True a supercell 3x3x3 is created.')
        spec.input('equiv_peak_threshold', valid_type=orm.Float, default=lambda: orm.Float(0.995), help='Threshold for selecting maxima peaks.')
        spec.input('hydrogen_pseudo', valid_type=UpfData)
        spec.input('clean_workdir', valid_type=orm.Bool, default=lambda: orm.Bool(True))
        spec.output('all_peaks', valid_type=orm.ArrayData, help='List of the maxima peaks') 
        spec.output('final_structure', valid_type=orm.StructureData, help='The final structure.')
        spec.output('initial_energy', valid_type=orm.Float, help='The energy of the input structure with H.')


        spec.outline(
            cls.setup,
            cls.run_initial_scf,
            cls.run_scf,
            cls.inspect_scf,
            cls.run_pp,
            cls.inspect_pp,
            cls.add_hydrogen,
            while_(cls.should_add_hydrogens)(
                cls.run_relax_hydrogens,
                cls.inspect_relax,
                cls.run_pp,
                cls.inspect_pp,
                cls.add_hydrogen,
                ),
            cls.pinball_is_needed, #this should be changed into an if statement, which now is included in the function
            cls.inspect_pinball,
            cls.run_relax_hydrogens, 
            cls.inspect_relax,
            cls.results
        )

        spec.exit_code(401, 'ERROR_SUB_PROCESS_FAILED_SCF',
            message='the `scf` PwBaseWorkChain sub process failed')
        spec.exit_code(402, 'ERROR_SUB_PROCESS_FAILED_PP',
            message='the `pp` PwBaseWorkChain sub process failed')
        spec.exit_code(403, 'ERROR_SUB_PROCESS_FAILED_RELAX',
            message='the `relax` PwBaseWorkChain sub process failed')
        spec.exit_code(501, 'WARNING_FINAL_STRUCTURE_NOT_COMPLETE', #this one shouldn't occur anymore
            message='the final obtained structure does not have the required number of hydrogen.')
        spec.exit_code(502, 'ERROR_SUB_PROCESS_FAILED_PINBALL',
            message='the `pinball` process failed')
        spec.exit_code(503, 'NEW_SITE_TOO_CLOSE_TO_EXISTING_ONE',
            message='the final obtained structure does not have the required number of hydrogen.')



    @classmethod
    def get_builder_from_protocol(
        cls,
        pw_code,
        pp_code,
        pinball_code,
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

        pinball_builder = PynballCalculation.get_builder()
        pinball_builder.code = pinball_code
        pinball_builder.metadata.options.resources = {
            'num_machines': 1,
            'num_mpiprocs_per_machine': 1,
            'num_cores_per_machine': 1,
        }
        pinball_builder.metadata.options.max_wallclock_seconds = 1800
        pinball_builder.hydrogen_pseudo = pseudo_family.get_pseudo('H')
        builder.pinball = pinball_builder

        return builder

    def setup(self):
        """Set up the initial context variables."""
        self.ctx.current_structure = self.inputs.structure
        self.ctx.current_folder = None
        self.ctx.failed_to_add_hydrogen = False
        self.ctx.all_peaks = None
        self.ctx.num_peaks = None

    def run_initial_scf(self):
        """Run the `PwBaseWorkChain` that calculations the energy for the reference structure."""
        structure_uuid = self.ctx.current_structure.extras['uuid_original_structure_withH']
        structure = orm.load_node(uuid=structure_uuid)

        inputs = AttributeDict(self.exposed_inputs(PwBaseWorkChain, namespace='scf'))
        inputs.pw.structure = structure
        inputs.pw.pseudos['H'] = self.inputs.hydrogen_pseudo

        parameters = inputs.pw.parameters.get_dict()
        
        inputs.pw.parameters = orm.Dict(parameters)

        pw_base_node = self.submit(PwBaseWorkChain, **inputs)
        self.report(f'launching PwBaseWorkChain<{pw_base_node.pk}> for scf on the reference structure.')

        return ToContext(workchain_scf_initialstructure=pw_base_node)        

    def run_scf(self):
        """Run the `PwBaseWorkChain` that calculates the initial potential."""
        structure = self.ctx.current_structure

        inputs = AttributeDict(self.exposed_inputs(PwBaseWorkChain, namespace='scf'))
        inputs.pw.structure = structure

        parameters = inputs.pw.parameters.get_dict()
        parameters['SYSTEM']['tot_charge'] = - (
            self.inputs.number_hydrogen.value - structure.get_pymatgen().composition['H']
        )
        inputs.pw.parameters = orm.Dict(parameters)

        pw_base_node = self.submit(PwBaseWorkChain, **inputs)
        self.report(f'launching PwBaseWorkChain<{pw_base_node.pk}> for initial scf.')

        return ToContext(workchain_scf=pw_base_node)

    def inspect_scf(self):
        """Inspect the results of the scf `PwBaseWorkChain`."""
        scf_workchain = self.ctx.workchain_scf

        if not scf_workchain.is_finished_ok:
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_SCF

        self.ctx.current_folder = scf_workchain.outputs.remote_folder

    def run_pp(self):
        """Run the `PwBaseWorkChain` that calculations the initial potential."""
        inputs = AttributeDict(self.exposed_inputs(PpCalculation, namespace='pp'))
        inputs.parent_folder = self.ctx.current_folder

        pp_calc_node = self.submit(PpCalculation, **inputs)
        self.report(f'launching pp.x <{pp_calc_node.pk}> to find electrostatic potential.')
        
        return ToContext(pp_calculation=pp_calc_node)

    def inspect_pp(self):
        """Inspect the results of the `PpCalculation`"""
        pp_calculation = self.ctx.pp_calculation
    
        if not pp_calculation.is_finished_ok:
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_PP

    def add_hydrogen(self):
        """Add hydrogen to the current structure."""
        structure = self.ctx.current_structure
        potential_array = self.ctx.pp_calculation.outputs.output_data

        try:
            results = add_hydrogens_to_structure(
                structure, 
                potential_array,
                self.inputs.do_supercell,
                self.inputs.equiv_peak_threshold,
                self.inputs.number_hydrogen
            )
        except ValueError:
            return self.exit_codes.NEW_SITE_TOO_CLOSE_TO_EXISTING_ONE
        
        if structure.get_pymatgen().composition['H'] == results['new_structure'].get_pymatgen().composition['H']:
            self.ctx.failed_to_add_hydrogen = True
            self.ctx.all_peaks = results['all_peaks']

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

    def pinball_is_needed(self):
        """Check if more hydrogens should be added to the structure."""
        if self.ctx.failed_to_add_hydrogen == True:

            inputs = AttributeDict(self.exposed_inputs(PynballCalculation, namespace='pinball'))
            inputs.parent_folder = self.ctx.current_folder # it should be the latest pw folder
            inputs.all_peaks = self.ctx.all_peaks
            inputs.number_hydrogen = orm.Int(self.inputs.number_hydrogen.value - self.ctx.current_structure.get_pymatgen().composition['H'])
            pinball_calc_node = self.submit(PynballCalculation, **inputs)
            self.report(f'launching pinball.x <{pinball_calc_node.pk}>.')
        
            return ToContext(pinball_calculation=pinball_calc_node)

        else:
            pass
    
    def inspect_pinball(self):
            
        if self.ctx.failed_to_add_hydrogen == True:
            """Inspect the results of the pinball calc"""
            pinball_calculation = self.ctx.pinball_calculation

            if not pinball_calculation.is_finished_ok:
                return self.exit_codes.ERROR_SUB_PROCESS_FAILED_PINBALL
            
            self.ctx.current_structure = pinball_calculation.outputs.final_structure
            self.ctx.current_folder = pinball_calculation.outputs.remote_folder 

        else: 
            pass

    def run_relax_hydrogens(self):
        """Run the relaxation for new structure."""

    # if self.ctx.failed_to_add_hydrogen == True or self.ctx.current_structure.get_pymatgen().composition['H'] == 0 : 
    #     pass
    # else:
        inputs = AttributeDict(self.exposed_inputs(PwBaseWorkChain, namespace='scf'))
        inputs.pw.structure = self.ctx.current_structure

        parameters = inputs.pw.parameters.get_dict()
        parameters['CONTROL']['calculation'] = 'relax'
        parameters['SYSTEM']['tot_charge'] = - (
            self.inputs.number_hydrogen.value - self.ctx.current_structure.get_pymatgen().composition['H']
        )
        parameters['CONTROL']['nstep'] = 250
        parameters['IONS'] = {'ion_dynamics': 'damp'}
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
            
        # if self.ctx.failed_to_add_hydrogen == True or self.ctx.current_structure.get_pymatgen().composition['H'] == 0 : 
        #     pass
        # else: 
        """Inspect the results of the relax calc"""
        workchain_relax = self.ctx.workchain_relax

        if not workchain_relax.is_finished_ok:
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_RELAX
        
        self.ctx.current_structure = workchain_relax.outputs.output_structure
        self.ctx.current_folder = workchain_relax.outputs.remote_folder

    def results(self):
        """Add the results to the outputs."""
        structure=self.ctx.current_structure
        all_peaks = self.ctx.all_peaks
        energy = self.ctx.workchain_scf_initialstructure.outputs.output_parameters.get_dict()['energy']
        initial_energy = get_energy(energy)

        self.ctx.enough_hydrogen = (
            self.ctx.current_structure.get_pymatgen().composition['H'] == self.inputs.number_hydrogen.value
        )
        self.out('all_peaks', all_peaks)
        self.out('final_structure', structure)
        self.out('initial_energy', initial_energy)

        if self.ctx.enough_hydrogen:
            self.report('Good job!')
        else:
            self.report('You need to change method.')
            return self.exit_codes.WARNING_FINAL_STRUCTURE_NOT_COMPLETE

    def on_terminated(self):
        """Clean the working directories of all child calculations if `clean_workdir=True` in the inputs."""
        super().on_terminated()

        if self.inputs.clean_workdir.value is False or not self.ctx.enough_hydrogen:
            self.report('remote folders will not be cleaned')
            return

        cleaned_calcs = []

        for called_descendant in self.node.called_descendants:
            if isinstance(called_descendant, orm.CalcJobNode):
                try:
                    called_descendant.outputs.remote_folder._clean()  # pylint: disable=protected-access
                    cleaned_calcs.append(called_descendant.pk)
                except (IOError, OSError, KeyError):
                    pass

        if cleaned_calcs:
            self.report(f"cleaned remote folders of calculations: {' '.join(map(str, cleaned_calcs))}")
