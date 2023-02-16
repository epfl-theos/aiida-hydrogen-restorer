# -*- coding: utf-8 -*-
"""Work chain to restore hydrogens to an inputs structure."""

from aiida.engine import ToContext, WorkChain, while_

from aiida import orm
from aiida.common import AttributeDict
from aiida_quantumespresso.workflows.pw.base import PwBaseWorkChain
from aiida_quantumespresso.calculations.pp import PpCalculation


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

    @classmethod
    def get_builder_from_protocol(
        cls,
        pw_code,
        pp_code,
        structure,
        number_hydrogen,
        #pp_builder,
        protocol=None,
        overrides=None,
        **kwargs
    ):
        base_scf = PwBaseWorkChain.get_builder_from_protocol(
            code=pw_code, structure=structure, protocol=protocol
        )
        base_scf['pw'].pop('structure', None)
        base_scf.pop('clean_workdir', None)

        builder = cls.get_builder()

        builder.scf = base_scf

        pp_builder = PpCalculation.get_builder()
        pp_builder.code = pp_code
        parameters = {
            'INPUTPP': {
                'plot_num': 11,  # electron charge density ===> what if I want also the number 11?
            },
            'PLOT': {
                'iflag': 3,  # 3D plot
                #'output_format': 6,  # format as gaussian cube file  (3D)
            },
        }
        pp_builder.metadata.options.resources = {
            'num_machines': 1,
            'num_mpiprocs_per_machine': 8,
            'num_cores_per_machine': 8,
        }
        pp_builder.metadata.options.max_wallclock_seconds = 1800
        pp_builder.parameters = orm.Dict(parameters)
        builder.pp = pp_builder
        builder.structure = structure
        builder.number_hydrogen = orm.Int(number_hydrogen)

        return builder

    def run_scf(self):
        """Run the `PwBaseWorkChain` that calculations the initial potential."""

        inputs = AttributeDict(self.exposed_inputs(PwBaseWorkChain, namespace='scf'))
        inputs.pw.structure = self.inputs.structure

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

    def add_hydrogen(self):
        """TODO"""
        pass

#         structure = self.inputs.structure
#         potential_array = self.ctx.pp_calculation.outputs.output_data

#         results = add_hydrogens_to_structure(
#             structure, 
#             potential_array,
#             self.inputs.do_supercell,
#             self.inputs.equiv_peak_threshold,
#             self.inputs.number_hydrogen
#         )

#         #structure = self.out('output_structure', results['new_structure'])
#         self.ctx.current_structure = results['new_structure']
#         self.out('final_structure', results['new_structure'])

     def more_hydrogen_needed(self):
        """Check if more hydrogen needs to be added to the structure."""
        pass
#         structure = self.ctx.current_structure
#         return structure.get_pymatgen().composition['H'] != self.inputs.number_hydrogen.value
#         #return False

#     def run_relax(self):
#         from aiida_quantumespresso.common.types import RelaxType
#         """Run the PwRelaxWorkChain to run a relax PwCalculation."""
#         #inputs = AttributeDict(self.exposed_inputs(PwRelaxWorkChain))#, namespace='relax')) 

#          #inputs.metadata.call_link_label = 'relax'
#          #inputs.structure = self.ctx.current_structure

#         overrides = {
#             'pseudo_family': 'PseudoDojo/0.4/PBE/SR/standard/upf'
#         }
#         code = orm.load_code('qe-7.0-pw@eiger')
#         #structure=self.ctx.current_structure,

#         builder = PwRelaxWorkChain.get_builder_from_protocol(
#             code=code,
#             structure=self.ctx.current_structure,
#             overrides=overrides,
#             relax_type=RelaxType.POSITIONS
#         )
#         builder.base.pw.parameters['SYSTEM']['tot_charge'] = -(self.inputs.number_hydrogen.value - builder.structure.get_pymatgen().composition['H'])
#         builder.base.pw.metadata.options.resources = {
#             'num_machines': 4,
#             'num_mpiprocs_per_machine': 1,
#             'num_cores_per_mpiproc': 12,
#         }
#         builder.clean_workdir = orm.Bool(False)
#         builder.base.pw.parallelization = orm.Dict({'npool': 4})

#         fixed_coords = []

#         for site in builder.structure.sites:
#             if site.kind_name == 'H':
#                 fixed_coords.append([False, False, False])
#             else:
#                 fixed_coords.append([True, True, True])

#  ## CREATE LIST THAT SETS TO 'FALSE' ALL Hs COORDINATES
#         builder.base.pw.settings = orm.Dict({'FIXED_COORDS': fixed_coords})
#         builder.base.pw.metadata.options['account'] = 'mr32'

#         builder.base.pw.parameters['CONTROL']['nstep'] = 250
#         builder.base.pw.parameters['IONS'] = {
#             'ion_dynamics': 'damp'
#         }
#         builder.base.max_iterations = orm.Int(1)


#         # inputs.base.pw.metadata.options.resources = {
#         #     'num_machines': 4,
#         #     'num_mpiprocs_per_machine': 1,
#         #     'num_cores_per_mpiproc': 12,
#         # }

#         running = self.submit(builder)

#         self.report(f'launching PwRelaxWorkChain<{running.pk}> to relax H positions.')

#         return ToContext(workchain_relax=running)

     def inspect_relax(self):
#         """Inspect the results of the relax calc"""
#         relax_workchain = self.ctx.workchain_relax
    
#         if not relax_workchain.is_finished_ok:
#             return self.exit_codes.ERROR_SUB_PROCESS_FAILED_RELAX
        
#         self.ctx.current_folder = relax_workchain.outputs.remote_folder
        pass 

    # def move(self):

    #     pass   

#     def run_final_relax(self):
#         from aiida_quantumespresso.common.types import RelaxType
#         """Run the PwRelaxWorkChain to run a relax PwCalculation."""
#         #inputs = AttributeDict(self.exposed_inputs(PwRelaxWorkChain))#, namespace='relax')) 

#          #inputs.metadata.call_link_label = 'relax'
#          #inputs.structure = self.ctx.current_structure

#         overrides = {
#             'pseudo_family': 'PseudoDojo/0.4/PBE/SR/standard/upf'
#         }
#         code = orm.load_code('qe-7.0-pw@daint-gpu')
#         structure=self.ctx.current_structure,

#         builder = PwRelaxWorkChain.get_builder_from_protocol(
#             code=code,
#             structure=structure,
#             overrides=overrides,
#             relax_type=RelaxType.POSITIONS
#         )
        
#         builder.base.pw.metadata.options.resources = {
#             'num_machines': 4,
#             'num_mpiprocs_per_machine': 1,
#             'num_cores_per_mpiproc': 12,
#         }
#         builder.clean_workdir = orm.Bool(False)
#         builder.base.pw.parallelization = orm.Dict({'npool': 4})

#         fixed_coords = []

#         for site in builder.structure.sites:
#             if site.kind_name == 'H':
#                 fixed_coords.append([False, False, False])
#             else:
#                 fixed_coords.append([True, True, True])

#  ## CREATE LIST THAT SETS TO 'FALSE' ALL Hs COORDINATES
#         builder.base.pw.settings = orm.Dict({'FIXED_COORDS': fixed_coords})

#         builder.base.pw.parameters['CONTROL']['nstep'] = 250
#         builder.base.pw.parameters['IONS'] = {
#             'ion_dynamics': 'damp'
#         }
#         builder.base.max_iterations = orm.Int(1)


#         # inputs.base.pw.metadata.options.resources = {
#         #     'num_machines': 4,
#         #     'num_mpiprocs_per_machine': 1,
#         #     'num_cores_per_mpiproc': 12,
#         # }

#         running = self.submit(builder)

#         self.report(f'launching PwRelaxWorkChain<{running.pk}> to relax H positions.')

#         return ToContext(workchain_relax=running)


    def results(self):
        structure=self.ctx.current_structure,
        self.out('final_structure', structure)
        #self.ctx.current_folder = relax_workchain.outputs.remote_folder

        #self.report('Sono alla fine, ma non credo')