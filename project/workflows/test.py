from __future__ import absolute_import, division, print_function, unicode_literals
from datetime import datetime
from pymatgen.analysis.elasticity.strain import Deformation
from pymatgen.io.vasp.sets import MPRelaxSet, MPStaticSet
from fireworks import Firework, Workflow
from atomate.vasp.fireworks.core import OptimizeFW, TransmuterFW
from atomate.vasp.firetasks.parse_outputs import GibbsFreeEnergyTask


def wf_gibbs_free_energy(structure, deformations, vasp_input_set_relax=None, vasp_input_set_static=None, vasp_cmd="vasp",
                             db_file=None, user_kpoints_settings=None, t_step=10, t_min=0, t_max=1000,
                             mesh=(20, 20, 20), eos="vinet", qha_type="debye_model", pressure=0.0,
                             poisson=0.25):
    """
    Returns quasi-harmonic gibbs free energy workflow.
    Note: phonopy package is required for the final analysis step if qha_type="phonopy"
    Args:
        structure (Structure): input structure.
        deformations (list): list of deformation matrices(list of lists).
        vasp_input_set_relax (VaspInputSet)
        vasp_input_set_static (VaspInputSet)
        vasp_cmd (str): vasp command to run.
        db_file (str): path to the db file.
        user_kpoints_settings (dict): example: {"grid_density": 7000}
        t_step (float): temperature step (in K)
        t_min (float): min temperature (in K)
        t_max (float): max temperature (in K)
        mesh (list/tuple): reciprocal space density
        eos (str): equation of state used for fitting the energies and the volumes.
            options supported by phonopy: "vinet", "murnaghan", "birch_murnaghan".
            Note: pymatgen supports more options than phonopy. see pymatgen.analysis.eos.py
        qha_type(str): quasi-harmonic approximation type: "debye_model" or "phonopy",
            default is "debye_model"
        pressure (float): in GPa
        poisson (float): poisson ratio
    Returns:
        Workflow
    """
    tag = datetime.utcnow().strftime('%Y-%m-%d-%H-%M-%S-%f')

    # get the input set for the optimization and update it if we passed custom settings
    vis_relax = vasp_input_set or MPRelaxSet(structure, force_gamma=True)
    if user_kpoints_settings:
        v = vis_relax.as_dict()
        v.update({"user_kpoints_settings": user_kpoints_settings})
        vis_relax = vis_relax.__class__.from_dict(v)

    # Structure optimization firework
    fws = [OptimizeFW(structure=structure, vasp_input_set=vis_relax, vasp_cmd=vasp_cmd,
                      db_file=db_file, name="{} structure optimization".format(tag))]

    # get the input set for the static calculations and update it if we passed custom settings
    uis_static = {"ISIF": 2, "ISTART":1}
    lepsilon = False # use IBRION = -1; don't update the ions
    if qha_type not in ['debye model']:
        uis_static = {'ISIF'}
        lepsilon = True # use IBRION = -8; DFPT
    vis_static = MPStaticSet(structure, force_gamma=True, lepsilon=lepsilon,
                         user_kpoints_settings=user_kpoints_settings,
                         user_incar_settings=uis_static)

    # create each deformation Firework and add them to the Fireworks list
    parents = fws[0]
    deformations = [Deformation(defo_mat) for defo_mat in deformations]
    for n, deformation in enumerate(deformations):
        fw = TransmuterFW(name="{} {} {}".format(tag, 'gibbs deformation', n), structure=structure,
                          transformations=['DeformStructureTransformation'],
                          transformation_params=[{"deformation": deformation.tolist()}],
                          vasp_input_set=vis_static, parents=parents,
                          vasp_cmd=vasp_cmd, db_file=db_file)
        fws.append(fw)

    parents = fws[1:] # all of the deformation Fireworks
    if qha_type not in ["debye_model"]:
        from phonopy import Phonopy
    fw_analysis = Firework(GibbsFreeEnergyTask(tag=tag, db_file=db_file, t_step=t_step, t_min=t_min,
                                               t_max=t_max, mesh=mesh, eos=eos, qha_type=qha_type,
                                               pressure=pressure, poisson=poisson),
                           name="gibbs free energy", parents=parents)
    fws.append(fw_analysis)

    # finally, create the workflow
    wf_gibbs = Workflow(fws)
    wf_gibbs.name = "{}:{}".format(structure.composition.reduced_formula, "gibbs free energy")

    return wf_gibbs
