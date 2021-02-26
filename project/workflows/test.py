from __future__ import absolute_import, division, print_function, unicode_literals
from datetime import datetime
from pymatgen.analysis.elasticity.strain import Deformation
from pymatgen.io.vasp.sets import MPRelaxSet, MPStaticSet
from fireworks import Firework, Workflow
from atomate.vasp.fireworks.core import OptimizeFW, TransmuterFW
from atomate.vasp.firetasks.parse_outputs import GibbsFreeEnergyTask
from ase.build import fcc111, add_adsorbate # ASE's utilities to build the surface
from clusterx.parent_lattice import ParentLattice
from clusterx.structures_set import StructuresSet
from clusterx.visualization import juview
from clusterx.super_cell import SuperCell
from random import randint

def wf_clusterx(structure_prim, vasp_input_set_relax=None, vasp_input_set_static=None, vasp_cmd="vasp", db_file=None, user_kpoints_settings=None, nstruc = 60,):
    """
    Args:
        structure_prim (Structure): input structure ASE Structure object.
        vasp_input_set_relax (VaspInputSet)
        vasp_input_set_static (VaspInputSet)
        vasp_cmd (str): vasp command to run.
        db_file (str): path to the db file.
        user_kpoints_settings (dict): example: {"grid_density": 7000}
	
    Returns:
        Workflow
    """
    tag = datetime.utcnow().strftime('%Y-%m-%d-%H-%M-%S-%f')
    # Build a 3-layer Al slab with vacancy in on-top configuration
    pri = fcc111('Al', size=(1,1,3))
    add_adsorbate(pri,'X',1.7,'ontop')
    pri.center(vacuum=10.0, axis=2)

    # Build the parent lattice and a 4x4 supercell
    platt = ParentLattice(pri, site_symbols=[['Al'],['Al'],['Al','Na'],['X','O']])
    platt.get_sublattice_types(pretty_print = True) # Print info regarding the ParentLattice object
    scell = SuperCell(platt,[4,4])

    # Collect 60 random structures in a StructuresSet object
    sset = StructuresSet(platt)
    nstruc = 60

    for i in range(nstruc):
        concentration = {0:[randint(0,4*4)],1:[randint(0,4*4)]} # Pick a random concentration of "Na" substitutions and "O" adsorbants
    	sset.add_structure(scell.gen_random(concentration)) # Generate and add a random structure to the StructuresSet

    print("\nRandom structures (first 3):")
    sset.write_to_db("sset.json") # Write JSON db file for visualization with ASE's GUI.
    


    # get the input set for the optimization and update it if we passed custom settings
    vis_relax = vasp_input_set or MPRelaxSet(structure, force_gamma=True)
    if user_kpoints_settings:
        v = vis_relax.as_dict()
        v.update({"user_kpoints_settings": user_kpoints_settings})
        vis_relax = vis_relax.__class__.from_dict(v)

    # Structure optimization firework
    fws = [OptimizeFW(structure=structure, vasp_input_set=vis_relax, vasp_cmd=vasp_cmd,
                      db_file=db_file, name="{} structure optimization".format(tag))]

    fws.append(fw_analysis)

    # finally, create the workflow
    wf_gibbs = Workflow(fws)
    wf_gibbs.name = "{}:{}".format(structure.composition.reduced_formula, "gibbs free energy")

    return wf_gibbs
