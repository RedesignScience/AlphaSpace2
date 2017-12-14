************************************
How to select entities in AlphaSpace
************************************



How to load files and select receptor and binder
================================================


This file gives some examples of atom selection in order to manipulate mdtraj objects. Note this is only an example,
and will not run unless you set the path string to a real structure file

Simply loading from files
-------------------------

First you need to initialize the universe object.

::

    import alphaspace
    import mdtraj
    universe = alphaspace.AS_Universe()

If the binder and receptor are in two different pdb file, you can use this command to load them into alphaspace

::

    receptor_structure = mdtraj.load("Path to receptor")
    binder_structure = mdtraj.load("Path to binder")
    universe.set_receptor(structure=receptor_structure)
    universe.set_binder(structure=binder_structure)

Stacking multiple snapshots from different files
------------------------------------------------

If you wishes to add multiple pdb file to the structure and binder as different snapshots, you can use a for loop to
append structures by setting the append option as True

::

    receptor_pdbs = ['pdb1', 'pdb2', 'pdb3']
    binder_pdb = ['pdb1', 'pdb2', 'pdb3']
    for i in range(3):
        receptor_structure = mdtraj.load(receptor_pdbs[i])
        binder_structure = mdtraj.load(binder_pdb[i])
        universe.set_receptor(receptor_structure, append=True)
        universe.set_binder(binder_structure, append=True)


Note that all the structures appended should have the same topology.

Loading from trajectory

If you receptor and binder is stored in a trajectory file, you can load them directly.

::

    receptor_traj = mdtraj.load("path to receptor trajectory", top="path to receptor topology")
    binder_traj = mdtraj.load("path to binder trajectory", top="path to binder topology")
    universe.set_receptor(receptor_traj)
    universe.set_binder(binder_traj)


In a lot of cases, the receptor and binder are in the same file, you can use two different ways to designate them,

1. Specify the residue number of the ligand,
2. Use the build in guess_receptor_binder method (this is less reliable but very convenient if your binder and receptor
    are very clearly seperated as two molecules


Defining receptor and binder by Specify residue number
------------------------------------------------------

::

    trajectory = mdtraj.load('path to trajectory', top='path topology')

If the ligand residue number is 100

::

    binder_traj = alphaspace.extractResidue(traj=trajectory, residue_numbers=[100], clip=True)

If the ligand residue name is 'p53'

::

    binder_traj = alphaspace.extractResidue(traj=trajectory, residue_names=['p53'], clip=True)
    receptor_traj = trajectory  # since the binder has been clipped
    universe.set_receptor(receptor_traj)
    universe.set_binder(binder_traj)



How to select pockets based on lining atoms and residues
========================================================

Here shows you how you can select pockets by screening them against a list of anchoring residues

By Lining Atom
--------------

If you want to screen out certain pockets not containing a given set of atoms, say from number 100 to number 300:

.. code-block:: python

    import alphaspace
    import mdtraj
    import sys

    receptor_path, binder_path = sys.argv[1], sys.argv[2]

    universe = alphaspace.AS_Universe()
    universe.set_receptor(structure=mdtraj.load(receptor_path), keepH=True)
    universe.set_binder(structure=mdtraj.load(binder_path))

    universe.run_alphaspace()


    for pocket in universe.pockets(snapshot_idx=0, active_only=False):
        if not set(pocket.lining_atoms_idx).intersection([100, 101, 102])

    for pocket in universe.pockets(snapshot_idx=0, active_only=False):
        if not set(pocket.lining_residues_idx).intersection(set(range(100,301))):
            pocket.deactivate()

Note this only deactivate the pockets that does not have the given lining residues. Now if you want to iterate through 
pockets that's active, you can:

.. code-block:: python
    
    for pocket in universe.pockets(snapshot_idx=0,active_only=True)

By Lining Residue
-----------------

Let's say you only want pockets that have lining atoms from residue 100,101,102, you can choose to deactivate all
that does not satisfy this criteria by simply calling the lining_residue_idx property .

.. code-block:: python

    for pocket in universe.pockets(snapshot_idx=0, active_only=False):
        if not set(pocket.lining_residues_idx).intersection({100, 101, 102}):
            pocket.deactivate()
    