import os
from os import listdir
from os.path import isfile, join
import numpy as np
import regex
import subprocess
from subprocess import DEVNULL
import sys

# Get input for CPU usage
cpus = sys.argv[1]
nodes = sys.argv[2]

# Get all structure files
cwd = os.getcwd()
folder = 'structure_files/slab_7'
structures = join(cwd, folder)
structure_files = [f for f in listdir(structures) if isfile(join(structures, f))]

for name in structure_files:
    # Create temp stucture files to pull from, allows for shorter relaxation times
    process = subprocess.Popen("cp {_folder}/{_name} {_folder}_relaxed/{_name}_temp".format(_folder=folder, _name=name), shell=True, stdout=DEVNULL)
    process.wait()
    
    # Loop through different external fields
    e_ext_array = np.arange(0.005, 0.085, 0.005)
    for e_ext in e_ext_array:

        # Remove data from structure file
        structure_file = open("{_folder}_relaxed/{_name}_temp".format(_folder=folder, _name=name), 'r')
        lines = structure_file.readlines()
        structure_file.close()
        
        # Extract data and process it into useful strings
        lattice_x = lines[2]
        lattice_y = lines[3]
        lattice_z = lines[4]
        
        lattice = ''.join([lattice_x, lattice_y, lattice_z])
 
        atomic_positions = lines[8:]
        _temp = []
        for i in atomic_positions:
            line = i.split()
            _temp.append(' '.join(line[3:] + line[:3]))
        atomic_positions = '\n'.join(_temp)
        
        typs = lines[5].split()
        number = lines[6].split()
        d_count = dict()
 
        for i,typ in enumerate(typs):
            if typ in d_count:
                d_count[typ] = str(int(d_count[typ]) + int(number[i]))
            else:
                d_count[typ] = number[i]
 
        ntyp = len(d_count)
        nat = sum(list(map(int, number)))
 
        masses_file = open('atomic_weights', 'r')
        lines = masses_file.readlines()
        masses_file.close()
        masses = []
        for line in lines:
            split_line = line.split()
            if split_line[0] in d_count:
                masses.append(split_line[0] + ' ' + split_line[1] + ' ' + split_line[0] + '.pbesol-nc-sr.upf')
        masses = '\n'.join(masses)


        # Write data to input file
        text = """
&CONTROL
  calculation=vc-relax
  prefix={_prefix}
  outdir='./out/'
  pseudo_dir='./pseudo/'

  forc_conv_thr=1.0d-3

  tefield=.TRUE.
  dipfield=.TRUE.
/

&SYSTEM
  occupations='smearing'
  smearing='gaussian'
  degauss=0.005

  ibrav=0
  nat={_nat}
  ntyp={_ntyp}

  ecutwfc=80
  ecutrho=400

  edir=3
  emaxpos=0.95
  eopreg=0.1
  eamp={_e_ext}
/

&ELECTRONS
  conv_thr=1.0d-9
  electron_maxstep=200
  mixing_beta=0.3
/

&IONS
  trust_radius_min=1d-5
/

&CELL
  cell_dofree='2Dxy'
/

ATOMIC_SPECIES
{_masses}

ATOMIC_POSITIONS {{crystal}}
{_atomic_positions}

CELL_PARAMETERS {{angstrom}}
{_lattice}

K_POINTS {{automatic}}
8 8 1 0 0 0
""".format(_e_ext=e_ext, _prefix=name, _masses=masses, _nat=nat, _ntyp=ntyp, _atomic_positions=atomic_positions, _lattice=lattice)
            
        # Write text to the input file
        input_file = 'in/{_name}_{_e_ext}.relax.in'.format(_name=name, _e_ext=round(e_ext*1000))
        output_file = 'out/{_name}_{_e_ext}.relax.out'.format(_name=name, _e_ext=round(e_ext*1000))
        with open(input_file, 'w') as f:
            f.write(text)
 
        # Run the input file
        process = subprocess.Popen("mpirun -n {_cpus} /gscratch/cmt/software/qe-6.4.1_cmt_icc/PW/src/pw.x -npool {_nodes} < {_in} > {_out}".format(_cpus=cpus, _nodes=nodes, _in=input_file, _out=output_file), shell=True, stdout=DEVNULL)
        process.wait()
 
        # Get structure and energy data and write to file
        record = False
        converged = False
        lines = []
        energy = 0
        with open(output_file) as f:
            for line in f:
                if '!' in line:
                    energy = [float(s) for s in regex.findall(r'-?\d+\.?\d*', line)][0]
                if 'Begin final coordinates' in line:
                    record = True
                    converged = True
                elif 'End final coordinates' in line:
                    record = False
 
                if record == True:
                    lines.append(line.rstrip())
        if converged == False:
            print("Broke at: " + str(e_ext))
            break

        cell_parameters = '\n'.join([line.replace('\t', ' ') for line in lines[5:8]])
        atom_positions = lines[10:]
        temp_atom_positions = []
        atoms = dict()
        for line in atom_positions:
            line = line.split()
            atoms[line[0]] = atoms.get(line[0], 0) + 1
            line = [line[1], line[2], line[3], line[0]]
            temp_atom_positions.append(' '.join(line))
        atom_positions = '\n'.join(temp_atom_positions)
        atom_types = ' '.join(list(atoms.keys()))
        atom_numbers = ' '.join([str(int) for int in list(atoms.values())])

        text = """Ba1 Ti1 O3
1.0
{_cell_parameters}
{_atom_types}
{_atom_numbers}
direct
{_atom_positions}""".format(_cell_parameters=cell_parameters, _atom_types=atom_types, _atom_numbers=atom_numbers, _atom_positions=atom_positions)
        with open("{_folder}_relaxed/{_name}_{_e_ext}".format(_folder=folder, _name=name, _e_ext=round(e_ext*1000)), 'w') as f:
            f.write(text)

        # Write to temp file
        with open("{_folder}_relaxed/{_name}_temp".format(_folder=folder, _name=name), 'w') as f:
            f.write(text)

        # Post Proccessing
        text = """
&INPUTPP
  prefix={_name}
  outdir='./out'
  filplot='filplot.vpot'
  plot_num=0
/
        """.format(_name=name)
        input_file = 'in/{_name}.pp.in'.format(_name=name)
        output_file = 'out/{_name}.pp.out'.format(_name=name)
        with open(input_file, 'w') as f:
            f.write(text)
        process = subprocess.Popen("/gscratch/cmt/software/qe-6.4.1_cmt_icc/PP/src/pp.x < {_in} > {_out}".format(_in=input_file, _out=output_file), shell=True, stdout=DEVNULL)
        process.wait()
 
        text = """
1
filplot.vpot
1.0D0
300
3
{_awin}
        """.format(_awin=lattice_x.split()[0])
        input_file = 'in/{_name}.average.in'.format(_name=name)
        output_file = 'out/{_name}.average.out'.format(_name=name)
        with open(input_file, 'w') as f:
            f.write(text)
        process = subprocess.Popen("/gscratch/cmt/software/qe-6.4.1_cmt_icc/PP/src/average.x < {_in} > {_out}".format(_in=input_file, _out=output_file), shell=True, stdout=DEVNULL)
        process.wait()
 
        process=subprocess.Popen("mv avg.dat rho_{_name}_{_e_ext}.dat".format(_name=name, _e_ext=round(e_ext*1000)), shell=True, stdout=DEVNULL)
        process.wait()
    

    # Run scf with structure relaxed without electric field

    # Remove data from structure file
    try:
        structure_file = open("{_folder}_relaxed/{_name}".format(_folder=folder, _name=name), 'r')
    except:
        continue
    lines = structure_file.readlines()
    structure_file.close()
    
    # Extract data and process it into useful strings
    lattice_x = lines[2]
    lattice_y = lines[3]
    lattice_z = lines[4]
    
    lattice = ''.join([lattice_x, lattice_y, lattice_z])

    atomic_positions = lines[8:]
    _temp = []
    for i in atomic_positions:
        line = i.split()
        _temp.append(' '.join(line[3:] + line[:3]))
    atomic_positions = '\n'.join(_temp)
    
    typs = lines[5].split()
    number = lines[6].split()
    d_count = dict()

    for i,typ in enumerate(typs):
        if typ in d_count:
            d_count[typ] = str(int(d_count[typ]) + int(number[i]))
        else:
            d_count[typ] = number[i]

    ntyp = len(d_count)
    nat = sum(list(map(int, number)))

    masses_file = open('atomic_weights', 'r')
    lines = masses_file.readlines()
    masses_file.close()
    masses = []
    for line in lines:
        split_line = line.split()
        if split_line[0] in d_count:
            masses.append(split_line[0] + ' ' + split_line[1] + ' ' + split_line[0] + '.pbesol-nc-sr.upf')
    masses = '\n'.join(masses)


    # Write data to input file
    text = """
&CONTROL
  calculation=scf
  prefix={_prefix}
  outdir='./out/'
  pseudo_dir='./pseudo/'

  forc_conv_thr=1.0d-3

  tefield=.TRUE.
  dipfield=.TRUE.
/

&SYSTEM
  occupations='smearing'
  smearing='gaussian'
  degauss=0.005

  ibrav=0
  nat={_nat}
  ntyp={_ntyp}

  ecutwfc=80
  ecutrho=400

  edir=3
  emaxpos=0.95
  eopreg=0.1
  eamp=0
/

&ELECTRONS
  conv_thr=1.0d-9
  electron_maxstep=200
  mixing_beta=0.3
/

&IONS
  trust_radius_min=1d-5
/

&CELL
  cell_dofree='2Dxy'
/

ATOMIC_SPECIES
{_masses}

ATOMIC_POSITIONS {{crystal}}
{_atomic_positions}

CELL_PARAMETERS {{angstrom}}
{_lattice}

K_POINTS {{automatic}}
8 8 1 0 0 0
""".format(_prefix=name, _masses=masses, _nat=nat, _ntyp=ntyp, _atomic_positions=atomic_positions, _lattice=lattice)
            
    # Write text to the input file
    input_file = 'in/{_name}_reference.scf.in'.format(_name=name)
    output_file = 'out/{_name}_reference.scf.out'.format(_name=name)
    with open(input_file, 'w') as f:
        f.write(text)

    # Run the input file
    process = subprocess.Popen("mpirun -n {_cpus} /gscratch/cmt/software/qe-6.4.1_cmt_icc/PW/src/pw.x -npool {_nodes} < {_in} > {_out}".format(_cpus=cpus, _nodes=nodes, _in=input_file, _out=output_file), shell=True, stdout=DEVNULL)
    process.wait()

    # Post Processing
    text = """
&INPUTPP
  prefix={_name}
  outdir='./out'
  filplot='filplot.vpot'
  plot_num=0
/
    """.format(_name=name)
    input_file = 'in/{_name}.pp.in'.format(_name=name)
    output_file = 'out/{_name}.pp.out'.format(_name=name)
    with open(input_file, 'w') as f:
        f.write(text)
    process = subprocess.Popen("/gscratch/cmt/software/qe-6.4.1_cmt_icc/PP/src/pp.x < {_in} > {_out}".format(_in=input_file, _out=output_file), shell=True, stdout=DEVNULL)
    process.wait()

    text = """
1
filplot.vpot
1.0D0
300
3
{_awin}
    """.format(_awin=lattice_x.split()[0])
    input_file = 'in/{_name}.average.in'.format(_name=name)
    output_file = 'out/{_name}.average.out'.format(_name=name)
    with open(input_file, 'w') as f:
        f.write(text)
    process = subprocess.Popen("/gscratch/cmt/software/qe-6.4.1_cmt_icc/PP/src/average.x < {_in} > {_out}".format(_in=input_file, _out=output_file), shell=True, stdout=DEVNULL)
    process.wait()

    process=subprocess.Popen("mv avg.dat avg_{_name}_reference.dat".format(_name=name), shell=True, stdout=DEVNULL)
    process.wait()

