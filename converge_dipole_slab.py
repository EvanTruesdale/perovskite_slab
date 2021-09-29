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
folder = 'structure_files/slab'
structures = join(cwd, folder)
structure_files = [f for f in listdir(structures) if isfile(join(structures, f))]
calculation = 'scf'

for name in structure_files:
    print("NAME: ")
    print(name)

    # Remove data from structure file
    structure_file = open(join(folder, name), 'r')
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
    
    dipole_length_array = np.array([0.2, 0.1, 0.05, 0.025, 0.0125])
    for dipole_length in dipole_length_array:
        emaxpos = 1-dipole_length
        eopreg = dipole_length

        # Write data to input file
        text = """
&CONTROL
  calculation={_calculation}
  prefix={_prefix}
  outdir='./out/'
  pseudo_dir='./pseudo/'

  tefield=.TRUE.
  dipfield=.TRUE.
/

&SYSTEM
  occupations='smearing'
  smearing='gaussian'
  degauss=0.001

  ibrav=0
  nat={_nat}
  ntyp={_ntyp}

  ecutwfc=80

  edir=3
  emaxpos={_emaxpos}
  eopreg={_eopreg}
  eamp=0
/

&ELECTRONS
  conv_thr=1.0d-6
/

&IONS
/

&CELL
/

ATOMIC_SPECIES
{_masses}

ATOMIC_POSITIONS {{crystal}}
{_atomic_positions}

CELL_PARAMETERS {{angstrom}}
{_lattice}

K_POINTS {{automatic}}
8 8 1 0 0 0
""".format(_emaxpos=emaxpos, _eopreg=eopreg, _calculation=calculation, _prefix=name, _masses=masses, _nat=nat, _ntyp=ntyp, _atomic_positions=atomic_positions, _lattice=lattice)
            
        # Write text to the input file
        input_file = 'in/{_name}.{_calculation}.in'.format(_calculation=calculation, _name=name)
        output_file = 'out/{_name}.{_calculation}.out'.format(_calculation=calculation, _name=name)
        with open(input_file, 'w') as f:
            f.write(text)

        # Run the input file
        process = subprocess.Popen("mpirun -n {_cpus} /gscratch/cmt/software/qe-6.4.1_cmt_icc/PW/src/pw.x -npool {_nodes} < {_in} > {_out}".format(_cpus=cpus, _nodes=nodes, _in=input_file, _out=output_file), shell=True, stdout=DEVNULL)
        process.wait()
 
        # Post Proccessing
        text = """
&INPUTPP
  prefix={_name}
  outdir='./out'
  filplot='filplot.vpot'
  plot_num=11
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
100
3
{_awin}
        """.format(_awin=lattice_x.split()[0])
        input_file = 'in/{_name}.average.in'.format(_name=name)
        output_file = 'out/{_name}.average.out'.format(_name=name)
        with open(input_file, 'w') as f:
            f.write(text)
        process = subprocess.Popen("/gscratch/cmt/software/qe-6.4.1_cmt_icc/PP/src/average.x < {_in} > {_out}".format(_in=input_file, _out=output_file), shell=True, stdout=DEVNULL)
        process.wait()
 
        process=subprocess.Popen("mv avg.dat avg_{_dipole_length}.dat".format(_dipole_length=int(dipole_length*100)), shell=True, stdout=DEVNULL)
        process.wait()

    process = subprocess.Popen("tar -cvzf data.tar.gz avg*", shell=True, stdout=DEVNULL)
    process.wait()
