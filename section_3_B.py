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

    # Run scf on with electric field without field-induced relaxations
    e_ext_array = np.arange(0.005, 0.085, 0.005)
    for e_ext in e_ext_array:

        # Remove data from structure file
        try:
            structure_file = open("{_folder}_relaxed/{_name}".format(_folder=folder, _name=name), 'r')
        except:
            print("No relaxed structure for {_name}".format(_name=name))
            break
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
  emaxpos=0.9875
  eopreg=0.0125
  eamp={_e_ext}
/

&ELECTRONS
  conv_thr=1.0d-9
  electron_maxstep=200
  mixing_beta=0.3
/

&IONS
  trust_radius_min=1d-4
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
        input_file = 'in/{_name}_{_e_ext}.scf.in'.format(_name=name, _e_ext=round(e_ext*1000))
        output_file = 'out/{_name}_{_e_ext}.scf.out'.format(_name=name, _e_ext=round(e_ext*1000))
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
 
        process=subprocess.Popen("mv avg.dat avg_rho_{_e_ext}.dat".format(_e_ext=round(e_ext*1000)), shell=True, stdout=DEVNULL)
        process.wait()

