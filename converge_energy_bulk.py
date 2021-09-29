import os
from os import listdir
from os.path import isfile, join
import numpy as np
import regex
import subprocess
import sys

# Get input for CPU usage
cpus = sys.argv[1]
nodes = sys.argv[2]

# Get all structure files
cwd = os.getcwd()
folder = 'structure_files/bulk'
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

    # Generate ecutwfc, ecutrho values
    row = np.arange(30, 90, 5)
    ecutwfc_array = np.stack((row, row, row), axis=0)
    ecutrho_array = np.copy(ecutwfc_array)
    ecutrho_array[0,:] *= 4
    ecutrho_array[1,:] *= 5
    ecutrho_array[2,:] *= 6
    energies = np.zeros(ecutwfc_array.shape)    

    # Run through all values of ecutwfc, ecutrho
    for i, row in enumerate(ecutwfc_array):
        for j, col in enumerate(row):
            ecutwfc = ecutwfc_array[i,j]
            ecutrho = ecutrho_array[i,j]
            
            # Write data to input file
            text = """
&CONTROL
  calculation={_calculation}
  prefix={_prefix}
  outdir='./out/'
  pseudo_dir='./pseudo/'
/

&SYSTEM
  occupations='smearing'
  smearing='gaussian'
  degauss=0.001

  ibrav=0
  nat={_nat}
  ntyp={_ntyp}

  ecutwfc={_ecutwfc}
  ecutrho={_ecutrho}
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
4 4 4 0 0 0
""".format(_calculation=calculation, _ecutwfc=ecutwfc, _ecutrho=ecutrho, _prefix=name, _masses=masses, _nat=nat, _ntyp=ntyp, _atomic_positions=atomic_positions, _lattice=lattice)
            
            # Write text to the input file
            input_file = 'in/{_name}.{_calculation}.in'.format(_calculation=calculation, _name=name)
            output_file = 'out/{_name}.{_calculation}.out'.format(_calculation=calculation, _name=name)
            with open(input_file, 'w') as f:
                f.write(text)

            # Run the input file
            process = subprocess.Popen("mpirun -n {_cpus} /gscratch/cmt/software/qe-6.4.1_cmt_icc/PW/src/pw.x -npool {_nodes} < {_in} > {_out}".format(_cpus=cpus, _nodes=nodes, _in=input_file, _out=output_file), shell=True, stdout=subprocess.PIPE)
            process.wait()

            # Get the output data
            energy_steps = np.array([])
            with open(output_file) as f:
                for line in f:
                    if '!' in line:
                        energy_steps = np.append(energy_steps, [float(s) for s in regex.findall(r'-?\d+\.?\d*', line)])
            if len(energy_steps) != 0:
                energy = energy_steps[len(energy_steps)-1]
            else:
                energy = 0
            energies[i,j] = energy
            print("{wfc} {rho} {e}".format(wfc=ecutwfc, rho=ecutrho, e=energies[i,j]))
            with open('data/energy', 'a') as f:
                f.write("{wfc} {rho} {e}\n".format(wfc=ecutwfc, rho=ecutrho, e=energies[i,j]))

