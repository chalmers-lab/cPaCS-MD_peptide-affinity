# This Python scripts executes PaCS MD to smash in small ligands to GPCRs

#run using conda environment containing MDTraj and Numpy
import mdtraj as md
import numpy as np
import os
import time
import math
import heapq
from ast import literal_eval

#define PaCS parameters here
total_batches = 500 #insert total no of batch here
batch_size = 3 #insert no of runs per cycle here
cycle_size = 3 #insert no of cycles per batch here
max_wait = 1000
slurm_mode = False

#load initial file
init = 'init.gro'

#make a pdb
#os.popen('gmx trjconv -f init.gro -s init.gro -o init.pdb << EOF \n 0 \n EOF')
#time.sleep(2)
topology = md.load(init, top='init.gro').topology

#Define ligand and Contacts
ligand = topology.select('protein and residue < 20 and mass > 2')
ligand_BB = topology.select('protein and residue < 20 and backbone')
receptor = topology.select('protein and residue > 20 and mass > 2')
receptor_BB = topology.select('protein and residue > 20 and backbone')
protlig = topology.select('protein')
print(ligand)
atom_index_selection = md.compute_neighbors(md.load(init), cutoff=0.5, query_indices=ligand_BB, haystack_indices=receptor)[0]
print('\n \n \n \n \n ')
print(atom_index_selection)
best_atom = 1
#compute pairs
atom_pairs = []



if os.path.exists('./contacts'):
  c_file = open('contacts', "r")
  c_lines = c_file.readlines()
  print(str(list(ligand)))
  print(str(c_lines[0]))
  if not c_lines:
    os.popen('rm contacts').read()
    for i in ligand_BB:
      print('begining search for atom #  ' + str (i))
      deviation = 999
      
      lig_atom_com = md.compute_center_of_mass(md.load(init, top='init.gro').atom_slice([i]))
      for atoms in atom_index_selection:
        receptor_atom_com = md.compute_center_of_mass(md.load(init, top='init.gro').atom_slice([atoms]))
        best_atom
        #print('lig com is: ' + str(lig_atom_com[0]))
        #print('rec com is: ' + str(receptor_atom_com[0]))
        if np.sqrt(np.sum(np.square(lig_atom_com[0]-receptor_atom_com[0]))) < deviation:
          tentative_rec_atom = atoms
          deviation = np.sqrt(np.sum(np.square(lig_atom_com[0]-receptor_atom_com[0])))
          best_atom = atoms
          print(atoms)
          print(deviation)
        
      atom_pairs.append([i,best_atom,deviation])
    print(ligand)
    print(atom_pairs)
    os.popen('touch contacts')
    c_file = open('contacts', "w")
    c_file.writelines([str(list(ligand))+ '\n'])
    c_file.writelines(str(list(atom_pairs)))
    
  elif str(c_lines[0].strip()) == str(list(ligand)).strip():
    atom_pairs = literal_eval(c_lines[1])
    print('sucess!!')
    print(atom_pairs)
    print(type(atom_pairs))
  
  else:
    os.popen('rm contacts').read()
    for i in ligand_BB:
      print('begining search for atom #  ' + str (i))
      deviation = 999
        
      lig_atom_com = md.compute_center_of_mass(md.load(init, top='init.gro').atom_slice([i]))
      for atoms in atom_index_selection:
        receptor_atom_com = md.compute_center_of_mass(md.load(init, top='init.gro').atom_slice([atoms]))
        best_atom
        #print('lig com is: ' + str(lig_atom_com[0]))
        #print('rec com is: ' + str(receptor_atom_com[0]))
        if np.sqrt(np.sum(np.square(lig_atom_com[0]-receptor_atom_com[0]))) < deviation:
          tentative_rec_atom = atoms
          deviation = np.sqrt(np.sum(np.square(lig_atom_com[0]-receptor_atom_com[0])))
          best_atom = atoms
          print(atoms)
          print(deviation)
          
      atom_pairs.append([i,best_atom,deviation])
    print(ligand)
    print(atom_pairs)
    os.popen('touch contacts')
    c_file = open('contacts', "w")
    c_file.writelines([str(list(ligand)) + '\n'])
    c_file.writelines(str(list(atom_pairs)))
      
else:
    print('contact file doesnt exist')
    for i in ligand_BB:
      print('begining search for atom #  ' + str (i))
      deviation = 999
        
      lig_atom_com = md.compute_center_of_mass(md.load(init, top='init.gro').atom_slice([i]))
      for atoms in atom_index_selection:
        receptor_atom_com = md.compute_center_of_mass(md.load(init, top='init.gro').atom_slice([atoms]))
        best_atom
        #print('lig com is: ' + str(lig_atom_com[0]))
        #print('rec com is: ' + str(receptor_atom_com[0]))
        if np.sqrt(np.sum(np.square(lig_atom_com[0]-receptor_atom_com[0]))) < deviation:
          tentative_rec_atom = atoms
          deviation = np.sqrt(np.sum(np.square(lig_atom_com[0]-receptor_atom_com[0])))
          best_atom = atoms
          print(atoms)
          print(deviation)
          
      atom_pairs.append([i,best_atom,deviation])
    print(ligand)
    print(atom_pairs)
    os.popen('touch contacts')
    c_file = open('contacts', "w")
    c_file.writelines([str(list(ligand))+ '\n'])
    c_file.writelines(str(list(atom_pairs)))

  

print(np.array(atom_pairs))

#input("Press Enter to continue...")

print(ligand)
print(receptor)


if os.path.exists('./test_index.ndx'):
    print('index exists')
else:
    os.popen('cp index.ndx test_index.ndx').read()
    os.popen('echo \"\" >> test_index.ndx').read()
    os.popen('echo \"[center]\" >> test_index.ndx').read()
    os.popen('echo ' + str(ligand[0]) + ' >> test_index.ndx').read()

#define different measurements
def center_of_mass(trj_name):
    load = md.load(trj_name, top=topology).center_coordinates().image_molecules()
    atompairs = np.array(atom_pairs)
    dist_array = []
    for entry in atompairs:
      #print(entry[0])
      lig_load = load.atom_slice([int(entry[0])])
      lig_ctr = md.compute_center_of_mass(lig_load)
      rec_load = load.atom_slice([int(entry[1])]) #modified to include only close contacts
      rec_ctr = md.compute_center_of_mass(rec_load)
      dist = np.sqrt(np.sum(np.square(lig_ctr-rec_ctr)))-entry[2]
      print(entry[0], entry[1], entry[2], dist)
      dist_array.append(dist)
    dist_avg = np.mean(dist_array)
    return dist_avg

#define detecting if things had finished.
def batch_detect_slurm(batch_name, cycle_size, batch_size):
    for x in range(cycle_size):        
        for i in range(batch_size):
            head_name ="{}_{}_{}".format(batch_name, x, i) 
            run_name = "{}.gro".format(head_name)
            print('detecting', run_name)
            count = 0
            while not os.path.exists(run_name):
                time.sleep(1)
                count = count+1
                if os.path.isfile(run_name):
                    break
                elif count == max_wait:
                    print('Woopsie it died mate so sorry!')
                    exit()

def batch_detect_loc(batch_name, cycle_no, batch_no):
    head_name ="{}_{}_{}".format(batch_name, cycle_no, batch_no) 
    run_name = "{}.gro".format(head_name)
    print('detecting', run_name)
    count = 0
    while not os.path.exists(run_name):
        time.sleep(1)
        count = count+1
        if os.path.isfile(run_name):
            break
        elif count == max_wait:
            print('Woopsie it died mate so sorry!')
            exit()

#define running a batch of simulations in a slurm cluster
def batch_run_slurm(batch_name, cycle_no, init_gro, batch_size):
    print(batch_name)
    print(cycle_no)
    print(init_gro)
    print(batch_size)
    if os.path.exists('grompp_file.tpr'):
        print('grompp found')
    else:
        os.system("/apps/linux/gromacs/gromacs-2023.1_Intel_CUDA11.8/bin/gmx grompp -f prod_100ps.mdp -p topol.top -n index.ndx -maxwarn 4 -c init.gro -o grompp_file.tpr -quiet")
    for i in range(batch_size):
        grompp_command = "/apps/linux/gromacs/gromacs-2023.1_Intel_CUDA11.8/bin/gmx grompp -f prod_100ps.mdp -p topol.top -n index.ndx -maxwarn 4 -c {} -o {}_{}_{}.tpr -quiet".format(init_gro, batch_name, cycle_no, i)
        os.system(grompp_command)
        print(grompp_command)
    for i in range(batch_size):
        slurm_name = "{}_{}_{}.slurm".format(batch_name, cycle_no, i)
        cmd_run = "/apps/linux/gromacs/gromacs-2023.1_Intel_CUDA11.8/bin/gmx mdrun -deffnm {}_{}_{} -v -notunepme".format(batch_name, cycle_no, i)
        with open(slurm_name, 'w') as f:
            f.writelines('#!/bin/bash')
            f.write('\n') 
            f.writelines('#SBATCH --partition=CLUSTER')
            f.write('\n')
            f.writelines('#SBATCH --nodes=1')
            f.write('\n')
            f.writelines('#SBATCH --ntasks=1')
            f.write('\n')
            f.writelines('#SBATCH --cpus-per-task=24')
            f.write('\n')
            f.writelines('#SBATCH --gres=gpu:1')
            f.write('\n')
            f.writelines('#SBATCH --mem=0')
            f.write('\n')
            f.writelines('#SBATCH --time=00-00:05:00') #add approximate time to run an MD here, will need adjusting
            f.write('\n')
            f.writelines(cmd_run)
            f.write('\n')


def batch_run_loc(batch_name, cycle_no, init_gro, batch_size):
    print(batch_name)
    print(cycle_no)
    print(init_gro)
    print(batch_size)
    if os.path.exists('grompp_file.tpr'):
        print('grompp found')
    else:
        os.system("/apps/linux/gromacs/gromacs-2023.1_Intel_CUDA11.8/bin/gmx grompp -f prod_100ps.mdp -p topol.top -n index.ndx -maxwarn 4 -c init.gro -o grompp_file.tpr -quiet")
    for i in range(batch_size):
        
        grompp_command = "/apps/linux/gromacs/gromacs-2023.1_Intel_CUDA11.8/bin/gmx grompp -f prod_100ps.mdp -p topol.top -n index.ndx -maxwarn 2 -c {} -o {}_{}_{}.tpr -quiet".format(init_gro, batch_name, cycle_no, i)
        os.system(grompp_command)
        print(grompp_command)
    for i in range(batch_size):
        cmd_run = "/apps/linux/gromacs/gromacs-2023.1_Intel_CUDA11.8/bin/gmx mdrun -deffnm {}_{}_{} -ntomp 12 -cpi -v -gpu_id 0".format(batch_name, cycle_no, i)
        os.system(cmd_run)
        batch_detect_loc(batch_name=batch_name, cycle_no=cycle_no, batch_no=i)
        trjconv_command = "/apps/linux/gromacs/gromacs-2023.1_Intel_CUDA11.8/bin/gmx trjconv -f {}_{}_{}.xtc -s grompp_file.tpr -o {}_{}_{}.gro -dump 100 -pbc mol -center -n test_index.ndx << EOF \n 3 \n 2 \n".format(batch_name, cycle_no, i, batch_name, cycle_no, i)
        os.system(trjconv_command)



#test measure function


dist_rms = center_of_mass(init)
dist_pen = dist_rms
print(dist_pen)

#define analyzing a batch of simulations
def batch_analyze(batch_name, cycle_size, batch_size):
    batch_trajs = []
    batch_rmsds = []
    for x in range(cycle_size):
        for i in range(batch_size):
            trj_name = "{}_{}_{}.gro".format(batch_name, x, i)
            batch_trajs.append(trj_name)
            dist_rms = center_of_mass(trj_name)
            dist_pen = dist_rms
            batch_rmsds.append(dist_pen) 
    min_rmsd = heapq.nlargest(cycle_size, batch_rmsds)
#    min_rmsd = heapq.nsmallest(cycle_size, batch_rmsds)
    print(min_rmsd)
    output_log = "{}_out.txt".format(batch_name)
    with open(output_log, 'w') as f:
        for name in batch_trajs:
                f.writelines(name)
                f.writelines(' ')
                f.writelines(str(batch_rmsds[batch_trajs.index(name)]))
                f.write('\n')
    output_measures = "{}_measures.txt".format(batch_name)
    with open(output_measures, 'w') as f:
        for name in batch_trajs:
                f.writelines(str(batch_rmsds[batch_trajs.index(name)]))
                f.write('\n')
    output_names = "{}_names.txt".format(batch_name)
    with open(output_names, 'w') as f:
        for name in batch_trajs:
                f.writelines(name)
                f.write('\n')
    min_no = []
    for i in min_rmsd:
        n = batch_rmsds.index(i)
        min_no.append(n)
    min_traj = []
    for i in min_no:
        name = batch_trajs[i]
        min_traj.append(name)
    out_txt = "{}_tops.txt".format(batch_name)
    with open(out_txt, 'w') as f:
        for name in min_traj:
                f.writelines(name)
                f.write('\n')

#run PaCS MD
batch_tops = []
batch_list = []
prog_rmsds = []
for n in range(total_batches):
    
    batch_name = "b{0:0=4d}".format(n)
    top_name ="{}_tops.txt".format(batch_name)
    if os.path.isfile(top_name) or os.path.isfile(batch_name + '.tar.bz2'):
        pass
    else:
        start = n
        rm_cmd = "rm b{0:0=4d}*".format(n)
        os.system(rm_cmd)
        break
for n in range(start, total_batches, 1):
    if n == 0:
        batch_name = "b0000"
        if slurm_mode==True:
            for tops in last_tops:
                batch_run_slurm(batch_name=batch_name, cycle_no=last_tops.index(tops), init_gro=tops.strip(), batch_size=batch_size)
            sbatch_cmd = "ls b{0:0=4d}_*.slurm | xargs -n 1 sbatch".format(n)
            print(sbatch_cmd)
            os.system(sbatch_cmd)
            batch_detect_slurm(batch_name=batch_name, cycle_size=cycle_size, batch_size=batch_size)
            print("Finish detecting......")
        else:
            for tops in last_tops:
                batch_run_loc(batch_name=batch_name, cycle_no=last_tops.index(tops), init_gro=tops.strip(), batch_size=batch_size)
            print("Finish detecting.....")
        batch_analyze(batch_name=batch_name, cycle_size=cycle_size, batch_size=batch_size)
        if n == 1:
            pass
        else:
            prev_batch_del = "rm b{0:0=4d}*.tpr b{0:0=4d}*.cpt b{0:0=4d}*.slurm b{0:0=4d}*.edr *out *#".format(n-2, n-2, n-2,n-2)
            prev_batch_tar = "tar cfvj b{0:0=4d}.tar.bz2 --remove-files b{0:0=4d}* & ".format(n-2, n-2, n-2)
            os.system(prev_batch_del)
            os.system(prev_batch_tar)

        top_values = open(f'{batch_name}_measures.txt', 'r').readlines()
        top_floats = np.array([float(x) for x in top_values])
        mean_top = np.mean(top_floats)
        print(f'for cycle {n} the distance is {mean_top}')
        if mean_top > 8.5:
            print('distance exceeded, breaking')
            exit()

print('finished PaCS MD')

