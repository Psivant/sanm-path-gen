import MDAnalysis as mda
from MDAnalysis.analysis import rms
from MDAnalysis.analysis import align
import numpy as np
import matplotlib.pyplot as plt
import os,sys
from sys import exit, stdout, stderr

path_start_pdb = str(sys.argv[1])
old_path_pdb   = str(sys.argv[2])
traj_dcd       = str(sys.argv[3])
new_path_pdb   = str(sys.argv[4])
new_path_pdb2  = new_path_pdb.replace("models", "plumed")

old_path = mda.Universe(path_start_pdb, old_path_pdb)
all_traj = mda.Universe(path_start_pdb, traj_dcd)

n_path_frames = len(old_path.trajectory)
n_traj_frames = len(all_traj.trajectory)
prmsd = np.zeros( (2,               # z-axis
                   n_path_frames,   # y-axis
                   n_traj_frames) ) # x-axis
new_path_indx = np.zeros( len(old_path.trajectory), dtype=np.int64 )
new_path_rmsd = np.zeros( len(old_path.trajectory) )

for i, frame_path in enumerate(old_path.trajectory):

    r = rms.RMSD(all_traj, old_path, select='name CA', ref_frame=i).run()
    rrr = r.results.rmsd
    rrr = rrr[rrr[:, -1].argsort()] # sort all columns by the RMSD
    prmsd[0][i] = rrr[:,  0]  # 1st column traj index
    prmsd[1][i] = rrr[:, -1]  # 3rd column RMSD value

#for i in range(n_path_frames):
    #print(i, int(prmsd[0][i][0]), prmsd[1][i][0])
rmsd_cutoff = np.amax(prmsd[1,:,0]) + 0.1
print ("RMSD cutoff= ", rmsd_cutoff)

# plt.imshow(prmsd[1], cmap='viridis', aspect='auto')
# plt.xlabel('Frame (all_traj)')
# plt.ylabel('Frame (old_path)')
# plt.colorbar(label=r'RMSD ($\AA$)');
# plt.show()

rmsd_min = 1e10
for i in range(n_traj_frames):

    if( prmsd[1][0][i] > rmsd_cutoff ):
        break

    for j in range(n_traj_frames):
    
        if( prmsd[1][1][j] > rmsd_cutoff ):
            break

        mobile_indx = int(prmsd[0][1][j])
        target_indx = int(prmsd[0][0][i])
        if( mobile_indx == target_indx ):
            continue
        atom_selection = all_traj.select_atoms('name CA')
        all_traj.trajectory[mobile_indx]
        mobile_coords = atom_selection.positions.copy()
        all_traj.trajectory[target_indx]
        target_coords = atom_selection.positions.copy()     
        rmsd = rms.rmsd( mobile_coords, target_coords, superposition=True )
        if( rmsd < rmsd_min ):
            rmsd_min = rmsd
            new_path_indx[0] = target_indx
            new_path_indx[1] = mobile_indx
            new_path_rmsd[1] = rmsd_min
            #print( rmsd_min, target_indx, mobile_indx );
#print("--------------------")

for i in range(2, n_path_frames):

    rmsd_min = 1e10
    for j in range(n_traj_frames):
    
        if( prmsd[1][i][j] > rmsd_cutoff ):
            break

        mobile_indx = int(prmsd[0][i][j])
        skip = False
        for k in range(i):
            if( mobile_indx == new_path_indx[k] ):
                skip = True
                break
        if( skip ):
            continue
        target_indx = new_path_indx[i-1]
        if( mobile_indx == target_indx ):
            continue
        atom_selection = all_traj.select_atoms('name CA')
        all_traj.trajectory[mobile_indx]
        mobile_coords = atom_selection.positions.copy()
        all_traj.trajectory[target_indx]
        target_coords = atom_selection.positions.copy()     
        rmsd = rms.rmsd( mobile_coords, target_coords, superposition=True )
        if( rmsd < rmsd_min ):
            rmsd_min = rmsd
            new_path_indx[i] = mobile_indx
            new_path_rmsd[i] = rmsd_min
            #print( rmsd_min, target_indx, mobile_indx );
    #print("--------------------",i)

rmsd_avg = np.mean( new_path_rmsd[1:] )
lam4plumed = 230. / rmsd_avg**2

atom_selection = old_path.select_atoms('name CA')
old_path.trajectory[0]
refA_coords = atom_selection.positions.copy()
old_path.trajectory[n_path_frames -1]
refB_coords = atom_selection.positions.copy()

print( 'traj index      RMSD    RMSD-A    RMSD-B' )
print( '----------------------------------------' )
for i in range(n_path_frames):

    atom_selection = all_traj.select_atoms('name CA')
    all_traj.trajectory[new_path_indx[i]]
    frame_coords = atom_selection.positions.copy()
    rmsdA = rms.rmsd( frame_coords, refA_coords, superposition=True )
    rmsdB = rms.rmsd( frame_coords, refB_coords, superposition=True )
    print( '{0:10d}  {1:8.3f}  {2:8.3f}  {3:8.3f}'.format( new_path_indx[i], new_path_rmsd[i], rmsdA, rmsdB))
print( '----------------------------------------' )
print( ' Suggested LAMBDA for PLUMED PATHMSD = {0:d}'.format( int(lam4plumed)))


with mda.Writer('tmp.pdb', multiframe=True) as pdb:
    for i in range(n_path_frames):

        atom_selection = all_traj.select_atoms('protein')
        all_traj.trajectory[new_path_indx[i]]
        pdb.write(all_traj)

trj = mda.Universe('tmp.pdb')
ref = mda.Universe(path_start_pdb)

align.AlignTraj(trj,  # trajectory to align
                ref,  # reference
                select='name CA',       # selection of atoms to align
                filename=new_path_pdb,  # file to write the trajectory to
                match_atoms=True,       # whether to match atoms based on mass
               ).run()

if os.path.isfile('tmp.pdb'):
    os.remove('tmp.pdb')

# PLUMEDify
command_line = "gawk '{if($1==\"ATOM\" || $1==\"ENDMDL\") print $0}' "+new_path_pdb+" > tmp.pdb"
print( command_line )
os.system( command_line )
command_line = "gawk -i inplace '{if($1==\"ATOM\") printf \"%s%s%s\\n\", substr($0,1,54), \"  1.00  1.00\", substr($0,67,length($0)-66); else print $0;}' tmp.pdb"
print( command_line )
os.system( command_line )
command_line = "sed -i 's/ENDMDL/END/g' tmp.pdb"
print( command_line )
os.system( command_line )
command_line = "mv tmp.pdb "+new_path_pdb2
print( command_line )
os.system( command_line )
