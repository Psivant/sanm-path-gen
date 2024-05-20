# sanm-path-gen
# Smooth anisotropic network model (SANM) protein conformational transition path generator 

## Released under MIT License

## Copyright (c) 2023 Psivant Therapeutics, LLC.

## Copyright (c) 2023, 2024 Istvan B Kolossvary.

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

### Please cite this publication in any work where you use our software or its derivatives.
Istvan B Kolossvary and Woody Sherman, Comprehensive Approach to Simulating Large Scale Conformational Changes in Biological Systems Utilizing a Path Collective Variable and New Barrier Restraint, https://pubs.acs.org/doi/10.1021/acs.jpcb.3c02028

### This README file provides basic documentation through a complex example. 
Following the instructions below will regenerate the test files and familiarize the user with the use of the software. One needs two PDB files representing two different conformations of a protein (or protein complex), in the example these files are STING-metin-AA.pdb and STING-metout-AA.pdb, respectively. The files should contain a full protein system. Steps 2-4 will build a minimum energy path "morphing" one endpoint structure into the other. Note that the computations in Step 4 require AmberTools. There will be three types of outputs, a CA path, a backbone plus CB (BBCB) path, and an all-atom (AA) path. Moreover, each output file has two versions, "models" that is a standard multi-PDB file for visualization, and "plumed" that is formatted to be used with the PLUMED software.

### 1. Compile the transition path generator and backbone builder programs:
```
gcc -o CA-transition-path.exe -O3 CA-transition-path.c common.c -lm
./CA-transition-path.exe
CA-transition-path.c(2551): ERROR
	Usage:  ./CA-transition-path.exe  n_CA_atoms  endpoints_confs_pdb_fname  CA_transition_path_fname(no extension)  rmsd_gap_between_output_frames_Angs

gcc -o ca2bb.exe -O3 ca2bb.c -lm
./ca2bb.exe
Usage:  ./ca2bb.exe  CA-pdb-fname  backbone-pdb-fname  chain_gap
```
### 2. Concatenate the endpoint structures:
```
cat STING-metin-AA.pdb STING-metout-AA.pdb > STING-endpoints-AA.pdb
grep -c CA STING-metin-AA.pdb
368
```
### 3. Generate the CA transtion path:
```
./CA-transition-path.exe 368 STING-endpoints-AA.pdb STING-metinout-path-CA 0.3 >& STING-metinout-path-CA.log
```
The logfile shows RMSD values between consecutive frames and also with respect to the endpoint structures. Note that there are two paths generated one termed "forward" and another "reverse", see details in the publication. The logfile also inlcudes a template of PLUMED input for running path-CV meta-eABF simulation using the generated path.

### 4. Generate the backbone and all-atom transition-paths: (This bash script requires AmberTools v2018 or later.)
```
./Run_AA-transition-path.sh
Usage:  ./Run_AA-transition-path.sh CA-transition-path-pdb-fname(in)  BBCB-transition-path-pdb-fname(out)  AA-transition-path-pdb-fname(out)  AA-frame#first-pdb-fname(in)  AA-frame#last-pdb-fname(in)  AA-Amber-parm-file(in)  num_cpu
^C
./Run_AA-transition-path.sh STING-metinout-path-CA-forward-models.pdb STING-metinout-path-BBCB-forward-models.pdb STING-metinout-path-AA-forward-models.pdb STING-metin-AA.pdb STING-metout-AA.pdb STING-AA.parm7 16 >& STING-metinout-path-AA-forward.log
```
### 5. Prune the all-atom path to only include (besides CA) the Met-loop MET SD atoms:
```
awk '{if($1!="ATOM" || $3=="CA" || ($3=="SD" && $6==61) || ($3=="SD" && $6==245)) print $0}' STING-metinout-path-AA-forward-plumed.pdb > STING-metinout-path-CA+MET-SD-forward-plumed.pdb
```
### 6. Reconstruct the path using trajectory frames from an existing simulation:
We provide a Python script that---given an all-atom (and preferably explicit-solvent) simulation trajectory based on any of the paths above---will construct a new path by slecting frames from the trajectory such that the trajectory frames are < 0.5 RMSD from their respective path nodes, and consecutive new frames are as close to each other as possible.
```
python find-path-in-traj.py STING-metin-AA.pdb STING-metinout-path-AA-forward-models.pdb STING-traj.dcd STING-metinout-traj-based-path-AA-forward-models.pdb 0.5 >& find-path-in-traj.log
```
