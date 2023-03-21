#!/bin/bash

# Released under MIT License

# Copyright (c) 2023 Psivant Therapeutics, LLC.

# Copyright (c) 2023 Istvan B Kolossvary.

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

echo "Usage:  ./Run_AA-transition-path.sh CA-transition-path-pdb-fname(in)  BBCB-transition-path-pdb-fname(out)  AA-transition-path-pdb-fname(out)  AA-frame#first-pdb-fname(in)  AA-frame#last-pdb-fname(in)  AA-Amber-parm-file(in)  num_cpu"

# Variables
CA_path_pdb=$1
BBCB_path_pdb=$2
AA_path_pdb=$3
AA_path_pdb_relax=$(echo $3 | cut -d'.' -f 1)'_relax.pdb'
AA1_pdb=$4
AA2_pdb=$5
AA_parm7=$6
n_cpu=$7

# Add the backbone atoms  (we assume that CA-CA distance > 4.5 Angs represents a chain gap)
./ca2bb.exe $CA_path_pdb $BBCB_path_pdb 4.5
n_frame=$(grep -c MODEL $BBCB_path_pdb)
n_atom=$(grep -c 'ATOM' $BBCB_path_pdb)
((n_atom=n_atom/n_frame))

# Fix the atom numbering in $BBCB_path_pdb to correspond to that of $AA1_pdb
rm -f tmp1.pdb
touch tmp1.pdb
for (( i=0; i<$n_frame; i++ ))
do
    start=$((i*(n_atom+2)+1))
    end=$((start+n_atom+1))

    gawk -v v1="$start" -v v2="$end" 'NR==FNR{if(FNR>=v1 && FNR<=v2)if($3=="CA" || $3=="CB" || $3=="N" || $3=="C" || $3=="O" || $3=="OXT")a[substr($0,13,14)]=substr($0,31,24);next;}{if (a[substr($0,13,14)]) printf("%s%s%s\n", substr($0,1,30), a[substr($0,13,14)], substr($0,55,length($0)-54)); else print $0;}' $BBCB_path_pdb $AA1_pdb > tmp2.pdb

    echo "MODEL"  >> tmp1.pdb
    cat tmp2.pdb  >> tmp1.pdb
    echo "ENDMDL" >> tmp1.pdb

done
# Save only the backbone
gawk '{if($1=="MODEL" || $1=="ENDMDL" || $3=="CA" || $3=="CB" || $3=="N" || $3=="C" || $3=="O" || $3=="OXT") print $0;}' tmp1.pdb > $BBCB_path_pdb

# Change the format to PLUMED
gawk '{if($1=="ATOM" || $1=="HETATM")printf "%s%s%s\n", substr($0,1,54), "  1.00  1.00", substr($0,67,length($0)-66); else print $0;}' $BBCB_path_pdb > tmp1.pdb
gawk '{if($1!="MODEL")print $0;}' tmp1.pdb > tmp2.pdb
sed -i 's/ENDMDL/END/g' tmp2.pdb
mv tmp2.pdb ${BBCB_path_pdb/models/plumed}

# Build the all-atom path:
# Minimize the all-atom left endpoint structure and initialize the path
cpptraj -p $AA_parm7 -y $AA1_pdb -x tmp1.rst7
./Run_xmin.sh tmp1.rst7 tmp1.rst7 $AA_parm7 tmp3.rst7 '@CA,CB' $n_cpu  # restrain CA,CB
cpptraj -p $AA_parm7 -y tmp3.rst7 -x frame1.pdb
mv tmp3.rst7 tmp1.rst7
./Run_xmin.sh tmp1.rst7 tmp1.rst7 $AA_parm7 tmp3.rst7 '!@*' $n_cpu     # unrestrained
cpptraj -p $AA_parm7 -y tmp3.rst7 -x frame1_relax.pdb

# Add up to half of the frames in ascending order
half=$((n_frame/2))
for (( i=1; i<$half; i++ ))
do
    start=$((i*(n_atom+2)+1))
    end=$((start+n_atom+1))

    # Replace (iteratively from left to right) the CA and CB atoms of the all-atom structure with those of the consecutive frames of the BBCB path and then minimize keeping the CA and CB atoms in place
    gawk -v v1="$start" -v v2="$end" 'NR==FNR{if(FNR>=v1 && FNR<=v2)if($3=="CA" || $3=="CB")a[substr($0,1,30)]=substr($0,31,24);next;}{if (a[substr($0,1,30)]) printf("%s%s%s\n", substr($0,1,30), a[substr($0,1,30)], substr($0,55,length($0)-54)); else print $0;}' $BBCB_path_pdb frame${i}.pdb > tmp2.pdb

    cpptraj -p $AA_parm7 -y tmp2.pdb -x tmp2.rst7
    cpptraj -p $AA_parm7 -y frame${i}.pdb -x tmp1.rst7
    ./Run_xmin.sh tmp1.rst7 tmp2.rst7 $AA_parm7 tmp3.rst7 '@CA,CB' $n_cpu  # restrain CA,CB
    cpptraj -p $AA_parm7 -y tmp3.rst7 -x frame$((i+1)).pdb
    mv tmp3.rst7 tmp1.rst7
    ./Run_xmin.sh tmp1.rst7 tmp1.rst7 $AA_parm7 tmp3.rst7 '!@*' $n_cpu     # unrestrained
    cpptraj -p $AA_parm7 -y tmp3.rst7 -x frame$((i+1))_relax.pdb

done

# Minimize the all-atom right endpoint structure
cpptraj -p $AA_parm7 -y $AA2_pdb -x tmp1.rst7
./Run_xmin.sh tmp1.rst7 tmp1.rst7 $AA_parm7 tmp3.rst7 '@CA,CB' $n_cpu  # restrain CA,CB
cpptraj -p $AA_parm7 -y tmp3.rst7 -x frame${n_frame}.pdb
mv tmp3.rst7 tmp1.rst7
./Run_xmin.sh tmp1.rst7 tmp1.rst7 $AA_parm7 tmp3.rst7 '!@*' $n_cpu     # unrestrained
cpptraj -p $AA_parm7 -y tmp3.rst7 -x frame${n_frame}_relax.pdb

# Add the remainder of the frames in descending order
for (( i=$((n_frame-2)); i>=$half; i-- ))
do
    start=$((i*(n_atom+2)+1))
    end=$((start+n_atom+1))

    # Replace (iteratively from right to left) the CA and CB atoms of the all-atom structure with those of the consecutive frames of the BBCB path and then minimize keeping the CA and CB atoms in place
    gawk -v v1="$start" -v v2="$end" 'NR==FNR{if(FNR>=v1 && FNR<=v2)if($3=="CA" || $3=="CB")a[substr($0,1,30)]=substr($0,31,24);next;}{if (a[substr($0,1,30)]) printf("%s%s%s\n", substr($0,1,30), a[substr($0,1,30)], substr($0,55,length($0)-54)); else print $0;}' $BBCB_path_pdb frame$((i+2)).pdb > tmp2.pdb

    cpptraj -p $AA_parm7 -y tmp2.pdb -x tmp2.rst7
    cpptraj -p $AA_parm7 -y frame$((i+2)).pdb -x tmp1.rst7
    ./Run_xmin.sh tmp1.rst7 tmp2.rst7 $AA_parm7 tmp3.rst7 '@CA,CB' $n_cpu  # restrain CA,CB
    cpptraj -p $AA_parm7 -y tmp3.rst7 -x frame$((i+1)).pdb
    mv tmp3.rst7 tmp1.rst7
    ./Run_xmin.sh tmp1.rst7 tmp1.rst7 $AA_parm7 tmp3.rst7 '!@*' $n_cpu     # unrestrained
    cpptraj -p $AA_parm7 -y tmp3.rst7 -x frame$((i+1))_relax.pdb

done

# Concatenate frames into single all-atom path
rm -f $AA_path_pdb
touch $AA_path_pdb
rm -f $AA_path_pdb_relax
touch $AA_path_pdb_relax
for (( i=1; i<=$n_frame; i++ ))
do
    echo "MODEL"    >> $AA_path_pdb
    cat frame${i}.pdb >> $AA_path_pdb
    echo "ENDMDL"   >> $AA_path_pdb
    echo "MODEL"    >> $AA_path_pdb_relax
    cat frame${i}_relax.pdb >> $AA_path_pdb_relax
    echo "ENDMDL"   >> $AA_path_pdb_relax

done

# Post-process all-atom output file
gawk -i inplace '{if($1!="TER" && $1!="END") print $0}' $AA_path_pdb
gawk -i inplace '{if($1!="TER" && $1!="END") print $0}' $AA_path_pdb_relax

# Change the format to PLUMED
gawk '{if($1=="ATOM" || $1=="HETATM")printf "%s%s%s\n", substr($0,1,54), "  1.00  1.00", substr($0,67,length($0)-66); else print $0;}' $AA_path_pdb > tmp1.pdb
gawk '{if($1!="MODEL")print $0;}' tmp1.pdb > tmp2.pdb
sed -i 's/ENDMDL/END/g' tmp2.pdb
mv tmp2.pdb ${AA_path_pdb/models/plumed}

# Delete temporary files
rm -f frame* tmp?.* mdinfo xmin.in xmin.out
