#!/bin/bash

cat > xmin.in <<EOF
 input for sander for xmin minimization
 &cntrl
   cut     = 999,
   rgbmax  = 999,
   ntx     = 1,
   irest   = 0,
   ipol    = 0,
   ntb     = 0,
   igb     = 8,
   imin    = 1,
   maxcyc  = 5000,
   ntmin   = 3,
   drms    = 5e-2,
   ntpr    = 0,
   ntr     = 1,
   restraint_wt = 100.0,
   restraintmask = "$5",
 /
 &lmod
   lbfgs_memory_depth = 3,
   xmin_method = 'LBFGS',
   xmin_verbosity = 0,
 /
EOF

mpirun -np $6 --oversubscribe ${AMBERHOME}/bin/sander.MPI -O -i xmin.in -o xmin.out \
   -c $1 -ref $2 -p $3 -r $4
