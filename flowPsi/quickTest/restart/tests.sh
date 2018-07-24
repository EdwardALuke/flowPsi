#!/bin/bash
export LD_LIBRARY_PATH=$LOCI_BASE/lib:$FLOWPSI_BASE/lib:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=$LOCI_BASE/lib:$FLOWPSI_BASE/lib:$DYLD_LIBRARY_PATH

BC_TESTS=$*

for ii in $BC_TESTS ; do
  i=`echo $ii | sed s/.vars//`
  cat $i.txt
  cat $i.txt > $i.test
  mkdir TEST_$i
  EXTRACT_ARGS=`cat $i.xtr.dat`
  TOL_ARGS=`cat $i.tol.dat`
  GRID_DAT=`cat $i.vog.dat`
  EXTRACT_TIME_1=`cat $i.info.dat | sed "s/.*start=\([^;]*\).*/\1/"`
  EXTRACT_TIME_2=`cat $i.info.dat | sed "s/.*end=\([^;]*\).*/\1/"`
  FLOWPSI_ARG=`cat $i.info.dat | sed "s/.*flowarg=\([^;]*\).*/\1/"`
#  echo EXTRACT_TIME_1=$EXTRACT_TIME_1 EXTRACT_TIME_2=$EXTRACT_TIME_2 FLOWPSI_ARG=$FLOWPSI_ARG
  cd TEST_$i
  ln -s ../$i.vars .
  ln -s ../$GRID_DAT $i.vog

  $FLOWPSI $i > flow.out1
  $EXTRACT -ascii $i $EXTRACT_TIME_1 $EXTRACT_ARGS > $i.1.dat
  rm -fr output
  $FLOWPSI $i $FLOWPSI_ARG > flow.out2
  $EXTRACT -ascii $i $EXTRACT_TIME_2 $EXTRACT_ARGS > $i.2.dat
  cd .. ;
  if $NDIFF TEST_$i/$i.1.dat TEST_$i/$i.2.dat $TOL_ARGS; then 
     echo Test Passed `pwd` $i
     echo Passed `pwd` $i >> $i.test
     rm -fr TEST_$i
  else 
     echo Test `pwd` $i FAILED !!!!!!!
     echo Test `pwd` $i FAILED >> $i.test
  fi

done
