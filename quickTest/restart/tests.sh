#!/bin/bash
export LD_LIBRARY_PATH=$LOCI_BASE/lib:$FLOWPSI_BASE/lib:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=$LOCI_BASE/lib:$FLOWPSI_BASE/lib:$DYLD_LIBRARY_PATH
export PATH=$FLOWPSI_BASE/bin:$PATH

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
  status1=$?
  $EXTRACT -ascii $i $EXTRACT_TIME_1 $EXTRACT_ARGS > $i.1.dat
  status2=$?
  rm -fr output
  $FLOWPSI $i $FLOWPSI_ARG > flow.out2
  status3=$?
  $EXTRACT -ascii $i $EXTRACT_TIME_2 $EXTRACT_ARGS > $i.2.dat
  status4=$?
  
  cd ..
  TESTNAME=`echo $PWD | sed s/.*quickTest[/]//`
  if [ $status1 -ne 0 -o $status2 -ne 0 -o $status3 -ne 0 -o $status4 -ne 0 ]; then
     echo Test $TESTNAME $i FAILED !!!!!!!
     echo Test $TESTNAME $i FAILED >> $i.test
  else
      if $NDIFF TEST_$i/$i.1.dat TEST_$i/$i.2.dat $TOL_ARGS; then 
	  cat $i.txt
	  echo Test $TESTNAME $i Passed  
	  echo Passed $TESTNAME $i >> $i.test
	  rm -fr TEST_$i
      else 
	  echo Test $TESTNAME $i FAILED !!!!!!!
	  echo Test $TESTNAME $i FAILED >> $i.test
      fi
  fi
done
