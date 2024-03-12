# FR file
#transform='/gm2data/tbarrett/gm2fr/results/Run4/Nominal/transform.root'
transform='/gm2data/tbarrett/gm2fr/studies/run23_bunch_corrections/CorrectedFourier_Run3b.root'
# types of formats
formatA='./formats/formatA_sBin_constraintOff.txt'
formatB='./formats/formatB_sBin_constraintOn.txt'
formatC='./formats/formatC_mBin_constraintOff.txt'
formatD='./formats/formatD_mBin_constraintOn.txt'
formatE='./formats/formatE_mBin_mop_constraintOn_wcboOn.txt'

run4='/gm2data/cornell/histograms/aMethod/histogram_a_run4.root'
run5='/gm2data/cornell/histograms/aMethod/histogram_a_run5.root'

for dat in noRF ;  do # dat = run number

    # full-fit on calo sum, kevin-style 
    cE=0.0;
    seedNo=0;
    for caloNum in {0..24} ; do
      ./fitplot \
      -i /gm2data/cornell/histograms/aMethod/histogram_a_${dat}.root \
      -o sBin_constraintOn_cE${cE}_seed${seedNo}_${dat}_calo${caloNum}.root \
      -n $caloNum \
      -c ${cE} \
      -s ${seedNo} \
      -a 1 \
      -b $transform \
      -f $formatB | tee sBin_constraintOn_cE${cE}_seed${seedNo}_${dat}_calo${caloNum}.log 
    done
    wait

done

# full-fit on calo sum, kevin-style but no constraint
cE=0.0;
seedNo=0;

# run 4 only
#./fitplot -i $run4 -o sBin_constraintOff_cE${cE}_seed${seedNo}_run4_calo${caloNum}.root -n $caloNum -c ${cE} -s ${seedNo} -a 1 -b $transform -f $formatA 2>&1 > sBin_constraintOff_cE${cE}_seed${seedNo}_run4_calo${caloNum}.log &
# run 5 only
#./fitplot -i $run5 -o sBin_constraintOff_cE${cE}_seed${seedNo}_run5_calo${caloNum}.root -n $caloNum -c ${cE} -s ${seedNo} -a 1 -b $transform -f $formatA 2>&1 > sBin_constraintOff_cE${cE}_seed${seedNo}_run5_calo${caloNum}.log &



## full-fit on calo sum, zep-style with no constraint and no cE
#
#cE=0.0;
#seedNo=0;
#
#./fitplot -i $run4a -o mBin_constraintOff_cE${cE}_seed${seedNo}_run4_calo${caloNum}.root -n $caloNum -c ${cE} -s ${seedNo} -a 1 -b $run4transform -f $formatC 2>&1 > mBin_constraintOff_cE${cE}_seed${seedNo}_run4_calo${caloNum}.log &
#
#
## full-fit on calo sum, zep-style with no constraint and cE
#
#cE=1.0;
#seedNo=0;
#
#./fitplot -i $run4a -o mBin_constraintOff_cE${cE}_seed${seedNo}_run4_calo${caloNum}.root -n $caloNum -c ${cE} -s ${seedNo} -a 1 -b $run4transform -f $formatC 2>&1 > mBin_constraintOff_cE${cE}_seed${seedNo}_run4_calo${caloNum}.log &
#
#
# full-fit on calo sum, zep-style with constraint and no cE

cE=0.0;
seedNo=0;

#./fitplot -i $run4a -o mBin_constraintOn_cE${cE}_seed${seedNo}_run4_calo${caloNum}.root -n $caloNum -c ${cE} -s ${seedNo} -a 1 -b $run4transform -f $formatD 2>&1 > mBin_constraintOn_cE${cE}_seed${seedNo}_run4_calo${caloNum}.log &


# full-fit on calo sum, zep-style with constraint and cE

cE=1.0;
seedNo=0;

#./fitplot -i $run4a -o mBin_constraintOn_cE${cE}_seed${seedNo}_run4_calo${caloNum}.root -n $caloNum -c ${cE} -s ${seedNo} -a 1 -b $run4transform -f $formatD 2>&1 > mBin_constraintOn_cE${cE}_seed${seedNo}_run4_calo${caloNum}.log &


# full-fit on calo sum, zep-style with constraint and cE and wCBO

cE=1.0;
seedNo=0;

#./fitplot -i $run4a -o mBin_constraintOn_cE${cE}_wcboOn_seed${seedNo}_run4_calo${caloNum}.root -n $caloNum -c ${cE} -s ${seedNo} -a 1 -b $run4transform -f $formatE 2>&1 > mBin_constraintOn_cE${cE}_wcboOn_seed${seedNo}_run4_calo${caloNum}.log 


