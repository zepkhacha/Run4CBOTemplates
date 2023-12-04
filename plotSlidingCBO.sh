transform='/gm2data/tbarrett/gm2fr/studies/run23_bunch_corrections/CorrectedFourier_Run3b.root'

runlist="
4
"

for run in ${runlist};
do

    for calo in {1..4};
    do
        echo "run $run calo $calo" 
        ./cboAmplitude -s slidingFits/run${run}_calo${calo}_windowFits.root \
        -f fullFits/run${run}/sBin_constraintOn_cE0.0_seed0_run${run}_calo${calo}.root \
        -o run${run}_calo${calo}_cboPlots
        
    done #end for loop over calo
    wait

done # end for loop over run
