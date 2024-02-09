transform='/gm2data/tbarrett/gm2fr/studies/run23_bunch_corrections/CorrectedFourier_Run3b.root'

runlist="
noRF
"

for run in ${runlist};
do

  for calo in 0; 
  do
      echo "run $run calo $calo" 
      ./dataDivFcn_noCBO \
      -i /gm2data/zkhechad/cboTemplates/Run4CBOTemplates/fullFits/sBin_constraintOn_cE0.0_seed0_${run}_calo${calo}.root \
      -o ${run}_calo${calo} \
      -b $transform 
      
  done #end for loop over calo
  wait

done # end for loop over run
