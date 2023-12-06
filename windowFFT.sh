fileList="
run4_calo1_windowFits.root
run4_calo2_windowFits.root
run4_calo3_windowFits.root
run4_calo4_windowFits.root
run4_calo5_windowFits.root
run4_calo6_windowFits.root
run4_calo7_windowFits.root
run4_calo8_windowFits.root
run4_calo9_windowFits.root
run4_calo10_windowFits.root
run4_calo11_windowFits.root
run4_calo12_windowFits.root
run4_calo13_windowFits.root
run4_calo14_windowFits.root
run4_calo15_windowFits.root
run4_calo16_windowFits.root
run4_calo17_windowFits.root
run4_calo18_windowFits.root
run4_calo19_windowFits.root
run4_calo20_windowFits.root
run4_calo21_windowFits.root
run4_calo22_windowFits.root
run4_calo23_windowFits.root
run4_calo24_windowFits.root
"

for file in $fileList; do
  runcalo=${file%_window*}
  run=${runcalo%_calo*}
  calo=${runcalo#run*_}
  fullfitfilename=sBin_constraintOn_cE0.0_seed0_${runcalo}.root
  #fullfitfilename=sBin_constraintOff_cE0.0_seed0_${run}_calo.root
  #echo "file slidingFits/${file} fullfitfilename fullFits/$run/$fullfitfilename run ${run} calo ${calo}"
  ./windowFFT -i slidingFits/${file} -o fft_${file%.root} -f fullFits/${run}/$fullfitfilename
done

