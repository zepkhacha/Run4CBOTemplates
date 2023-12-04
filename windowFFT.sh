fileList="
run4_calo1_windowFits.root
run4_calo2_windowFits.root
run4_calo3_windowFits.root
run4_calo4_windowFits.root
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

