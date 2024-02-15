fileList="
noRF_calo1_windowFits.root
noRF_calo2_windowFits.root
noRF_calo3_windowFits.root
noRF_calo4_windowFits.root
noRF_calo5_windowFits.root
noRF_calo6_windowFits.root
noRF_calo7_windowFits.root
noRF_calo8_windowFits.root
noRF_calo9_windowFits.root
noRF_calo10_windowFits.root
noRF_calo11_windowFits.root
noRF_calo12_windowFits.root
noRF_calo13_windowFits.root
noRF_calo14_windowFits.root
noRF_calo15_windowFits.root
noRF_calo16_windowFits.root
noRF_calo17_windowFits.root
noRF_calo18_windowFits.root
noRF_calo19_windowFits.root
noRF_calo20_windowFits.root
noRF_calo21_windowFits.root
noRF_calo22_windowFits.root
noRF_calo23_windowFits.root
noRF_calo24_windowFits.root
"

for file in $fileList; do
  runcalo=${file%_window*}
  run=${runcalo%_calo*}
  calo=${runcalo#run*_}
  fullfitfilename=sBin_constraintOn_cE0.0_seed0_${runcalo}.root
  #fullfitfilename=sBin_constraintOff_cE0.0_seed0_${run}_calo.root
  ./windowFFT -i slidingFits2/${file} -o fft_${file%.root} -f templateFits2/${run}/$fullfitfilename
done

