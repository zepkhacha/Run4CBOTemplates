for calo in {1..4}; do
  for run in 4 ; do
    hadd -f run${run}_calo${calo}_windowFits.root run${run}/run${run}_calo${calo}_window{0000..1300..20}.root
  done
done
