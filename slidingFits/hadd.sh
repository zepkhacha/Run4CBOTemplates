for calo in {1..24}; do
  for run in 4 ; do
    hadd -f run${run}_calo${calo}_windowFits.root run${run}/run${run}_calo${calo}_window{0000..2600..40}.root 
  done
done

hadd -f run${run}_windowFits.root run${run}/run${run}_calo*_window{0000..2600..40}.root
