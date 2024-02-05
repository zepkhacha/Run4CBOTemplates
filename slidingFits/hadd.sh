for calo in {1..24}; do
  for run in noRF ; do
    hadd -f ${run}_calo${calo}_windowFits.root ${run}/${run}_calo${calo}_window{0000..2600..40}.root 
  done
done

hadd -f ${run}_windowFits.root ${run}/${run}_calo{1..24}_window{0000..2600..40}.root 
