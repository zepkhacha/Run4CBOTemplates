for calo in {0..24}; do
  for run in noRF ; do
    hadd -f -O ${run}_calo${calo}_windowFits.root \
        ${run}/${run}_calo${calo}_window0{000..648..18}.root \
        ${run}/${run}_calo${calo}_window{0648..2600..70}.root
  done
done

run=noRF
hadd -f -O ${run}_windowFits.root ${run}/${run}_calo{0..24}_window0{000..648..18}.root ${run}/${run}_calo{0..24}_window{0648..2600..70}.root
