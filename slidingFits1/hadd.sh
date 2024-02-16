#for calo in {1..24}; do
#  for run in noRF ; do
#    hadd -f -O ${run}_calo${calo}_windowFits.root \
#        ${run}/${run}_calo${calo}_window00{00..70..35}.root \
#        ${run}/${run}_calo${calo}_window{0140..2600..70}.root
#  done
#done

run=noRF
hadd -f -O ${run}_windowFits.root \
    ${run}/${run}_calo{1..24}_window00{00..70..35}.root \
    ${run}/${run}_calo{1..24}_window{0140..2600..70}.root
