for CALO in {0..24}; do
  screen -dmS calo${CALO}_noRF bash -c "source ~/.gm2Setup; ./runSlidingFits.sh noRF ${CALO}"
done

