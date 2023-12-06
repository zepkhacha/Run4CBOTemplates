for CALO in {1..24}; do
  screen -dmS calo${CALO}_run4 bash -c "source ~/.gm2Setup; ./runSlidingFits.sh 4 ${CALO}"
  #screen -dmS calo${CALO}_run5 bash -c "source ~/.gm2Setup; ./runSlidingFits.sh 5 ${CALO}"
done

