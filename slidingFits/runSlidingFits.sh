run=${1}
calo=${2}

if [[ -z "${1}" || -z "${2}" ]]; then
  echo "usage: ./runSlidingFits.sh [run] [calo]"
  exit 1
fi

# FR file
frFile='/gm2data/tbarrett/gm2fr/studies/run23_bunch_corrections/CorrectedFourier_Run3b.root'

# types of formats
formatA="/gm2data/zkhechad/run45/envelope/formats/formatA_sBin_constraintOff.txt"
formatB="/gm2data/zkhechad/run45/envelope/formats/formatB_sBin_constraintOn.txt"
formatC="/gm2data/zkhechad/run45/envelope/formats/formatC_mBin_constraintOff.txt"
formatD="/gm2data/zkhechad/run45/envelope/formats/formatD_mBin_constraintOn.txt"

mainDirectory="/gm2data/zkhechad/cboTemplates/Run4CBOTemplates/"
mkdir -p ${mainDirectory}/slidingFits/run${run}
outputDirectory="${mainDirectory}/slidingFits/run${run}/"

# -p points to kevin-style fits
#-q ${mainDirectory}fullFits/run${run}/sBin_constraintOn_cE0.0_seed0_run${run}_calo${calo}.root \

echo "performing window 0"
./slidingwindowfitplot -p ${mainDirectory}/cboIsolate/run${run}_calo${calo}.root \
-a 1 \
-W ${outputDirectory}/run${run}_calo${calo}_window1300.root \
-w 0 -n ${calo} -q ${mainDirectory}fullFits/run${run}/sBin_constraintOn_cE0.0_seed0_run${run}_calo${calo}.root \
-i /gm2data/cornell/histograms/aMethod/histogram_a_run${run}.root \
-o ${outputDirectory}/run${run}_calo${calo}_window0000.root \
-c 0 -s 0 -a 0 -b ${frFile} -f ${formatB} 2>&1 > ${outputDirectory}/run${run}_calo${calo}_window0000.log 
echo "done performing window 0" 

step=20
start=${step}
end=2600

for window in $(seq ${start} ${step} ${end}); do
   
   printf -v windowLabel "%04d" ${window}
   printf -v prevWindowLabel "%04d" $((window - step))
   echo "windowNo ${window}"

#-p ${mainDirectory}/cboIsolate/run${run}.root \
#-q ${mainDirectory}fullFits/run${run}/sBin_constraintOn_cE0.0_seed0_run${run}_calo${calo}.root \
   ./slidingwindowfitplot \
      -w ${window}  -n ${calo} -p ${mainDirectory}/cboIsolate/run${run}_calo${calo}.root \
      -a 1 \
      -q ${mainDirectory}fullFits/run${run}/sBin_constraintOn_cE0.0_seed0_run${run}_calo${calo}.root \
      -i /gm2data/cornell/histograms/aMethod/histogram_a_run${run}.root \
      -o ${outputDirectory}/run${run}_calo${calo}_window${windowLabel}.root \
      -c 0 \
      -s 0 \
      -a 0 \
      -b ${frFile} \
      -f ${formatB} \
      -W ${outputDirectory}/run${run}_calo${calo}_window1300.root \
      2>&1 > ${outputDirectory}/run${run}_calo${calo}_window${windowLabel}.log &

   wait
done 

#hadd -f ${mainDirectory}/slidingFits/run${run}_calo${calo}_windowFits.root ${outputDirectory}/run${run}_calo${calo}_window{0000..3050..10}.root
