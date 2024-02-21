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
mkdir -p ${mainDirectory}/slidingFits0/${run}
outputDirectory="${mainDirectory}/slidingFits0/${run}/"

#-p ${mainDirectory}/cboIsolate/${run}_calo${calo}.root \
echo "performing window 0"
./slidingwindowfitplot -a 1 -w 0 -n ${calo} -q ${mainDirectory}fullFits/${run}/sBin_constraintOn_cE0.0_seed0_${run}_calo${calo}.root -i /gm2data/cornell/histograms/aMethod/histogram_a_${run}.root -o ${outputDirectory}/${run}_calo${calo}_window0000.root -c 0 -s 0 -b ${frFile} -f ${formatB} 2>&1 > ${outputDirectory}/${run}_calo${calo}_window0000.log 
echo "done performing window 0" 

step=70
start=${step}
end=2600

for window in $(seq ${start} ${step} ${end}); do
   
   printf -v windowLabel "%04d" ${window}
   printf -v prevWindowLabel "%04d" $((window - step))
   echo "windowNo ${window}"

#-p ${mainDirectory}/cboIsolate/run${run}.root \
#-q ${mainDirectory}fullFits/run${run}/sBin_constraintOn_cE0.0_seed0_run${run}_calo${calo}.root \
#-W ${outputDirectory}/run${run}_calo${calo}_window1300.root \
#-p ${mainDirectory}/cboIsolate/${run}_calo${calo}.root \

./slidingwindowfitplot -w ${window} -n ${calo} -a 1 -q ${mainDirectory}fullFits/${run}/sBin_constraintOn_cE0.0_seed0_${run}_calo${calo}.root -i /gm2data/cornell/histograms/aMethod/histogram_a_${run}.root -o ${outputDirectory}/${run}_calo${calo}_window${windowLabel}.root -c 0 -s 0 -b ${frFile} -f ${formatB} 2>&1 > ${outputDirectory}/${run}_calo${calo}_window${windowLabel}.log 

done 

#hadd -f ${mainDirectory}/slidingFits0/${run}_calo${calo}_windowFits.root ${outputDirectory}/${run}_calo${calo}_window{0000..2600..40}.root
