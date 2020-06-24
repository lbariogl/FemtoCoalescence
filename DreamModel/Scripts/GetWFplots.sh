#!/bin/bash

WorkDirectory=${1}
InputDir=${2}
echo Working here ${WorkDirectory} and getting stuff  from ${InputDir} 
myArray=(5 0)
declare -a arr=(D Q)
for s in $(seq 0 1); do
   for wave in $(seq 0 1); do
		for kBin in $(seq 1 26); do 
		lval=$((${myArray[${wave}]}))
		OutDir=${WorkDirectory}/Spin_${arr[${s}]}/Wave_$((${myArray[${wave}]}))
		if [ ! -d ${OutDir} ]; then
		    mkdir -p ${OutDir}
		fi
		FileName=${InputDir}/Wf-pd-${arr[${s}]}/Wf-pd-${arr[${s}]}-k$((${kBin}*5))-L400-l${wave}p$((${myArray[${wave}]})).dat

		${HOME}/Coalescence/install/DreamModel/excutePdWF $((${kBin}*5)) ${wave} $((${myArray[${wave}]})) ${FileName} ${WorkDirectory} ${arr[${s}]} 
		mv WF_${arr[${s}]}_k$((${kBin}*5))_l${wave}p$((${myArray[${wave}]})).png ${OutDir}/WF_${arr[${s}]}_k$((${kBin}*5))_l${wave}p$((${myArray[${wave}]})).png 
		#echo 0 1 0 0 0 ${FileName} hCk_RebinnedMeV_1 ${OutDir}
		done 
  done

done
