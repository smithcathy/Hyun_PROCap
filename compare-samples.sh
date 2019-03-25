#!/usr/bin/env bash

date
#for dist in {50..150..5}
#do
	#for file in 48646_R1_cap.clip.TGACCA.hs37d5.bowtie.bed.gz 51644_R1_cap.clip.CAGATC.hs37d5.bowtie.bed.gz 48641_R1_cap.clip.ACAGTG.hs37d5.bowtie.bed.gz
	#do
	
		#gzip -cd $file | Rscript CompareSamples.r $dist ${file%%_*} &
	#done
#wait
#echo first loop $dist
#done
#wait
parallel gzip -d ::: SamplesToCompare/*
for dist in {50..150..5}
do
	for filei in SamplesToCompare/ReadsAtLoci_48646_$dist SamplesToCompare/ReadsAtLoci_48641_$dist 
	do
		for filej in SamplesToCompare/ReadsAtLoci_48641_$dist SamplesToCompare/ReadsAtLoci_51644_$dist
		do
			if [ "${filei##*/}" == "${filej##*/}" ]
			then
				continue
			fi
			Rscript MergeSamples.r $filei $filej &
		done
	done
#wait
#echo second loop $dist
done
wait
parallel gzip ::: SamplesToCompare/*
wait
date
