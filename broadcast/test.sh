#!/bin/bash
for i in 2
do
#	for j in ../testcases/inp_$i"_"*
	for j in ../testcases/inp_$i"_400"
	do
		start=`date +%s%N`
		./prog < $j > ../testcases/result/"out_"$i"_"${j##*_}
		end=`date +%s%N`
		echo $j" finished in "$[($end  - $start)/1000000]" ms by METHOD "$i
	done
done
