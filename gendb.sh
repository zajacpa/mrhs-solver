#!/bin/bash

#params: $1 - number of ANDs, $2 - number of key bits, $3 - density, $4 - number of instances

i=1
a=1

while [ $i -le $4 ]
do
  echo Number: $i
  lst=(`bin/mrhs-v0 -P -m $1 -l $2 -k $2 -d $3 -s $a -o tmp.mrhs`)
  
  #check whether solution exists
  ok=${lst[6]}
  echo $ok
  
  #ok, solution found, move system to db
  if [ $ok = '1' ]; then
    cp tmp.mrhs db/$1.$2.$3.$i.mrhs
    ((i++))
  fi
  
  rm tmp.mrhs
  ((a++))
done


