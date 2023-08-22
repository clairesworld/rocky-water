#!/bin/bash

for d in */ ; do
 cd $d
 #for i in ?????*bar; do
  # mv -i "$i" "${i%?????}bar"
 #done
 for i in ?????*bar; do
  mv $i "${i:0:5}bar"
 done
 #for j in ?????*bar; do
  # mv -- "$j" "${j//,/}"
 #done
 cd ../
 
done
