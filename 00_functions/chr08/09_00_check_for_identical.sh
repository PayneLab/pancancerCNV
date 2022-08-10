#!/bin/bash

# Given a directory, see if any of the files in the directory are identical to each other
# Usage: comp.sh [dir_to_test]
# E.g. comp.sh run5/tmp000/ 

for i in "$@"*
do
  for j in "$@"*
  do
    if [ "$i" \< "$j" ]
    then
     cmp -s "$i" "$j" && echo "Files $i and $j are identical"
    fi
  done
done

