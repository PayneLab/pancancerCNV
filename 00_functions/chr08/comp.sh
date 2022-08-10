#!/bin/bash

for i in "$@"/*
do
  for j in "$@"/*
  do
    if [ "$i" \< "$j" ]
    then
     cmp -s $i $j && echo "Files $i and $j are identical"
    fi
  done
done

