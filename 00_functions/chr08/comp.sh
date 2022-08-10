#!/bin/bash

for i in "$@"/*
do
  for j in "$@"/*
  do
    if [ "$i" \< "$j" ]
    then
     diff -s $i $j | grep -v '^Only in '
    fi
  done
done

