#!/bin/bash
if [ "X$1" == "X" ]; then
  echo "provide directory for new job"
  exit 1
fi
mkdir $1
cp -r $LOCAL_MCatNLODIR/MCatNLO_4_1_0_Base/* $1

