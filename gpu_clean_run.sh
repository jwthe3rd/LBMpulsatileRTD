#!/bin/bash

make clean
make -j 6
./arterialRTD | tee run_log
