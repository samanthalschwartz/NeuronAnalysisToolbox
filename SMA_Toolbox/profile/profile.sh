#!/bin/bash
OUTFILE=callgrind.out.2
valgrind -v --tool=callgrind --dump-instr=yes --callgrind-out-file=$OUTFILE $@ \
    && kcachegrind $OUTFILE &
