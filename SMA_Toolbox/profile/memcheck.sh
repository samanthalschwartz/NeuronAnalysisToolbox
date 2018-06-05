#!/bin/bash

valgrind --tool=memcheck --suppressions=./.valgrind.supp --gen-suppressions=yes --leak-check=full --show-leak-kinds=all -v  $@
