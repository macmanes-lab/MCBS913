#!/bin/bash

# Guanwen Wang
# help Junhong to kill all of his pids
# while working with Mosaik 

PIDS=$(ps aux | grep junhong | awk '{print $2}')
kill $PIDS
