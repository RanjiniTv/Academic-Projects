#!/bin/bash
#Command line arguments include Simulator file, Memory trace file, Cache Size, Cache line size in bytes, # ways
#configured for Cache Size=1MB ,Cache line size=64 ,#ways=16
python3 cacheProject.py $1 1048576 64 16
