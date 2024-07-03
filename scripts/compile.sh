#!/bin/bash

g++ ./bktree/bktree.cpp -Ofast -march=native -shared -Wall -fPIC $(python -m pybind11 --includes) -o ./bktree/bktree$(python3-config --extension-suffix) -std=c++17
