#!/bin/bash

pushd ../
julia --project -e "using Pkg; Pkg.test()"
