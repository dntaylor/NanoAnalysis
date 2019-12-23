#!/usr/bin/env bash
tar -zxf columnar.tar.gz
source columnar/bin/activate

echo "Running command:" $@
eval $@ || exit $?
