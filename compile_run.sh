#!/usr/bin/env bash

echo "Compilation:"
echo "============"
time (stack install)

echo ""

echo "Execution:"
echo "=========="
time diffeq-hs
