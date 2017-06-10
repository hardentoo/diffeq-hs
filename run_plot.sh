#!/usr/bin/env bash

echo "Emptying plots/"
rm -R ./plots
mkdir ./plots

echo Running "$1"

TIME="$(time (stack runghc "$1"))"

open ./plots/*
