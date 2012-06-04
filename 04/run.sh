#!/bin/bash

rm files/*

mpirun -np 6 ./sim
