#!/bin/bash

f2py -c LB2D.f95 -m LB2D --f90flags="-O2"
