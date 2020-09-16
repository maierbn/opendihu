#!/bin/bash

if [[ $@ == *"03_generator_vc"* ]]; then
  echo using g++ on $@
  g++ $@
else
  echo using scorep on $@
  scorep --user --instrument-filter=scripts/compile-time-filter.filt g++ $@
fi
