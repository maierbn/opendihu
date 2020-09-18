#!/bin/bash

#g++ $@

if true; then
if [[ $@ == *"03_generator_vc"* ]]; then
  echo using g++ on $@
  g++ $@
else
  echo using scorep on $@
  scorep-g++ $@
fi
fi
