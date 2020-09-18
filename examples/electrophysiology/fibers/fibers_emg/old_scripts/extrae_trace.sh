#!/bin/sh

# Extrae configuration file
export EXTRAE_CONFIG_FILE=./extrae.xml

# Select tracing library (depends on runtime & lang)

# For MPI/C apps
export LD_PRELOAD=$EXTRAE_HOME/lib/libmpitrace.so # eems like this is also needed if extae is statically linked
# For MPI/Fortran apps
#export LD_PRELOAD=$EXTRAE_HOME/lib/libmpitracef.so
# For OpenMP apps
#export LD_PRELOAD=$EXTRAE_HOME/lib/libomptrace.so
# For MPI+OpenMP/C apps
#export LD_PRELOAD=$EXTRAE_HOME/lib/libompitrace.so
# For MPI+OpenMP/Fortran apps
#export LD_PRELOAD=$EXTRAE_HOME/lib/libompitracef.so

# Run the application
$@
