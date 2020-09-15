

EXAMPLE_HOME=$OPENDIHU_HOME/examples/electrophysiology/fibers/fibers_contraction/no_precice
#mpirun -n 2 ./biceps_contraction ../settings_biceps_contraction.py ramp.py --end_time=20

module load maqao

maqao oneview create-report=one config=maqao_config.lua xp=ov1
