

EXAMPLE_HOME=$OPENDIHU_HOME/examples/electrophysiology/fibers/fibers_emg

#binary=$EXAMPLE_HOME/build_release/fast_fibers_emg
#$binary $EXAMPLE_HOME/settings_fibers_emg.py ramp_emg.py --end_time=10

module load maqao

maqao oneview create-report=one config=maqao_config.lua xp=ov
