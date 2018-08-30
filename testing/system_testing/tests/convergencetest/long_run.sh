cd Euler/scripts && . long_actualize_settings.sh &&  cd ../../
cd Euler/build_debug && . long_run_convergence_test_euler.sh && cd ../../Heun/build_debug && . long_run_convergence_test_heun.sh && cd ../..
cd Euler/scripts && pyod long_measure_error.py && yes | cp -rf long_fig.pdf ../../ && cd ../..
# if you get an error here, replace 'pyod' by 'python' or set an alias "pyod='<YOUR_SYSTEM_PYTHON>'" in your .bashrc file.
