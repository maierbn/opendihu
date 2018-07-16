cd Euler/scripts && . actualize_settings.sh &&  cd ../../
cd Euler/build_debug && . run_convergence_test_euler.sh && cd ../../Heun/build_debug && . run_convergence_test_heun.sh && cd ../..
cd Euler/scripts && pyod measure_error.py && yes | cp -rf fig.pdf ../../ && cd ../..
# if you get an error here, replace 'pyod' by 'python' or set an alias "pyod='<YOUR_SYSTEM_PYTHON>'" in your .bashrc file.