#sudo python setup.py install --record file.txt
rm -rf dist/scons_config-0.1-py2.7.egg
zip -r dist/scons_config-0.1-py2.7.egg sconsconfig dist/EGG-INFO
rm -rf dist/scons_config-0.1-py3.5.egg
zip -r dist/scons_config-0.1-py3.5.egg sconsconfig dist/EGG-INFO
