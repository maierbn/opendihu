
# This script declares to SCons how to compile the example.
# It has to be called from a SConstruct file.
# The 'env' object is passed from there and contains further specification like directory and debug/release flags.
#
# Note: If you're creating a new example and copied this file, adjust the desired name of the executable in the 'target' parameter of env.Program.


Import('env')     # import Environment object from calling SConstruct

# create the main executable
env.Program(target = 'muscle', source = "src/muscle.cpp")
env.Program(target = 'only_mechanics_muscle', source = "src/only_mechanics_muscle.cpp")
env.Program(target = 'tendon', source = "src/tendon.cpp")
env.Program(target = 'linear_tendon', source = "src/linear_tendon.cpp")

# env.Program(target = 'linear_tendon_precice', source = "src/linear_tendon_precice.cpp")

