-- ================================================================================================ --
--                                         Basic parameters                                         --
-- ================================================================================================ --
-- Template of configuration file for oneview module

-- Name of the experiment specified in report main pages
experiment_name = nil

-- Name of the binary file to analyze
binary         = "/lustre/cray/ws9/2/ws/icbbnmai-opendihu1/opendihu-hawk-gnu-openmpi/examples/electrophysiology/fibers/fibers_contraction/no_precice/build_release/biceps_contraction"

-- List of external libraries to analyze
external_libraries = {
  "build_release/lib/hodgkin_huxley_1952_fast_monodomain.so", "libpetsc.so.3.12.2", "libc-2.28.so", "libpetsc.so"
-- "lib.so", "lib.so"
}
-- List of external categories to analyze in LPROF.
-- If nil, one category is created for each entry in external_libraries.
-- If not nil, only categories described by the user are analyzed.
-- A category is a subtable with following fields:
--   + type - "library" to create a category corresponding to a library,
--            "library_group" to create a category that accumulate several libraries,
--            "all-external-libraries" to reproduce the default behavior with a single entry
--   + value - "library name" if type is "library",
--             A table containg libraries names to accumulate: {"lib1.so", "lib2.so"},
--             nil if type is "all-external-libraries"
--   + name - Name of the category when type is "library_group"
-- /!\ Warning /!\ If a single library is used in several library groups, some errors will appear in reports
-- custom_categories = {
--   {type = "all-external-libraries"},
--   {type = "library", value = "lib.so"},
--   {type = "library_group", value = {"lib.so", "lib1.so"}, name = "my_libs"},
-- }

-- Path to the dataset directory used by the application.
-- The dataset is considered as a directory containing all files
-- needed to run the application. It is copied / linked into
-- the experiment directory according to dataset_handler value.-- Use an empty string or remove the declaration if no 
-- dataset is required by the application
dataset        = ""

-- How dataset is handled. "link" (default) to create a symbolic link
-- from the experiment directory to the specified dataset, "copy" to
-- duplicate the dataset directory into the experiment directory
dataset_handler= "link"

-- Command used to run the binary. 
-- + <binary> will be replaced by the path to the binary file to run 
-- If some parameters are used, specify them in the string.
--   example: "<binary> -n4"
-- Use an empty string or remove the declaration if no 
-- specific command is required to run the application
run_command    = "<binary> /lustre/cray/ws9/2/ws/icbbnmai-opendihu1/opendihu-hawk-gnu-openmpi/examples/electrophysiology/fibers/fibers_contraction/no_precice/settings_biceps_contraction.py ramp.py --end_time=2"

-- Script to use with job scheduler. 
-- If your application must be run on a cluster using a job scheduler (ex. SLURM), fill
-- this field with the path to a script which will be run using
-- the command detailed in 'batch_command'. The script must have been modified to replace
-- the application executable and its arguments with keyword <run_command>.
-- This configuration file must be filled with options to run your application.
-- The number of processes to use can be referred as <number_processes>.
-- The number of nodes to use can be referred as <number_nodes>.
-- The number of tasks per node to use can be referred as <number_tasks_nodes>.
-- It can be set to an empty string or removed if it not used.
-- batch_script    = ""

-- Command to use to submit job to the job scheduler. 
-- If your application must be run on a cluster using a job scheduler (ex. SLURM), fill
-- this field with the submission command. The script passed to the job scheduler
-- must be referred as <batch_script>.
-- Example using sbatch command (SLURM): batch_command = "sbatch <batch_script>"
-- batch_command   = ""

-- Number of MPI processes.
-- The variable can be referred as <number_processes> in field 'mpi_command' and
-- in the batch script.
number_processes = 2

-- Number of nodes to uses in an MPI application.
-- The variable can be referred as <number_nodes> in field 'mpi_command' and
-- in the batch script.
-- number_nodes = 1

-- Number of tasks per node to uses in an MPI application.
-- The variable can be referred as <number_tasks_nodes> in field 'mpi_command' and
-- in the batch script.
-- number_tasks_nodes = 1

-- Part of the command used to run the binary.
-- It will be added at the beginning of the final command. It is used to specify
-- a launcher for the application, such as mpirun for MPI applications
--   example: "mpirun -n <number_processes>"
-- Use an empty string or remove the declaration if no 
-- specific command prefix is required to run the application
mpi_command    = "mpirun -n <number_processes>"

-- Define the corresponding OpenMP variable to set the maximal number of threads.
-- If the application uses OpenMP, the value must be greater than 1.
-- If the application does not use OpenMP, the value can be equal to 1, nil
-- or the field can be removed from the configuration file.
-- omp_num_threads= 1

-- Part of the command used to run the binary.
-- It will be added at the begining of the final command, before mpi_command.
-- It is used to specify how to pin the binary and is added untouched in the command.
-- pinning_command    = ""

-- Directory where the binary must be run
-- + <dataset> will be replaced by the dataset directory located into the
--   experiment directory.
-- It is needed when the application must be run from a specific directory in
-- the dataset directory.
-- Use an empty string or remove the declaration if no 
-- specific directory is required to run the application
run_directory  = "/zhome/academic/HLRS/cbm/icbbnmai/workshop_2020/opendihu/examples/electrophysiology/fibers/fibers_contraction/no_precice/build_release"

-- ================================================================================================ --
--                            Filter used to select loops to analyze                                --
-- ================================================================================================ --
-- !! Uncomment the filter you want to use !!
-- !! If no filter is specified, the ten hottest loops are dynamically analyzed !!

-- This filter uses the first <value> loops, ordered by coverage
filter = {
   type = "number",
   value = 10,
}

-- This filter uses all loops whose coverage is greater than <value> (in percentage)
--filter = {
--   type = "coverage",
--   value = 1,
--}

-- This filter uses all loops while the cumulated coverage is lower than <value> (in percentage)
--filter = {
--   type = "cumulated_coverage",
--   value = 1,
--}

-- This filter uses all loops.
--filter = {
--   type = "all",
--}

-- ================================================================================================ --
--                            Specify when the profiling should start                               --
-- ================================================================================================ --
-- !! Uncomment the table you want to use !!
-- !! If no table is specified, the default table has 'p' as unit and 30 as value !!
-- !! Report ONE always analyzes all loops !!

-- If the profiling should begin after a given percentage of the application time.
-- A first run of the application is automatically performed to time the application.
profile_start = {
   unit = "p",      -- 'p' for 'percentage'
   value = 30,         -- delay in percentage of the total application time
}
-- If the profiling should begin when the application is started
--profile_start = {
--   unit = "none",   -- Specify that no delay is needed
--   value = 0,         -- Useless value
--}

-- If the profiling should begin after a given time (in second)
--profile_start = {
--   unit = "s",      -- 's' for 'seconds'
--   value = 0,         -- delay in seconds
--}


-- ================================================================================================ --
--                                      Additional parameters                                       --
-- ================================================================================================ --
-- Frequencies to use in some dynamic analysis
-- Use an empty table if no frequency must be used.
-- !! Need permissions to change cpufreq files !!
-- !! Do not specify frequencies if MAQAO can not modify these files !!
-- frequencies    = {}

-- Table describing additional hardware counters analyzed.
-- Each entry is a list of hardware counters analyzed during a run.
-- Hardware counters must be separated with a comma ','.
-- Either hardware counter codes or names can be used.
-- If a set of hardware counters needs sudo permissions, set sudo_mode at true
-- !! Check that all hardware counters in a same list can be analyzed in a single run !!
--additional_hwc = {
--   {names="<names_1>", sudo_mode = <boolean>},
--   {names="<names_2>"},
--}

-- Exclude some areas (loops or blocks) from the experiment
-- type defines if the area is a loop ("loop") or a basic block ("block")
-- id is the MAQAO internal identifier
-- module represents the binary file containing the area to exclude.
-- The main binary is referred as "binary" (or a nil / empty string).
-- Another source must be referred by its name and must be declared in the external_libraries table.
--excluded_areas = {
--   {type = "loop", id = <MAQAO loop identifier>, module = ""},
--}

-- Include some areas (loops or blocks) from the experiment in addition of what is analyzed using filters.
-- type defines if the area is a loop ("loop") or a basic block ("block")
-- id is the MAQAO internal identifier
-- module represents the binary file containing the area to include.
-- The main binary is referred as "binary" (or a nil / empty string).
-- Another source must be referred by its name and must be declared in the external_libraries table.
--included_areas = {
--   {type = "loop", id = <MAQAO loop identifier>, module = ""},
--}

-- DECAN analyzes several variants in a single run.
-- If decan_multi_variant is set to false, only one variant will be analyzed per run.
-- The experiment will be longer.
-- decan_multi_variant = false

-- If false, DECAN analyzes variants in several sets and, what takes more time 
-- decan_all_variants = true

-- Specify the location of the source code.
-- ONE-View looks for source code in directory specified in debug data. However when the
-- application is not compiled on the same machine than the run, source code is not located
-- in the same directory. This variable can be used to specify the source code of the application.
-- When used, the last subdirectory given must appears in the initial location. For example,
-- if sources where located in "/apps/my_app/src/", source_code_location can be
-- "/home/user/sources/my_app" as "my_app" appears also in the original location,
-- but it can not be "/home/user/sources/my_new_app" as "my_new_app" does not appears
-- in debug data. However, it can be "/home/user/sources/my_new_app/src".
-- source_code_location = ""

-- Specify maximal number of path a loop can have.
-- If a loop has more paths, it will not be analyzed by ONEVIEW. To analyze loops regardless of
-- the number of paths, the value must be 0. Default value is 4.
maximal_path_number = 4

-- Specify if sudo mode can be used during experiments.
-- Sudo mode is needed to analyze some hardware counters in report THREE or change the frequency
-- in reports TWO and THREE.
is_sudo_available = false

-- Specify additional options for LPROF.
-- If not empty, everything specified in this field is passed to LPROF through the command line
lprof_params = ""

-- Specify additional options for LPROF post processing phase.
-- If not empty, everything specified in this field is passed to LPROF through the lprof_opts table
lprof_post_process_params = {
}

-- Specify additional options for VPROF.
-- If not empty, everything specified in this field is passed to VPROF through the command line
vprof_params = ""

-- Specify additional options for DECAN.
-- If not empty, everything specified in this field is passed to DECAN through the command line
decan_params = ""

-- Specify additional options for CQA.
-- If not empty, everything specified in this field is passed to CQA through the CQA context
cqa_params = {}

-- If true, it specify that the binary should not be copied into the experiment directory but keep in place
-- keep_binary_location = true

-- Specify runs parameters to use when the scalability report or the compare report must be generated.
-- Each entry in the table describes values to override for a run:
--   + number_processes - its number of processes (default is 1),
--   + omp_num_threads - its number of threads (default is 1),
--   + number_nodes - its number of nodes (default is 1),
--   + number_tasks_nodes - its number of tasks per node (default is 1),
--   + run_command - a command specific for the run (default is general value),
--   + run_directory - a directory where the application must be run (default is general value),
--   + mpi_command - a mpi command specific for the run (default is general value),
--   + dataset - a specific dataset (default is general value),
--   + script_variables - some user defined variables used for batch-script,
--   + profile_start - when the profiling should start,
--   + name - the run name (default is 'run_<index>', with <index> the index of the entry in the table.
multiruns_params = {
--  {},
}

-- Specify which run will be used as reference when scalability reports are generated.
-- If two entries match the filter, the first one described is used as refences (main entry, then
-- the entry with the lowest index in multirun-params table). Availables values are:
--   + main (default) - Uses the run specified by main parameters (and not in multirun-params) as reference
--   + lowest-time - Uses the shortest run in time as reference
--   + highest-time - Uses the longest run in time as reference
--   + lowest-threads - Uses the run with the lowest number of threads analyzed using LPROF as reference
--   + highest-threads - Uses the run with the highest number of threads analyzed using LPROF as reference
--   + lowest-efficiency - Uses the run with the worst efficiency as reference
--   + highest-efficiency - Uses the run with the better efficiency as reference
--   + <number> - Uses the entry at index <number> in the multirun-param table as reference.
-- scalability_reference = "main"

-- Specify environment variables to set before running the application.
environment_variables = {
--  {name = <name>, value = <value>},
}

-- Specify user defined variables that are replaced in the batch_script file.
script_variables = {
--  <name> = <value>,
}

-- Specify the path of QPlot directory containing the file parse.js. Qplot is run with the command
-- $ nodejs <dir>/parse.js <oneview-report.xlsx>
qplot_script = nil

-- Specify a threshold to apply to threads to consider them as valid threads and
-- use them in reports. Value is a string composed of a number and an unit: s for seconds
-- and % for percentage of the total walltime (default unit). For example "5s" or "20%" or 10 (="10%").
thread_filter_threshold = nil

-- Specify a threshold to apply to objects coverages (loops anf functions) to include them in static analysis.
-- All objects whose coverage is lower than the given value will not be statically and dynamically analized
object_coverage_threshold = nil

