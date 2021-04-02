from python_settings.python_settings import *

# each dict entry corresponds to a cpp-template
# each template-dict can have a ordered list with template_arguments (assuming there are no template_arguments if omitted)

# possible outer templates (runnables) have the key runnable set to True (assuming False if ommited)

# templates that are discretizable in time have discretizableInTime set to True (assuming False if omited)
# "discretizableInTime" in template_arguments will get expanded to all classes, which have discretizableInTime == True

# templates that are a "TimeSteppingScheme" (e.g. all TimeSteppingScheme:: and OperatorSplitting::) have timeSteppingScheme set to True (assuming False if omited)
# "timeSteppingScheme" in template_arguments will get expanded to all classes, which have timeSteppingScheme == True

# templates with optional template_arguments can have the key template_arguments_needed set to the minimal required argument count
# template_arguments_needed is assumed to be len(template_arguments) if template_arguments_needed is omitted
# e.g. in BasisFunction::LagrangeOfOrder and PrescribedValues

# the keyword "Integer" can be used in template_arguments where an integer is expected (e.g. in CellmlAdapter)

# lists of the form [ "Mesh::" ] get auto expanded to [ "Mesh::StructuredRegularFixedOfDimension", "Mesh::Str..", ... ]

# templates added so far:
# TODO add postprocessing Postprocessing::ParallelFiberEstimation
# Postprocessing::StreamlineTracer
# PreciceAdapter::ContractionDirichletBoundaryConditions
# PreciceAdapter::ContractionNeumannBoundaryConditions
# PreciceAdapter::PartitionedFibers
# PreciceAdapter::MuscleContraction MuscleContractionSolver FastMonodomainSolver
# SpatialDiscretization::HyperelasticitySolver
# Control::MultipleInstances
# Control::Coupling
# Control::LoadBalancing
# Control::MapDofs
# OperatorSplitting::
# CellmlAdapter
# ModelOrderReduction::POD
# ModelOrderReduction::LinearPart
# ModelOrderReduction::ExplicitEulerReduced
# ModelOrderReduction::ImplicitEulerReduced
# FunctionSpace::
# OutputWriter::OutputSurface
# PrescribedValues
# TimeSteppingScheme::
# TimeSteppingScheme::DynamicHyperelasticitySolver
# TimeSteppingScheme::StaticBidomainSolver
# TimeSteppingScheme::MultidomainSolver
# TimeSteppingScheme::MultidomainWithFatSolver
# TimeSteppingScheme::QuasiStaticNonlinearElasticitySolverFebio
# TimeSteppingScheme::NonlinearElasticitySolverFebio
# TimeSteppingScheme::QuasiStaticLinearElasticitySolver
# TimeSteppingScheme::QuasiStaticNonlinearElasticitySolverChaste
# SpatialDiscretization::FiniteElementMethod
# Mesh::
# BasisFunction::
# Quadrature::
# Equation::
# Dummy

solver_common = [
    SettingsDictEntry("solverType", '"gmres"', 'the KSPType of the solver, i.e. which solver to use', 'solver.html#solvertype'),
    SettingsDictEntry("preconditionerType", '"none"', 'the preconditioner type of PETSc to use', 'solver.html#preconditionertype'),
    SettingsDictEntry("relativeTolerance", '1e-5', 'the relative tolerance of the residuum after which the solver is converged', 'solver.html#relativetolerance'),
    # undocumented
    SettingsDictEntry("absoluteTolerance", '0', 'absolute tolerance of the residual of the linear solver'),
    SettingsDictEntry("maxIterations", '1e4', 'the maximum number of iterations after which the solver aborts and states divergence', 'solver.html#maxiterations'),
    SettingsDictEntry("dumpFilename", '""',
                      "if this is set to a non-empty string, the system matrix and right hand side vector will be dumped before every linear solve", 'solver.html#dumpfilename'),
    SettingsDictEntry("dumpFormat", '"default"', 'the format in which to export/dump data of matrices and vectors in the file', 'solver.html#dumpformat')
]

solver_linear = SettingsSolver(
    solver_common
)

solver_nonlinear = SettingsSolver(
    solver_common + [
        SettingsDictEntry("snesMaxFunctionEvaluations", '1e3', 'maximum number of function iterations', 'hyperelasticity.html#python-settings'),
        SettingsDictEntry("snesMaxIterations", '50', 'maximum number of iterations in the nonlinear solver', 'hyperelasticity.html#python-settings'),
        SettingsDictEntry("snesRelativeTolerance", '1e-10', 'relative tolerance of the nonlinear solver', 'hyperelasticity.html#python-settings'),
        SettingsDictEntry("snesLineSearchType", '"l2"', 'type of linesearch, possible values: "bt" "nleqerr" "basic" "l2" "cp" "ncglinear"', 'hyperelasticity.html#python-settings'),
        SettingsDictEntry("snesAbsoluteTolerance", '1e-10', 'absolute tolerance of the nonlinear solver', 'hyperelasticity.html#python-settings'),
        SettingsDictEntry("snesRebuildJacobianFrequency", '5', 'how often the jacobian should be recomputed, -1 indicates NEVER rebuild, 1 means rebuild every time the Jacobian is computed within a single nonlinear solve, 2 means every second time the Jacobian is built etc. -2 means rebuild at next chance but then never again', 'hyperelasticity.html#python-settings'),
    ]
)

outputwriter = SettingsDictEntry("OutputWriter", SettingsList([
    SettingsListEntry(
        SettingsDict([
            SettingsDictEntry("format", '"Paraview"', 'one of Paraview, PythonFile, ExFile, MegaMol, PythonCallback', 'output_writer.html#outputwriter'),
            SettingsDictEntry("filename", '"out/filename"', 'the file name of the output file to write', 'output_writer.html#filename'),
            SettingsDictEntry("outputInterval", '1', 'the interval in which timesteps an actual file should be written', 'output_writer.html#outputinterval'),
            SettingsDictEntry("fileNumbering", '"incremental"', 'incremental or timeStepIndex', 'output_writer.html#filenumbering'),
            SettingsDictEntry("binary", 'True', 'whether to produce binary data files', 'output_writer.html#binary'),
            SettingsDictEntry("fixedFormat", 'True', None, 'output_writer.html#fixedformat'),
            SettingsDictEntry("combineFiles", 'False', None, 'output_writer.html#combinefiles'),
            SettingsChoice([],[
                SettingsDictEntry("onlyNodalValues", 'True', None, None),
            ]),
            SettingsChoice([],[
                SettingsDictEntry("sphereSize", '"0.005*0.005*0.01"', 'ExFile: defines how spheres, used to visualize nodes, will be rendered. The format is x*y*z', 'output_writer.html#exfile'),
            ]),
            SettingsChoice([],[
                SettingsDictEntry("callback", 'callback', 'PythonCallback: python-function to call back to', 'output_writer.html#pythoncallback'),
            ]),
        ])
    )
]), 'specifies a list of output writers that can be used to output geometry field variables in various formats', 'output_writer.html#outputwriter')

timestepping_schemes_common = [
    SettingsChoice([
        SettingsDictEntry("numberTimeSteps", '10', None, 'timestepping_schemes_ode.html#endtime-numbertimesteps-and-timestepwidth')
    ], [
        SettingsDictEntry("timeStepWidth", '0.001', None, 'timestepping_schemes_ode.html#endtime-numbertimesteps-and-timestepwidth')
    ]),
    SettingsChoice([], [
        SettingsDictEntry("logTimeStepWidthAsKey", '"timestep_width"', 'the time step width of this scheme will be stored under this key in logs/log.csv', 'timestepping_schemes_ode.html#logtimestepwidthaskey-lognumbertimestepsaskey-and-durationlogkey')
    ]),
    SettingsChoice([], [
        SettingsDictEntry("logNumberTimeStepsAsKey", '"timesteps_number"', 'the number of time steps of this scheme will be stored under this key in logs/log.csv', 'timestepping_schemes_ode.html#logtimestepwidthaskey-lognumbertimestepsaskey-and-durationlogkey')
    ]),
    SettingsChoice([], [
        SettingsDictEntry("durationLogKey", '"duration"', 'the total time that has passed for the computation will be stored under this key in logs/log.csv', 'timestepping_schemes_ode.html#logtimestepwidthaskey-lognumbertimestepsaskey-and-durationlogkey')
    ]),
    SettingsDictEntry("timeStepOutputInterval", '100', 'a positive integer value that specifies the interval in which timesteps are printed to standard output', 'timestepping_schemes_ode.html#timestepoutputinterval'),
]

timestepping_schemes_ode_common = timestepping_schemes_common + [
    SettingsDictEntry("endTime", '1', 'run() method performs the simulation for t∈[0,endTime]', 'timestepping_schemes_ode.html#endtime-numbertimesteps-and-timestepwidth'),
    SettingsDictEntry("initialValues", '[]', 'list of double values to use as initial values. The solution is set to these values upon initialization', 'timestepping_schemes_ode.html#initialvalues'),
    SettingsDictEntry("inputMeshIsGlobal", 'True', 'the degrees of freedom are interpreted in global numbering, if inputMeshIsGlobal is set to True, or in local numbering of the process, if inputMeshIsGlobal is False', 'timestepping_schemes_ode.html#dirichletboundaryconditions-and-inputmeshisglobal'),
    SettingsDictEntry("dirichletBoundaryConditions", '{}', 'dictionary with degrees of freedom as key and the value as value (i.e. {"dof": value, ...}', 'timestepping_schemes_ode.html#dirichletboundaryconditions-and-inputmeshisglobal'),
    SettingsDictEntry("dirichletOutputFilename", 'None', 'write Dirichlet Boundary conditions to .vtp file', 'boundary_conditions.html#dirichlet-output-filename'),
    SettingsDictEntry("checkForNanInf", 'False', 'check if the solution vector contains nan or +/-inf values, if yes, an error is printed. This is a time-consuming check'),
    SettingsChoice([], [
        outputwriter
    ]),
    SettingsChoice([], [
        SettingsDictEntry("nAdditionalFieldVariables", '0', 'number of additional field variables that will be created', 'timestepping_schemes_ode.html#nadditionalfieldvariables'),
        SettingsDictEntry("additionalSlotNames", '["connector_slot_1"]', 'list of strings, names for of connector slots for the additional field variables', 'timestepping_schemes_ode.html#additionalslotnames')
    ])
]

operator_splitting_common = timestepping_schemes_ode_common + [
    SettingsDictEntry("connectedSlotsTerm1To2", '[0]', 'list of slots of term 2 that are connected to the slots of term 1', 'output_connector_slots.html#connectedslotsterm1to2-and-connectedslotsterm2to1'),
    SettingsDictEntry("connectedSlotsTerm2To1", '[0]', 'list of slots of term 1 that are connected to the slots of term 2', 'output_connector_slots.html#connectedslotsterm1to2-and-connectedslotsterm2to1'),
    SettingsDictEntry("Term1", SettingsDict([
        SettingsChildPlaceholder(0)
    ])),
    SettingsDictEntry("Term2", SettingsDict([
        SettingsChildPlaceholder(1)
    ]))
]

multidomain_solver_common = timestepping_schemes_ode_common + [
    SettingsDictEntry("nCompartments", '1', 'number of compartments', 'multidomain_solver.html#python-settings'),
    SettingsDictEntry("compartmentRelativeFactors", '[]', 'list of lists of (the factors for all dofs), if "inputIsGlobal": True, this contains the global dofs', 'multidomain_solver.html#python-settings'),
    SettingsDictEntry("inputIsGlobal", 'True', 'if values and dofs correspond to the global numbering', 'multidomain_solver.html#python-settings'),
    SettingsDictEntry("am", '500.0', 'am parameter for every motor unit (ration of surface to volume of fibers)', 'multidomain_solver.html#python-settings'),
    SettingsDictEntry("cm", '0.58', 'cm parameter for every motor unit (capacitance of the cellular membrane)', 'multidomain_solver.html#python-settings'),
    # TODO maybe add a special SettingsSolver()
    SettingsDictEntry("solverName", '"multidomainLinearSolver"', 'reference to the solver used for the global linear system of the multidomain eq.', 'multidomain_solver.html#python-settings'),
    SettingsDictEntry("alternativeSolverName", '"multidomainAlternativeLinearSolver"', 'reference to the alternative solver, which is used when the normal solver diverges', 'multidomain_solver.html#python-settings'),
    SettingsDictEntry("subSolverType", '"gamg"', 'sub solver when block jacobi preconditioner is used', 'multidomain_solver.html#python-settings'),
    SettingsDictEntry("subPreconditionerType", '"none"', 'sub preconditioner when block jacobi preconditioner is used', 'multidomain_solver.html#python-settings'),
    SettingsDictEntry("gamgType", '"classical"', 'one of agg, geo, or classical', 'multidomain_solver.html#python-settings'),
    SettingsDictEntry("cycleType", '"cycleV"', 'either cycleV or cycleW', 'multidomain_solver.html#python-settings'),
    SettingsDictEntry("nLevels", '25', None, 'multidomain_solver.html#python-settings'),
    SettingsDictEntry("hypreOptions", '"-pc_hypre_boomeramg_strong_threshold 0.7"', 'additional options if a hypre preconditioner is selected', 'multidomain_solver.html#python-settings'),
    SettingsDictEntry("theta", '0.5', 'weighting factor of implicit term in Crank-Nicolson scheme, 0.5 gives the classic, 2nd-order Crank-Nicolson scheme, 1.0 gives implicit euler', 'multidomain_solver.html#python-settings'),
    SettingsDictEntry("useLumpedMassMatrix", 'True', 'which formulation to use, the formulation with lumped mass matrix (True) is more stable but approximative, the other formulation (False) is exact but needs more iterations', 'multidomain_solver.html#python-settings'),
    SettingsDictEntry("useSymmetricPreconditionerMatrix", 'True', 'if the diagonal blocks of the system matrix should be used as preconditioner matrix', 'multidomain_solver.html#python-settings'),
    SettingsDictEntry("initialGuessNonzero", 'True', 'if the initial guess for the 3D system should be set as the solution of the previous timestep, this only makes sense for iterative solvers', 'multidomain_solver.html#python-settings'),
    SettingsDictEntry("enableFatComputation", 'True', 'disabling the computation of the fat layer is only for debugging and speeds up computation. If set to False, the respective matrix is set to the identity', 'multidomain_solver.html#python-settings'),
    SettingsDictEntry("showLinearSolverOutput", 'True', 'if convergence information of the linear solver in every timestep should be printed, this is a lot of output for fast computations', 'multidomain_solver.html#python-settings'),
    SettingsDictEntry("updateSystemMatrixEveryTimestep", 'False', 'if this multidomain solver will update the system matrix in every first timestep, us this only if the geometry changed, e.g. by contraction', 'multidomain_solver.html#python-settings'),
    SettingsDictEntry("recreateLinearSolverInterval", '0', 'how often the Petsc KSP object (linear solver) should be deleted and recreated. This is to remedy memory leaks in Petsc\'s implementation of some solvers. 0 means disabled.', 'multidomain_solver.html#python-settings')
]

hyperelasticity_common = [
    SettingsDictEntry("materialParameters", '[]', 'list of material parameters, must match the number of parameters in the material', 'hyperelasticity.html#materialparameters'),
    SettingsDictEntry("density", '1.0', 'density of the material', 'hyperelasticity.html#python-settings'),
    SettingsDictEntry("displacementsScalingFactor", '1.0', 'scaling factor for displacements, only set to sth. other than 1 only to increase visual appearance for very small displacements', 'hyperelasticity.html#python-settings'),
    SettingsDictEntry("residualNormLogFilename", '"residual_norm.txt"', 'log file where residual norm values of the nonlinear solver will be written', 'hyperelasticity.html#python-settings'),
    SettingsDictEntry("useAnalyticJacobian", 'True', 'whether to use the analytically computed jacobian matrix in the nonlinear solver (fast)', 'hyperelasticity.html#python-settings'),
    SettingsDictEntry("useNumericJacobian", 'True', 'whether to use the numerically computed jacobian matrix in the nonlinear solver (slow), only works with non-nested matrices, if both numeric and analytic are enable, it uses the analytic for the preconditioner and the numeric as normal jacobian', 'hyperelasticity.html#python-settings'),
    # undocumented
    SettingsDictEntry("nNonlinearSolveCalls", '1', 'how often the nonlinear solve should be called'),
    # undocumented
    SettingsDictEntry("loadFactorGiveUpThreshold", '1e-5', 'a threshold for the load factor, when to abort the solve of the current time step. The load factors are adjusted automatically if the nonlinear solver diverged. If the load factors get too small, it aborts the solve'),
    SettingsDictEntry("dumpDenseMatlabVariables", 'False', 'whether to have extra output of matlab vectors, x,r, jacobian matrix (very slow)', 'hyperelasticity.html#python-settings'),
    # we use meshName directly here, because the Mesh::CompositeOfDimension is sometimes CHILD0 and sometimes CHILD1
    SettingsDictEntry("meshName", '"meshX"', 'mesh with quadratic Lagrange ansatz functions', 'hyperelasticity.html#python-settings'),
    SettingsDictEntry("inputMeshIsGlobal", 'True', 'boundary conditions are specified in global numberings, whereas the mesh is given in local numberings', 'hyperelasticity.html#python-settings'),
    SettingsChoice([
        SettingsDictEntry("fiberMeshNames", '[]', 'fiber meshes that will be used to determine the fiber direction', 'hyperelasticity.html#python-settings')
    ], [
        SettingsDictEntry("fiberDirection", '[0, 0, 1]', 'if fiberMeshNames is empty, directly set the constant fiber direction, in element coordinate system', 'hyperelasticity.html#python-settings')
    ]),
    SettingsDictEntry("loadFactors", '[]', 'if []: no load factors are used', 'hyperelasticity.html#python-settings'),
    SettingsDictEntry("nNonlinearSolveCalls", '1', 'how often the nonlinear solve should be called', 'hyperelasticity.html#python-settings'),
    SettingsDictEntry("dirichletBoundaryConditions", '{}', 'the initial Dirichlet boundary conditions that define values for displacements u', 'hyperelasticity.html#python-settings'),
    SettingsDictEntry("dirichletOutputFilename", '""', 'a filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable. This is for debugging the Dirichlet boundary condition nodes', 'hyperelasticity.html#dirichletoutputfilename'),
    SettingsDictEntry("neumannBoundaryConditions", '[]', 'neumann boundary conditions that define traction forces on surfaces of elements', 'hyperelasticity.html#python-settings'),
    SettingsDictEntry("divideNeumannBoundaryConditionValuesByTotalArea", 'False', 'if the given Neumann boundary condition values under "neumannBoundaryConditions" are total forces instead of surface loads and therefore should be scaled by the surface area of all elements where Neumann BC are applied', 'hyperelasticity.html#python-settings'),
    SettingsDictEntry("updateDirichletBoundaryConditionsFunction", 'None', 'function that updates the dirichlet BCs while the simulation is running', 'hyperelasticity.html#python-settings'),
    SettingsDictEntry("updateDirichletBoundaryConditionsFunctionCallInterval", '1', 'every which step the update function should be called, 1 means every time step', 'hyperelasticity.html#python-settings'),
    SettingsDictEntry("initialValuesDisplacements", '[[0.0,0.0,0.0] for _ in range(mx*my*mz)]', 'initial values for the displacements, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]', 'hyperelasticity.html#python-settings'),
    SettingsDictEntry("initialValuesVelocities", '[[0.0,0.0,0.0] for _ in range(mx*my*mz)]', 'initial values for the velocities, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]', 'hyperelasticity.html#python-settings'),
    SettingsDictEntry("extrapolateInitialGuess", 'True', 'if the initial values for the dynamic nonlinear problem should be computed by extrapolating the previous displacements and velocities', 'hyperelasticity.html#python-settings'),
    SettingsDictEntry("constantBodyForce", '[0.0, 0.0, 0.0]', 'a constant force that acts on the whole body, e.g. for gravity', 'hyperelasticity.html#python-settings'),
    SettingsChoice([],[
        outputwriter
    ]),
    SettingsChoice([],[
        SettingsDictEntry("pressure", SettingsDict([outputwriter]), 'additional output writer that writes also the hydrostatic pressure', 'hyperelasticity.html#python-settings')
    ]),
    SettingsChoice([],[
        SettingsDictEntry("LoadIncrements", SettingsDict([outputwriter]), 'output writer for debugging, outputs files after each load increment, the geometry is not changed but u and v are written', 'hyperelasticity.html#python-settings')
    ]),
]

# this is used in HyperelasticitySolver and MuscleContractionSolver
hyperelasticity_solver = SettingsDictEntry("HyperelasticitySolver", SettingsDict(
        [SettingsDictEntry("durationLogKey", '"duration_mechanics"', 'key to find duration of this solver in the log file', 'hyperelasticity.html#python-settings'),] +
        hyperelasticity_common +
        [solver_nonlinear]
    ))

# this is used in DynamicHyperelasticitySolver and MuscleContractionSolver
dynamic_hyperelasticity_solver = SettingsDictEntry("DynamicHyperelasticitySolver", SettingsDict(
        timestepping_schemes_common + [
            SettingsDictEntry("endTime", '1', 'run() method performs the simulation for t∈[0,endTime]', 'timestepping_schemes_ode.html#endtime-numbertimesteps-and-timestepwidth'),
        ] +
        hyperelasticity_common + [
            SettingsDictEntry("updateNeumannBoundaryConditionsFunction", 'None', 'function that updates the Neumann BCs while the simulation is running', 'dynamic_hyperelasticity.html#python-settings'),
            SettingsDictEntry("updateNeumannBoundaryConditionsFunctionCallInterval", '1', 'every which step the update function should be called, 1 means every time step', 'dynamic_hyperelasticity.html#python-settings'),
            SettingsChoice([],[
                SettingsDictEntry("dynamic", SettingsDict([outputwriter]), 'additional output writer that writes virtual work terms', 'dynamic_hyperelasticity.html#python-settings')
            ]),
        ] +
        [solver_nonlinear]
    ))


possible_solver_combinations = {
    "GLOBAL": {
        # template_arguments gets set by cpp_tree.py (to all runnables)
        "python_options": SettingsDict([
            SettingsDictEntry("scenarioName", '"test-scenario"'),
            SettingsDictEntry("logFormat", '"csv"', "csv or json"),
            SettingsDictEntry("meta", SettingsDict(), 'additional fields that will appear in the log'),
            SettingsDictEntry("solverStructureDiagramFile",
                              '"solver_structure.txt"', 'filename of file that will contain a visualization of the solver structure and data mapping', 'output_connector_slots.html#solverstructurediagramfile'),
            SettingsDictEntry("connectedSlots", '[]', None, 'output_connector_slots.html#using-global-slot-names'),
            SettingsDictEntry("mappingsBetweenMeshesLogFile",
                              '""', 'this is the name of a log file that will contain events during creation and mapping', 'mappings_between_meshes.html#mappingsbetweenmesheslogfile'),
            SettingsDictEntry("MappingsBetweenMeshes", '{}', None, 'mappings_between_meshes.html#mappingsbetweenmeshes'),
            SettingsChildPlaceholder(0)
        ])
    },

    "Dummy": {
        "timeSteppingScheme" : True
    },

    "Postprocessing::ParallelFiberEstimation": {
        "runnable": True,
        "template_arguments": [
            ('BasisFunction', ["BasisFunction::"])
        ]
    },
    "Postprocessing::StreamlineTracer": {
        "runnable": True,
        "template_arguments": [
            ('DiscretizableInTime', ["discretizableInTime"])
        ]
    },


    "PreciceAdapter::ContractionDirichletBoundaryConditions": {
        "runnable": True,
        "template_arguments": [
            ('TimeSteppingScheme', ["timeSteppingScheme"])
        ]
    },
    "PreciceAdapter::ContractionNeumannBoundaryConditions": {
        "runnable": True,
        "template_arguments": [
            ('TimeSteppingScheme', ["timeSteppingScheme"])
        ]
    },
    "PreciceAdapter::PartitionedFibers": {
        "runnable": True,
        "template_arguments": [
            ('TimeSteppingScheme', ["timeSteppingScheme"])
        ]
    },
    "PreciceAdapter::MuscleContraction": {
        "runnable": True,
        "template_arguments": [
            ('MuscleContractionSolver', ["MuscleContractionSolver"])
        ]
    },
    "MuscleContractionSolver": {
        "runnable": True,
        "timeSteppingScheme": True,
        "template_arguments_needed": 0,
        "template_arguments": [
            ('Mesh', ["Mesh::"]),
            ('Material', [
             "Equation::"]),
        ],
        "python_options": SettingsDict([
            SettingsDictEntry("MuscleContractionSolver", SettingsDict([
                SettingsDictEntry("numberTimeSteps", '1', 'number of timesteps to use per call', 'muscle_contraction_solver.html#python-settings'),
                SettingsDictEntry("timeStepOutputInterval", '100', 'how often the current timestep will be displayed, if this is >100 and numberTimeSteps is 1, nothing will be printed', 'muscle_contraction_solver.html#python-settings'),
                SettingsDictEntry("Pmax", '1.0', 'maximum PK2 active stress', 'muscle_contraction_solver.html#python-settings'),
                SettingsChoice([],[
                    outputwriter
                ]),
                SettingsDictEntry("mapGeometryToMeshes", '[]', 'the mesh names of the meshes that will get the geometry transferred', 'muscle_contraction_solver.html#python-settings'),
                SettingsDictEntry("slotNames", '[]', 'names of the connector slots, maximum 6 characters per name', 'muscle_contraction_solver.html#python-settings'),
                SettingsDictEntry("dynamic", 'True', 'if the dynamic solid mechanics solver should be used, else it computes the quasi-static problem', 'muscle_contraction_solver.html#python-settings'),
                SettingsChoice([
                    dynamic_hyperelasticity_solver
                ],[
                    hyperelasticity_solver
                ]),
            ]))
        ])
    },
    "FastMonodomainSolver": {
        "runnable": True,
        "timeSteppingScheme": True,
        "template_arguments": [
            ('Nested solver', ["Control::MultipleInstances"])
        ],
        "python_options" : SettingsDict([
            SettingsChildPlaceholder(0),
            SettingsDictEntry("fiberDistributionFile", '"MU_fibre_distribution_3780.txt"', 'this file contains the assignment of fibers to motor units', 'fast_monodomain_solver.html#fiberdistributionfile'),
            SettingsDictEntry("firingTimesFile", '"MU_firing_times_real.txt"', 'this file specifies when which motor unit fires', 'fast_monodomain_solver.html#firingtimesfile'),
            SettingsDictEntry("onlyComputeIfHasBeenStimulated", 'True', 'if True: disable computation of the Monodomain equation as long as the fiber has not been stimulated in therefore is in equilibrium', 'fast_monodomain_solver.html#onlycomputeifhasbeenstimulated'),
            SettingsDictEntry("disableComputationWhenStatesAreCloseToEquilibrium", 'True', 'similar to onlyComputeIfHasBeenStimulated, this checks whether the values have reached the equilibrium and then disables the computation', 'fast_monodomain_solver.html#disablecomputationwhenstatesareclosetoequilibrium'),
            SettingsDictEntry("valueForStimulatedPoint", '20.0', 'value that will be set for the transmembrane potential Vm when it is stimulated', 'fast_monodomain_solver.html#valueforstimulatedpoint'),
            SettingsDictEntry("neuromuscularJunctionRelativeSize", '0.0', 'relative range of the position of the neuromuscular junction', 'fast_monodomain_solver.html#neuromuscularjunctionrelativesize'),
        ])
    },
    "SpatialDiscretization::HyperelasticitySolver": {
        "runnable": True,
        # TODO can this be handled like a timeSteppingScheme?
        "timeSteppingScheme": True,
        "template_arguments_needed": 0,
        "template_arguments": [
            ('Material', [
             "Equation::SolidMechanics::"]),
            ('Mesh', ["Mesh::StructuredRegularFixedOfDimension"]),
            ('nDisplacementComponents', ["Integer"])
        ],
        "python_options": SettingsDict([
            hyperelasticity_solver
        ])
    },


    "Control::MultipleInstances": {
        "runnable": True,
        "timeSteppingScheme": True,
        "template_arguments": [
            ('TimeSteppingScheme', ["timeSteppingScheme"])
        ],
        "python_options": SettingsDict([
            SettingsDictEntry("MultipleInstances", SettingsDict([
                SettingsDictEntry("nInstances", '1', 'the number of instance to create and run in total', 'multiple_instances.html#ninstances'),
                SettingsDictEntry("instances", SettingsList([
                    SettingsListEntry(SettingsDict([
                        SettingsDictEntry("ranks", '[0]', 'list of MPI ranks on which the instance should be computed', 'multiple_instances.html#ranks'),
                        SettingsChildPlaceholder(0)
                    ]))
                ]), 'list of settings for the instances', 'multiple_instances.html#instances'),
                SettingsChoice([],[
                    outputwriter
                ])
            ]))
        ])
    },


    "Control::LoadBalancing": {
        "runnable": True,
        "timeSteppingScheme": True,
        "template_arguments": [
            ('TimeSteppingScheme', ["timeSteppingScheme"])
        ]
    },
    "Control::MapDofs": {
        "runnable": True,
        "timeSteppingScheme": True,
        "template_arguments": [
            ('FunctionSpace', ["FunctionSpace::"]),
            ('TimeSteppingScheme', ["timeSteppingScheme"])
        ],
        "python_options": SettingsDict([
            SettingsDictEntry("MapDofs", SettingsDict([
                # mesh from FunctionSpace
                SettingsChildPlaceholder(0),
                # nested solver
                SettingsChildPlaceholder(1),
                SettingsDictEntry("nAdditionalFieldVariables", '0', 'number of additional field variables that are defined by this object. They have 1 component, use the templated function space and mesh given by meshName', 'map_dofs.html#nadditionalfieldvariables'),
                SettingsDictEntry("additionalslotnames", '[]', 'list of names of the slots for the additional field variables', 'map_dofs.html#additionalslotnames'),
                # TODO maybe define options in beforeComputation and afterComputation
                SettingsDictEntry("beforeComputation", 'None', 'transfer/mapping of dofs that will be performed before the computation of the nested solver, can be None if not needed', 'map_dofs.html#beforecomputation-and-aftercomputation'),
                SettingsDictEntry("afterComputation", 'None', 'transfer/mapping of dofs that will be performed after the computation of the nested solver, can be None if not needed', 'map_dofs.html#beforecomputation-and-aftercomputation'),
            ]))
        ])
    },


    "OperatorSplitting::Godunov": {
        "runnable": True,
        "timeSteppingScheme": True,
        "template_arguments": [
            ('First TimeSteppingScheme', ["timeSteppingScheme"]),
            ('Second TimeSteppingScheme', ["timeSteppingScheme"])
        ],
        "python_options": SettingsDict([
            SettingsDictEntry("GodunovSplitting", SettingsDict(
                operator_splitting_common
            ))
        ])
    },
    "OperatorSplitting::Strang": {
        "runnable": True,
        "timeSteppingScheme": True,
        "template_arguments": [
            ('First TimeSteppingScheme', ["timeSteppingScheme"]),
            ('Second TimeSteppingScheme', ["timeSteppingScheme"])
        ],
        "python_options": SettingsDict([
            SettingsDictEntry("StrangSplitting", SettingsDict(
                operator_splitting_common
            ))
        ])
    },


    "Control::Coupling": {
        "runnable": True,
        "timeSteppingScheme": True,
        "template_arguments": [
            ('First TimeSteppingScheme', ["timeSteppingScheme"]),
            ('Second TimeSteppingScheme', ["timeSteppingScheme"])
        ],
        "python_options": SettingsDict([
            SettingsDictEntry("Coupling", SettingsDict(
                operator_splitting_common
            ))
        ])
    },
    "Control::MultipleCoupling": {
        "runnable": True,
        "timeSteppingScheme": True,
        "template_arguments_needed" : 2,
        "template_arguments": [
            ('First TimeSteppingScheme', ["timeSteppingScheme"]),
            ('Second TimeSteppingScheme', ["timeSteppingScheme"]),
            ('Third TimeSteppingScheme', ["timeSteppingScheme"]),
            ('Fourth TimeSteppingScheme', ["timeSteppingScheme"]),
            ('Fifth TimeSteppingScheme', ["timeSteppingScheme"])
        ],
        "python_options": SettingsDict([
            SettingsDictEntry("MultipleCoupling", SettingsDict(
                timestepping_schemes_ode_common + [
                    SettingsDictEntry("connectedSlotsTerm1To2", 'None', 'this would be the connected slots for a normal Coupling scheme, but here it should be set to None, use global option "connectedSlots" instead', 'output_connector_slots.html#connectedslotsterm1to2-and-connectedslotsterm2to1'),
                    SettingsDictEntry("connectedSlotsTerm2To1", 'None', 'this would be the connected slots for a normal Coupling scheme, but here it should be set to None, use global option "connectedSlots" instead', 'output_connector_slots.html#connectedslotsterm1to2-and-connectedslotsterm2to1'),
                    SettingsDictEntry("Term1", SettingsDict([
                        SettingsChildPlaceholder(0)
                    ])),
                    SettingsDictEntry("Term2", SettingsDict([
                        SettingsChildPlaceholder(1)
                    ])),
                    SettingsDictEntry("Term3", SettingsDict([
                        SettingsChildPlaceholder(2)
                    ])),
                    SettingsDictEntry("Term4", SettingsDict([
                        SettingsChildPlaceholder(3)
                    ])),
                    SettingsDictEntry("Term5", SettingsDict([
                        SettingsChildPlaceholder(4)
                    ]))
                ]
            ))
        ])
    },


    "CellmlAdapter": {
        "discretizableInTime": True,
        "template_arguments_needed": 2,
        "template_arguments": [
            ('Number of states', ["Integer"]),
            ('Number of algebraics', ["Integer"]),
            ('FunctionSpace', ["FunctionSpace::"])
        ],
        "python_options": SettingsDict([
            SettingsDictEntry("CellML", SettingsDict([
                SettingsChoice([
                    SettingsDictEntry("modelFilename", '"../../input/hodgkin_huxley_1952.c"', 'filename of the CellML model file (cellml-XML or c/c++)', 'cellml_adapter.html#modelfilename')
                ], [
                    SettingsDictEntry("libraryFilename", "../lib/model.so", 'filename of a shared object library (.so) that will be used to compute the model', 'cellml_adapter.html#libraryfilename')
                ]),
                SettingsDictEntry("statesInitialValues", '"CellML"', 'list with initial values for each state of the CellML model or "CellML" or "undefined"', 'cellml_adapter.html#statesinitialvalues'),
                SettingsDictEntry("initializeStatesToEquilibrium", 'True', 'if the equilibrium values of the states should be computed before the simulation starts', 'cellml_adapter.html#initializestatestoequilibrium-and-initializestatestoequilibriumtimestepwidth'),
                SettingsDictEntry("initializeStatesToEquilibriumTimestepWidth", '1e-4', 'timestep width to use to solve the equilibrium equation', 'cellml_adapter.html#initializestatestoequilibrium-and-initializestatestoequilibriumtimestepwidth'),
                SettingsDictEntry("additionalArgument", 'None', 'each callback function will get this as its last argument', 'cellml_adapter.html#callbacks'),
                SettingsChoice([],[
                    SettingsDictEntry("setSpecificParametersFunction", 'set_specific_parameters', 'function name', 'cellml_adapter.html#setspecificparametersfunction-and-setspecificparameterscallinterval')
                ]),
                SettingsChoice([],[
                    SettingsDictEntry("setSpecificParametersCallInterval", '0', None, 'cellml_adapter.html#setspecificparametersfunction-and-setspecificparameterscallinterval')
                ]),
                SettingsChoice([],[
                    SettingsDictEntry("setSpecificStatesFunction", 'set_specific_states', 'function name', 'cellml_adapter.html#setspecificstatesfunction-and-setspecificstatescallinterval')
                ]),
                SettingsChoice([],[
                    SettingsDictEntry("setSpecificStatesCallInterval", '1', None, 'cellml_adapter.html#setspecificstatesfunction-and-setspecificstatescallinterval')
                ]),
                SettingsChoice([],[
                    SettingsDictEntry("setSpecificStatesCallEnableBegin", '0', None, 'cellml_adapter.html#setspecificstatescallenablebegin-setspecificstatescallfrequency-and-setspecificstatesfrequencyjitter'),
                    SettingsDictEntry("setSpecificStatesCallFrequency", '0.1', None, 'cellml_adapter.html#setspecificstatescallenablebegin-setspecificstatescallfrequency-and-setspecificstatesfrequencyjitter'),
                    SettingsDictEntry("setSpecificStatesFrequencyJitter", '0', None, 'cellml_adapter.html#setspecificstatescallenablebegin-setspecificstatescallfrequency-and-setspecificstatesfrequencyjitter'),
                    SettingsDictEntry("setSpecificStatesRepeatAfterFirstCall", '0.0', '[ms] simulation time span for which the setSpecificStates callback will be called after a call was triggered')
                ]),
                SettingsChoice([], [
                    SettingsDictEntry("handleResultFunction", 'handle_result', 'function name', 'cellml_adapter.html#handleresultfunction-and-handleresultcallinterval')
                ]),
                SettingsChoice([], [
                    SettingsDictEntry("handleResultCallInterval", '1', None, 'cellml_adapter.html#handleresultfunction-and-handleresultcallinterval')
                ]),
                SettingsChoice([], [
                    SettingsDictEntry("parametersUsedAsAlgebraic", '[]', 'list of algebraic numbers that will be replaced by parameters', 'cellml_adapter.html#parametersusedasalgebraic')
                ]),
                SettingsChoice([], [
                    SettingsDictEntry("parametersUsedAsConstant", '[]', 'list of indices, which constants in the computation will be replaced by parameters', 'cellml_adapter.html#parametersusedasconstant')
                ]),
                SettingsChoice([], [
                    SettingsDictEntry("algebraicsForTransfer", '[]', 'list of algebraics that should be transferred to the other solver in either a Coupling, GodunovSplitting or StrangSplitting', 'cellml_adapter.html#algebraicsfortransfer-and-statesfortransfer')
                ]),
                SettingsChoice([], [
                    SettingsDictEntry("statesForTransfer", '[]', 'list of states that should be transferred to the other solver in either a Coupling, GodunovSplitting or StrangSplitting', 'cellml_adapter.html#algebraicsfortransfer-and-statesfortransfer')
                ]),
                SettingsChoice([], [
                    SettingsDictEntry("parametersInitialValues", '[]', 'list of values of the parameters. This also defines the number of parameters', 'cellml_adapter.html#parametersinitialvalues')
                ]),
                SettingsChoice([], [
                    SettingsDictEntry("mappings", '{}', None, 'cellml_adapter.html#mappings')
                ]),
                SettingsChoice([
                    SettingsDictEntry("nElements", '0', 'set the number of instances to be computed', 'cellml_adapter.html#meshname')
                ],[
                    SettingsDictEntry("meshName", '"meshX"', 'the mesh to use, to be defined under "Meshes"', 'cellml_adapter.html#meshname')
                ]),
                SettingsDictEntry("stimulationLogFilename", '"out/stimulation.log"', 'file name of an output file that will contain all firing times', 'cellml_adapter.html#stimulationlogfilename'),
                SettingsDictEntry("optimizationType", '"vc"', 'one of simd, vc, openmp', 'cellml_adapter.html#optimizationtype'),
                SettingsDictEntry("approximateExponentialFunction", 'True', 'if optimizationType is "vc", whether the exponential function exp(x) should be approximate by (1+x/n)^n with n=1024', 'cellml_adapter.html#optimizationtype'),
                SettingsDictEntry("maximumNumberOfThreads", '0', 'if optimizationType is "openmp", the maximum number of threads to use. Default value 0 means no restriction', 'cellml_adapter.html#optimizationtype'),
                SettingsChoice([], [
                    SettingsDictEntry("compilerFlags", '"-fPIC -finstrument-functions -ftree-vectorize -fopt-info-vec-optimized=vectorizer_optimized.log -shared"', 'additional compiler flags for the compilation of the source file', 'cellml_adapter.html#compilerflags')
                ]),
                SettingsChoice([], [
                    outputwriter
                ])
            ]))
        ])
    },


    "ModelOrderReduction::POD": {
        "discretizableInTime": True,
        "template_arguments": [
            ('DiscretizableInTime', ["discretizableInTime"]),
            ('LinearPart', ["ModelOrderReduction::LinearPart"])
        ]
    },
    "ModelOrderReduction::LinearPart": {},
    "ModelOrderReduction::ExplicitEulerReduced": {
        "timeSteppingScheme": True,
        "template_arguments": [
            ('ExplicitEuler', ["TimeSteppingScheme::ExplicitEuler"])
        ]
    },
    "ModelOrderReduction::ImplicitEulerReduced": {
        "timeSteppingScheme": True,
        "template_arguments": [
            ('ImplicitEuler', ["TimeSteppingScheme::ImplicitEuler"])
        ]
    },


    "FunctionSpace::FunctionSpace": {
        "template_arguments": [
            ('Mesh', ["Mesh::"]),
            ('BasisFunction', ["BasisFunction::"])
        ],
        "python_options": SettingsDict([
            # just proxy the mesh (this is used in PrescribedValues)
            SettingsChildPlaceholder(0)
        ])
    },


    "TimeSteppingScheme::ExplicitEuler": {
        "runnable": True,
        "timeSteppingScheme": True,
        "template_arguments": [
            ('DiscretizableInTime', ["discretizableInTime"])
        ],
        "python_options": SettingsDict([
            SettingsDictEntry("ExplicitEuler", SettingsDict(
                timestepping_schemes_ode_common + [SettingsChildPlaceholder(0)]
            ))
        ])
    },
    "TimeSteppingScheme::ImplicitEuler": {
        "runnable": True,
        "timeSteppingScheme": True,
        "template_arguments": [
            ('DiscretizableInTime', ["discretizableInTime"])
        ],
        "python_options": SettingsDict([
            SettingsDictEntry("ImplicitEuler", SettingsDict(
                timestepping_schemes_ode_common + [
                    SettingsChildPlaceholder(0),
                    solver_linear,
                    SettingsDictEntry("timeStepWidthRelativeTolerance", '1e-10', 'tolerance for the time step width which controls when the system matrix has to be recomputed', 'timestepping_schemes_ode.html#impliciteuler'),
                    SettingsChoice([], [
                        SettingsDictEntry("timeStepWidthRelativeToleranceAsKey", '"relative_tolerance"', 'timeStepWidthRelativeTolerance will be stored under this key in logs/log.csv', 'timestepping_schemes_ode.html#impliciteuler')
                    ]),
                    SettingsChoice([], [
                        SettingsDictEntry("durationInitTimeStepLogKey", '"duration_init_time_step"', 'duration for the time step initialization  will be stored under this key in logs/log.csv', 'timestepping_schemes_ode.html#impliciteuler')
                    ])
                ]
            ))
        ])
    },
    "TimeSteppingScheme::Heun": {
        "runnable": True,
        "timeSteppingScheme": True,
        "template_arguments": [
            ('DiscretizableInTime', ["discretizableInTime"])
        ],
        "python_options": SettingsDict([
            SettingsDictEntry("Heun", SettingsDict(
                timestepping_schemes_ode_common + [SettingsChildPlaceholder(0)]
            ))
        ])
    },
    "TimeSteppingScheme::HeunAdaptive": {
        "runnable": True,
        "timeSteppingScheme": True,
        "template_arguments": [
            ('DiscretizableInTime', ["discretizableInTime"])
        ],
        "python_options": SettingsDict([
            SettingsDictEntry("HeunAdaptive", SettingsDict(
                timestepping_schemes_ode_common + [
                    SettingsChildPlaceholder(0),
                    SettingsDictEntry("tolerance", '0.1', 'tolerance for the estimated error. It is guaranteed, that the error is always smaller than this value', 'timestepping_schemes_ode.html#tolerance'),
                    SettingsDictEntry("minTimeStepWidth", '1e-6', 'the minimum timestepwidth to use', 'timestepping_schemes_ode.html#mintimestepwidth'),
                    SettingsDictEntry("timeStepAdaptOption", '"regular"', 'method for the adaptive time step width computation (regular or modified)', 'timestepping_schemes_ode.html#timestepadaptoption'),
                    SettingsChoice([], [
                        SettingsDictEntry("lowestMultiplier", '1000', 'minimum number of timesteps to perform in the time span for the "modified" method', 'timestepping_schemes_ode.html#lowestmultiplier')
                    ])
                ]
            ))
        ])
    },
    "TimeSteppingScheme::CrankNicolson": {
        "timeSteppingScheme": True,
        "template_arguments": [
            ('DiscretizableInTime', ["discretizableInTime"])
        ],
        "python_options": SettingsDict([
            SettingsDictEntry("CrankNicolson", SettingsDict(
                timestepping_schemes_ode_common + [
                    SettingsChildPlaceholder(0),
                    solver_linear,
                    SettingsDictEntry("timeStepWidthRelativeTolerance", '1e-10', 'tolerance for the time step width which controls when the system matrix has to be recomputed', 'timestepping_schemes_ode.html#lowestmultiplier'),
                    SettingsChoice([], [
                        SettingsDictEntry("timeStepWidthRelativeToleranceAsKey", '"relative_tolerance"', 'timeStepWidthRelativeTolerance will be stored under this key in logs/log.csv', 'timestepping_schemes_ode.html#lowestmultiplier')
                    ]),
                    SettingsChoice([], [
                        SettingsDictEntry("durationInitTimeStepLogKey", '"duration_init_time_step"', 'duration for the time step initialization  will be stored under this key in logs/log.csv', 'timestepping_schemes_ode.html#lowestmultiplier')
                    ])
                ]
            ))
        ])
    },


    "TimeSteppingScheme::RepeatedCall": {
        "runnable": True,
        "timeSteppingScheme": True,
        "template_arguments": [
            ('TimeSteppingScheme', ["timeSteppingScheme"])
        ]
    },
    "TimeSteppingScheme::RepeatedCallStatic": {
        "runnable": True,
        "timeSteppingScheme": True,
        "template_arguments": [
            ('FiniteElementMethod', ["SpatialDiscretization::FiniteElementMethod"])
        ]
    },
    # specalizedSolvers:
    "TimeSteppingScheme::DynamicHyperelasticitySolver": {
        "runnable": True,
        "timeSteppingScheme": True,
        "template_arguments_needed": 0,
        "template_arguments": [
            ('Equation', ["Equation::"]),
            ('Mesh', ["Mesh::StructuredRegularFixedOfDimension"])
        ],
        "python_options": SettingsDict([
            dynamic_hyperelasticity_solver
        ])
    },
    "TimeSteppingScheme::StaticBidomainSolver": {
        "runnable": True,
        "timeSteppingScheme": True,
        "template_arguments": [
            ('FiniteElementMethod', ["SpatialDiscretization::FiniteElementMethod"]),
            ('FiniteElementMethod', ["SpatialDiscretization::FiniteElementMethod"])
        ],
        "python_options" : SettingsDict([
            SettingsDictEntry("StaticBidomainSolver", SettingsDict(
                timestepping_schemes_common + [
                SettingsDictEntry("solverName", 'activationSolver', None, 'static_bidomain_solver.html#python-settings'),
                SettingsDictEntry("initialGuessNonzero", 'True', None, 'static_bidomain_solver.html#python-settings'),
                SettingsDictEntry("slotNames", '[]', 'list of strings, names for the connector slots', 'static_bidomain_solver.html#slotnames'),
                SettingsDictEntry("PotentialFlow", SettingsDict([
                    SettingsChildPlaceholder(0)
                ])),
                SettingsDictEntry("Activation", SettingsDict([
                    SettingsChildPlaceholder(1)
                ])),
            ]))
        ])
    },
    "TimeSteppingScheme::MultidomainSolver": {
        "runnable": True,
        "timeSteppingScheme": True,
        "template_arguments": [
            ('FiniteElementMethod', ["SpatialDiscretization::FiniteElementMethod"]),
            ('FiniteElementMethod', ["SpatialDiscretization::FiniteElementMethod"])
        ],
        "python_options" : SettingsDict([
            SettingsDictEntry("MultidomainSolver", SettingsDict(
                multidomain_solver_common + [
                    SettingsDictEntry("PotentialFlow", SettingsDict([
                        SettingsChildPlaceholder(0)
                    ])),
                    SettingsDictEntry("Activation", SettingsDict([
                        SettingsChildPlaceholder(1)
                    ])),
                ]
            ))
        ])
    },
    "TimeSteppingScheme::MultidomainWithFatSolver": {
        "runnable": True,
        "timeSteppingScheme": True,
        "template_arguments": [
            ('FiniteElementMethod', ["SpatialDiscretization::FiniteElementMethod"]),
            ('FiniteElementMethod', ["SpatialDiscretization::FiniteElementMethod"]),
            ('FiniteElementMethod', ["SpatialDiscretization::FiniteElementMethod"])
        ],
        "python_options": SettingsDict([
            SettingsDictEntry("MultidomainSolver", SettingsDict(
                multidomain_solver_common + [
                    SettingsDictEntry("PotentialFlow", SettingsDict([
                        SettingsChildPlaceholder(0)
                    ])),
                    SettingsDictEntry("Activation", SettingsDict([
                        SettingsChildPlaceholder(1)
                    ])),
                    SettingsDictEntry("Fat", SettingsDict([
                        SettingsChildPlaceholder(2)
                    ]))
                ]
            ))
        ])
    },
    "TimeSteppingScheme::QuasiStaticNonlinearElasticitySolverFebio": {
        "runnable": True,
        "timeSteppingScheme": True
    },
    "TimeSteppingScheme::NonlinearElasticitySolverFebio": {
        "runnable": True,
        "timeSteppingScheme": True
    },
    "TimeSteppingScheme::QuasiStaticLinearElasticitySolver": {
        "runnable": True,
        "timeSteppingScheme": True,
        "template_arguments": [
            ('FiniteElementMethod', ["SpatialDiscretization::FiniteElementMethod"])
        ],
        "python_options" : SettingsDict([
            SettingsDictEntry("QuasiStaticLinearElasticitySolver", SettingsDict([
                SettingsChildPlaceholder(0),
                SettingsDictEntry("fiberDirection", '[0, 0, 1]', 'direction for anisotropy of elasticity formulation', 'quasi_static_linear_elasticity_solver.html#python-settings'),
                SettingsDictEntry("PotentialFlow", 'None', 'if fiberDirection is not set to a constant direction, a potential flow simulation can be used where the fiber direction is set to the streamlines of the flow through the volume. In this case, set "PotentialFlow" to the settings for the FEM for the potential flow.', 'quasi_static_linear_elasticity_solver.html#python-settings'),
                SettingsDictEntry("maximumActiveStress", '1.0', 'scaling factor to the active stress, σ_active = activation * anisotropyTensor * maximumActiveStress', 'quasi_static_linear_elasticity_solver.html#python-settings'),
                SettingsDictEntry("strainScalingCurveWidth", '1.0', 'parameter for strain-stress curve of active stress, has no effect, because strain-stress curve is commented out in the code', 'quasi_static_linear_elasticity_solver.html#python-settings'),
                SettingsDictEntry("scalingFactor", '1.0', 'scaling factor for displacements, to overrate them, if != 0 it is only for visualization purposes and not physical', 'quasi_static_linear_elasticity_solver.html#python-settings'),
                SettingsDictEntry("inputMeshIsGlobal", 'True', 'if boundary conditions are specified in global numbering', 'quasi_static_linear_elasticity_solver.html#python-settings'),
                SettingsDictEntry("slotNames", '[]', 'list of strings, names for the connector slots', 'quasi_static_linear_elasticity_solver.html#slotnames'),
                SettingsDictEntry("anisotropyTensor", 'anisotropy for active stress. The tensor is given in a local basis where the fiber direction is (1,0,0), one list item = same tensor for all elements, multiple list items = a different tensor for each element. The tensor has to be symmetric', '[1,0,0, 0,0,0, 0,0,0]', 'quasi_static_linear_elasticity_solver.html#python-settings'),
            ]))
        ])
    },
    "TimeSteppingScheme::QuasiStaticNonlinearElasticitySolverChaste": {
        "runnable": True,
        "timeSteppingScheme": True,
        "template_arguments": [
            ('Integer', ["Integer"])
        ]
    },


    "PrescribedValues": {"runnable": True,
                         "timeSteppingScheme": True,
                         "template_arguments_needed": 1,
                         "template_arguments": [
                             ('FunctionSpace', ["FunctionSpace::FunctionSpace"]),
                             ('nComponents1', ["1", "2", "3"]),
                             ('nComponents2', ["1", "2", "3"])
                         ],
        "python_options" : SettingsDict([
            SettingsDictEntry("PrescribedValues", SettingsDict([
                SettingsChildPlaceholder(0),
                #SettingsDictEntry("meshName", '"meshX"', 'the mesh to use for the field variables', 'prescribed_values.html#meshname'),
                SettingsDictEntry("numberTimeSteps", '1', 'number of timesteps to call the callback functions subsequently, this is usually 1 for prescribed values, because it is enough to set the reaction term only once per time step', 'prescribed_values.html#usage'),
                SettingsDictEntry("timeStepOutputInterval", '10', 'if the time step should be written to console, a value > 1 produces no output', 'prescribed_values.html#usage'),
                SettingsDictEntry("slotNames", '[]', 'list of names of the connector slots, maximum length is 6 characters per name', 'prescribed_values.html#slotnames'),
                SettingsDictEntry("fieldVariables1", '[]', 'list of field variables that will get values assigned in every timestep, by the provided callback function', 'prescribed_values.html#usage'),
                SettingsDictEntry("fieldVariables2", '[]', 'list of field variables that will get values assigned in every timestep, by the provided callback function', 'prescribed_values.html#usage'),
                SettingsDictEntry("additionalArgument", 'None', 'a custom argument to the fieldVariables callback functions, this will be passed on as the last argument', 'prescribed_values.html#additionalargument'),
                SettingsChoice([], [
                    outputwriter
                ]),
            ]))
        ])
    },



    "OutputWriter::OutputSurface": {
        "runnable": True,
        "timeSteppingScheme": True,
        "template_arguments": [
            ('Nested solver', ["timeSteppingScheme",
                      "SpatialDiscretization::FiniteElementMethod"])
        ],
        "python_options": SettingsDict([
            SettingsDictEntry("OutputSurface", SettingsDict([
                outputwriter,
                SettingsDictEntry("face", '["1-"]', 'the 2D faces of the 3D volume which will be extracted', 'output_surface.html#face'),
                SettingsDictEntry("samplingPoints", 'None', 'the electrode positions', 'output_surface.html#surface-sampling'),
                SettingsChoice([], [
                    SettingsDictEntry("updatePointPositions", 'False', 'the electrode points should be initialize in every timestep (set to False for the static case). This makes a difference if the muscle contracts, then True=fixed electrodes, False=electrodes moving with muscle.', 'output_surface.html#surface-sampling'),
                    SettingsDictEntry("filename", 'out/sampledPoints.csv', None, 'output_surface.html#surface-sampling'),
                    SettingsDictEntry("enableCsvFile", 'True', 'if the values at the sampling points should be written to csv files', 'output_surface.html#surface-sampling'),
                    SettingsDictEntry("enableVtpFile", 'True', 'if the values at the sampling points should be written to vtp files', 'output_surface.html#surface-sampling'),
                    SettingsDictEntry("xiTolerance", '0.3', 'tolerance for element-local coordinates xi, for finding electrode positions inside the elements. Increase or decrease this numbers if not all electrode points are found', 'output_surface.html#surface-sampling'),
                    SettingsDictEntry("enableGeometryInCsvFile", 'True', 'if the csv output file should contain geometry of the electrodes in every time step. This increases the file size and only makes sense if the geometry changed throughout time, i.e. when computing with contraction', 'output_surface.html#enablegeometryincsvfile'),
                ]),
                SettingsChildPlaceholder(0)
            ]))
        ])
    },


    "SpatialDiscretization::FiniteElementMethod": {
        "runnable": True,
        "discretizableInTime": True,
        "template_arguments": [
            ('Mesh', ["Mesh::"]),
            ('BasisFunction', ["BasisFunction::"]),
            ('Quadrature', ["Quadrature::"]),
            ('Equation', ["Equation::"])
        ],
        "python_options": SettingsDict([
            SettingsDictEntry("FiniteElementMethod", SettingsDict([
                SettingsChildPlaceholder(0),
                SettingsDictEntry("slotName", '""', 'specifies the name of the slot that contains the solution variable', 'finite_element_method.html#slotname'),
                SettingsDictEntry("prefactor", '1', 'a scalar multiplier of the Laplace operator term, i.e. c in c⋅Δu or c⋅∇⋅(A∇u)', 'finite_element_method.html#prefactor'),
                SettingsDictEntry("rightHandSide", '{}', 'right hand side vector f (either as list or as dict)', 'finite_element_method.html#righthandside'),
                SettingsDictEntry("dirichletBoundaryConditions", '{}', 'a dict with {<dof no>: <value>} entries', 'boundary_conditions.html#dirichlet-boundary-conditions'),
                SettingsDictEntry("dirichletOutputFilename", 'None', 'write Dirichlet Boundary conditions to .vtp file', 'boundary_conditions.html#dirichlet-output-filename'),
                SettingsDictEntry("neumannBoundaryConditions", '[]', None, 'boundary_conditions.html#neumann-boundary-conditions'),
                # this is undocumented:
                SettingsDictEntry("divideNeumannBoundaryConditionValuesByTotalArea", 'False', 'if the neumann boundary condition vectors should be divided by the total surface area where surface loads are applied, this allows to specify the total force that acts on the surface. If set to False (default), the given traction is a per-surface quantity'),
                # only used in quasi_static_linear_elasticity_solver
                SettingsChoice([], [
                    SettingsDictEntry("bulkModulus", '1.0', 'bulk modulus, K, material parameter for compressibility','quasi_static_linear_elasticity_solver.html#python-settings'),
                ]),
                SettingsChoice([], [
                    SettingsDictEntry("shearModulus", '1.0', 'shear modulus, μ','quasi_static_linear_elasticity_solver.html#python-settings'),
                ]),
                SettingsDictEntry(
                    "updatePrescribedValuesFromSolution", 'False', 'if set to true, the values that are initially set in the solution field variable are used as the prescribed values at the dofs in dirichletBoundaryConditions', 'finite_element_method.html#updateprescribedvaluesfromsolution'),
                SettingsDictEntry("inputMeshIsGlobal", 'True', 'together with rightHandSide it specifies whether the given values are interpreted as local values or global values in the context of a parallel execution on multiple processes', 'finite_element_method.html#inputmeshisglobal'),
                SettingsChoice([], [
                    SettingsDictEntry("diffusionTensor", '[]', 'for anisotropic diffusion, the diffusion tensor can be given as a list of double valus in row-major order', 'finite_element_method.html#diffusiontensor')
                ]),
                # undocumented in FiniteElementMethod
                SettingsChoice([], [
                    SettingsDictEntry("extracellularDiffusionTensor", '[]', 'sigma_e, one list item = same tensor for all elements, multiple list items = a different tensor for each element', 'static_bidomain_solver.html#python-settings')
                ]),
                solver_linear,
                SettingsChoice([], [
                    outputwriter
                ]),
            ]))
        ])
    },

    "Mesh::StructuredRegularFixedOfDimension": {
        "template_arguments": [
            ("Dimension", ["1", "2", "3"])
        ],
        "python_options": SettingsDict([
            SettingsMesh([
                SettingsDictEntry("nElements", SettingsList(['1', '1']), 'the number of elements of the mesh in the coordinate directions', 'mesh.html#nelements'),
                SettingsDictEntry("physicalExtent", '[1.0, 1.0]', 'the “size” of the mesh in physical units (e.g. meters if SI units are used), in the coordinate directions', 'mesh.html#physicalextent'),
                #SettingsDictEntry("physicalOffset", '[0, 0]', 'offsets all node positions by the given vector. The physicalOffset is always global, i.e. even when inputMeshIsGlobal is set to False, the physicalOffset specifies the offset to the front lower left corner of the global mesh.', 'mesh.html#physicalextent'),
                SettingsDictEntry("inputMeshIsGlobal", 'True', 'whether the values of nElements and physicalExtent describe the global domain (True) or the local subdomain (False) in parallel execution', 'mesh.html#inputmeshisglobal')
            ])
        ])
    },

    "Mesh::StructuredDeformableOfDimension": {
        "template_arguments": [
            ("Dimension", ["1", "2", "3"])
        ],
        "python_options": SettingsDict([
            SettingsMesh([
                SettingsDictEntry("nElements", SettingsList(['1', '1']), 'the number of elements of the mesh in the coordinate directions', 'mesh.html#id1'),
                SettingsDictEntry("inputMeshIsGlobal", 'True', 'whether the values of nElements and the nodePositions describe the global domain (True) or the local subdomain (False) in parallel execution', 'mesh.html#id2'),
                SettingsChoice([
                    SettingsDictEntry("physicalExtent", '[1.0, 1.0]', 'the “size” of the mesh in physical units (e.g. meters if SI units are used), in the coordinate directions', 'mesh.html#physicalextent'),
                    SettingsDictEntry("physicalOffset", '[0.0, 0.0]', 'offset all node positions by the given vector', 'mesh.html#structureddeformable')
                ], [
                    SettingsDictEntry("nodePositions", '[[0,0,0], [0,0,0]]', 'specify all node positions', 'mesh.html#nodepositions')
                ])
            ])
        ])
    },

    "Mesh::UnstructuredDeformableOfDimension": {
        "template_arguments": [
            ('Dimension', ["1", "2", "3"])
        ],
        "python_options": SettingsDict([
            SettingsMesh([
                SettingsChoice([
                    SettingsDictEntry("elements", '[[0,0]]', 'list of nodes, each node is [node no, version-at-that-node no] or just node-no then it assumes version no 0', 'mesh.html#elements'),
                    SettingsDictEntry("nodePositions", '[[0,0]]', 'list of positions of the nodes, each node position is a list with maximum three entries for the components in x,y and z direction', 'mesh.html#id5')
                ], [
                    SettingsDictEntry("exelem", "filename.exelem", 'the filename of the exelem file', 'mesh.html#exelem'),
                    SettingsDictEntry("exnode", "filename.exnode", 'the filename of the exnode file', 'mesh.html#exnode')
                ])
            ])
        ])
    },

    "Mesh::CompositeOfDimension": {
        "template_arguments": [
            ('Dimension', ["1", "2", "3"])
        ],
        "python_options" : SettingsDict([
            SettingsDictEntry("meshName", '["meshX", "meshY"]', 'a list containing all mesh names of the submeshes', 'mesh.html#compositeofdimension-d')
        ])
    },


    "BasisFunction::CompletePolynomialOfDimensionAndOrder": {
        "template_arguments": [
            ('Dimension', ["1", "2", "3"]),
            ('Order', ["0", "1", "2"])
        ]
    },
    "BasisFunction::Hermite": {},
    "BasisFunction::LagrangeOfOrder": {
        "template_arguments_needed": 0,
        "template_arguments": [
            ('Order', ["1", "2"])
        ]
    },
    # TODO are there BasisFunction::Mixed?


    "Quadrature::None": {},
    "Quadrature::ClenshawCurtis": {
        "template_arguments": [
            ('NumberIntegrationPoints', ["1", "2", "3", "4", "5", "6", "7", "64"])
        ]
    },
    "Quadrature::Gauss": {
        "template_arguments": [
            ('NumberGaussPoints', ["1", "2", "3", "4", "5", "6", "7",
                      "8", "10", "12", "16", "20", "24", "64"])
        ]
    },
    "Quadrature::NewtonCotes": {
        "template_arguments": [
            ('NumberIntegrationPoints', ["1", "2", "3", "4", "5", "6", "7", "8"])
        ]
    },
    "Quadrature::TensorProduct": {
        "template_arguments": [
            ('Dimension', ["1", "2", "3"]),
            ('Quadrature', ["Quadrature::"])
        ]
    },
    # TODO are there Quadrature::Mixed?
    # "Quadrature::Mixed" : {
    #    "template_arguments" : [
    #       [ ],
    #       [ "Quadrature::ClenshawCurtis", "Quadrature::Gauss", "Quadrature::NewtonCotes" ]
    #    ]
    # },


    "Equation::Dynamic::IsotropicDiffusion": {},
    "Equation::Dynamic::AnisotropicDiffusion": {},
    "Equation::Dynamic::DirectionalDiffusion": {},
    "Equation::Static::Laplace": {},
    "Equation::Static::GeneralizedLaplace": {},
    "Equation::Static::LinearElasticity": {},
    "Equation::Static::LinearElasticityActiveStress": {},
    # TODO are these correct?
    "Equation::SolidMechanics::MooneyRivlinIncompressible3D": {},
    "Equation::SolidMechanics::TransverselyIsotropicMooneyRivlinIncompressible3D": {},
    "Equation::SolidMechanics::TransverselyIsotropicMooneyRivlinIncompressibleActive3D": {},
    "Equation::SolidMechanics::HyperelasticTendon": {},
    "Equation::SolidMechanics::HyperelasticityBase": {},
    "Equation::Static::Poisson": {},
    "Equation::Static::GeneralizedPoisson": {},
}
