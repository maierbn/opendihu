# Overview
The working title of this software framework is "opendihu" - from the project name "Digital Human". It serves as a code base to solve static and dynamic problems, where the Finite Element Method is used for spatial discretization. Due to its modular nature it is supposed to be adaptible for future problems.

# Relationship to OpenCMISS
* Probably `*.ex` files will be supported for the beginning to be able to use cmgui. (Currently a reader/writer for exelem/exnode files for unstructured meshes is implemented.)

# Possible items for the feature list
* T. Heidlauf multi-scale electro-chemo-mechanical model of skeletal muscle
* Addition of Bidomain EMG model of Milena 
* Multi-compartment model by Thomas
* Addition of neural feedback model
* Dynamic mechanical formulation

# Roadmap
* continuum mechanics model
* multi-compartment model
* Heidlauf model (wtihout EMG)
* shared memory parallelisation
* distributed memory parallelisation
* integration of VIS file format

# Documentation

## Interface classes
Interface classes promise a specific behaviour of the classes that derive from them.

### Runnable
Class represents a standalone simulation problem that can be run. These have to be the top-level classes in main() because it can be executed.
**Important methods:** `run()`
**Derived classes:**
`TimeSteppingScheme::ExplicitEuler<DiscretizableInTime>`
`OperatorSplitting::Godunov`
`OperatorSplitting::Strang`
`SpatialDiscretization::FiniteElementMethod` (computes a static problem like `Δu = f`)

### DiscretizableInTime
A time dependent problem of the form u_t = f(t). A time stepping scheme needs a class of this type. 
**Important methods:** `evaluateTimesteppingRightHandSide()` computes f(t)
**Derived classes:**
`CellmlAdapter`
`ModelOrderReduction::POD<DiscretizableInTimeType, PartType>`
`SpatialDiscretization::FiniteElementMethod<MeshType, BasisFunctionType, QuadratureType, Term>` for `Equation::hasLaplaceOperatorWithTimeStepping<Term>` (computes a dynamic problem like u_t = Δu)

### TimeSteppingScheme
A time stepping scheme can compute time steps with a time step width. This class is not a pure interface without functionality but already implements parsing of time step widths and number from the config.
**Important methods:** `advanceTimeSpan`
**Derived classes:**
`TimeSteppingScheme::ExplicitEuler<DiscretizableInTime>`
`OperatorSplitting::Godunov`

## Classes with state
Most classes only store helper variables and not real data. Classes that store objects are:
`TimeSteppingScheme::TimeSteppingSchemeOde<DiscretizableInTime>` stores object of class `DiscretizableInTime` (and therefore all derived timestepping schemes e.g. `ExplicitEuler`).
`ModelOrderReduction::PODBase<DiscretizableInTime>` stores object of class `DiscretizableInTime`.
`OperatorSplitting::Godunov<TimeStepping1, TimeStepping2>` stores objects of classes `TimeStepping1`, `TimeStepping2`.

## Config
