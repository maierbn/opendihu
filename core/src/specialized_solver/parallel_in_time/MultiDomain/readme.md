In this section all necessary functions are implemented for MGRIT applied to the Multidomain solver.
MGRIT is implemented with the help of XBRAID.
multidomain_wrapper: implements the interface for Xbraid having access to the multidomain solver, which is an Opendihu solver.
PinT_fun: This file includes basic functions for XBraid and MGRIT. For example: How to create a vector or access a vector. The header file defines the main structures of XBraid: The Braid app and the Braid vector.
PinT_lib: This file includes helper methods. For example: printing output.
PinT_MD: This file is the main file. It creates the problem with the help of Opendihu and starts the solver.
PinT_MD_Braid: This file includes the two most important, user-defined XBraid functions: Braid Step and Braid Init.