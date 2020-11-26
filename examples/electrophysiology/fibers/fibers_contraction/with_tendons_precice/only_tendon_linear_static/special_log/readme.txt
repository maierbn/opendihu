/store/software/opendihu/examples/electrophysiology/fibers/fibers_contraction/with_tendons_precice/only_tendon_linear_static


./tendon_linear_precice_quasistatic settings_tendon_linear.py

This is opendihu 1.1, built Nov 10 2020, C++ 201402, GCC 7.5.0, current time: 2020/11/10 16:21:31, hostname: lapsgs05, n ranks: 1
Open MPI v2.1.1, package: Open MPI buildd@lcy01-amd64-009 Distribution, ident: 2.1.1, repo rev: v2.1.0-100-ga2fdb5b, May 10, 2017
File "settings_tendon_linear.py" loaded.
---------------------------------------- begin python output ----------------------------------------
Load: 100.0 N/cm^2 = 1.0 MPa
n fibers:              25 (5 x 5)
n points per fiber:    21
1 rank, partitioning: x1 x y1 x z1
5 x 5 = 25 fibers, per partition: 4 x 4 = 16
per fiber: 1D mesh    nodes global: 21, local: 21
  sampling 3D mesh with stride 1 x 1 x 1 
    linear 3D mesh    nodes global: 5 x 5 x 21 = 525, local: 5 x 5 x 21 = 525
    linear 3D mesh elements global: 4 x 4 x 20 = 320, local: 4 x 4 x 20 = 320
 quadratic 3D mesh    nodes global: 5 x 5 x 21 = 525, local: 5 x 5 x 21 = 525
 quadratic 3D mesh elements global: 2 x 2 x 10 = 40, local: 2 x 2 x 10 = 40
number of degrees of freedom:
                    1D fiber:         21  (per process: 21)
            0D-1D monodomain:         84  (per process: 84)
 all fibers 0D-1D monodomain:       2100  (per process: 1344)
                 3D bidomain:        525  (per process: 525)
                       total:       2625  (per process: 1869)
nRanks:  [1, 1, 1]
----------------------------------------- end python output -----------------------------------------
Read from file "../../../../input/left_biceps_brachii_tendon1.bin", 1075 collective chunks.
done.
You have specified the solver in-line and not under the extra key "Solvers". You could do so,  by defining "Solvers": {"<your custom solver name>": {<your solver parameters>}} at the beginning of the  config and "solverName": "<your custom solver name>" where you currently have specified the solver parameters.  This is required if you want to use the same solver for multiple objects.
Preallocation for matrix "combinedJacobian": diagonal nz: 262440, offdiagonal nz: 262440
Warning: Coupling is disabled (option "couplingEnabled": False).

  Nonlinear solver: iteration  0, residual norm     94.6979
  Nonlinear solver: iteration  1, residual norm      469488, e_new=e_old^c with c=2.86976
  Nonlinear solver: iteration  2, residual norm      126840, e_new=e_old^c with c=0.899787
  Nonlinear solver: iteration  3, residual norm     35448.2, e_new=e_old^c with c=0.891508
  Nonlinear solver: iteration  4, residual norm     8946.23, e_new=e_old^c with c=0.86857
  Nonlinear solver: iteration  5, residual norm     6578.12, e_new=e_old^c with c=0.966207
  Nonlinear solver: iteration  6, residual norm      3575.2, e_new=e_old^c with c=0.930646
  Nonlinear solver: iteration  7, residual norm     2916.38, e_new=e_old^c with c=0.975106
  Nonlinear solver: iteration  8, residual norm     6940.43, e_new=e_old^c with c=1.10868
  Nonlinear solver: iteration  9, residual norm     1457.78, e_new=e_old^c with c=0.823581
  Nonlinear solver: iteration 10, residual norm     1436.99, e_new=e_old^c with c=0.998027
  Nonlinear solver: iteration 11, residual norm     1323.32, e_new=e_old^c with c=0.988665
  Nonlinear solver: iteration 12, residual norm     1131.94, e_new=e_old^c with c=0.978268
  Nonlinear solver: iteration 13, residual norm     1038.16, e_new=e_old^c with c=0.987701
  Nonlinear solver: iteration 14, residual norm     738.121, e_new=e_old^c with c=0.950887
  Nonlinear solver: iteration 15, residual norm     572.288, e_new=e_old^c with c=0.961469
  Nonlinear solver: iteration 16, residual norm     470.648, e_new=e_old^c with c=0.969206
  Nonlinear solver: iteration 17, residual norm     400.247, e_new=e_old^c with c=0.973672
  Nonlinear solver: iteration 18, residual norm     347.222, e_new=e_old^c with c=0.976282
  Nonlinear solver: iteration 19, residual norm     304.944, e_new=e_old^c with c=0.977806
  Nonlinear solver: iteration 20, residual norm     269.925, e_new=e_old^c with c=0.978674
  Nonlinear solver: iteration 21, residual norm     240.247, e_new=e_old^c with c=0.979194
  Nonlinear solver: iteration 22, residual norm     214.865, e_new=e_old^c with c=0.979631
  Nonlinear solver: iteration 23, residual norm       193.2, e_new=e_old^c with c=0.980207
  Nonlinear solver: iteration 24, residual norm     174.828, e_new=e_old^c with c=0.981017
  Nonlinear solver: iteration 25, residual norm     159.307, e_new=e_old^c with c=0.981996
  Nonlinear solver: iteration 26, residual norm     146.157, e_new=e_old^c with c=0.98301
  Nonlinear solver: iteration 27, residual norm     134.924, e_new=e_old^c with c=0.983957
  Nonlinear solver: iteration 28, residual norm     125.223, e_new=e_old^c with c=0.984787
  Nonlinear solver: iteration 29, residual norm     116.749, e_new=e_old^c with c=0.985493
  Nonlinear solver: iteration 30, residual norm     109.263, e_new=e_old^c with c=0.986079
  Nonlinear solver: iteration 31, residual norm      102.58, e_new=e_old^c with c=0.986554
  Nonlinear solver: iteration 32, residual norm     96.5547, e_new=e_old^c with c=0.986927
  Nonlinear solver: iteration 33, residual norm     91.0711, e_new=e_old^c with c=0.987206
  Nonlinear solver: iteration 34, residual norm      86.037, e_new=e_old^c with c=0.987396
  Nonlinear solver: iteration 35, residual norm     81.3775, e_new=e_old^c with c=0.987501
  Nonlinear solver: iteration 36, residual norm     77.0312, e_new=e_old^c with c=0.987523
  Nonlinear solver: iteration 37, residual norm     72.9474, e_new=e_old^c with c=0.987461
  Nonlinear solver: iteration 38, residual norm     69.0836, e_new=e_old^c with c=0.987314
  Nonlinear solver: iteration 39, residual norm      65.404, e_new=e_old^c with c=0.987077
  Nonlinear solver: iteration 40, residual norm     61.8772, e_new=e_old^c with c=0.986741
  Nonlinear solver: iteration 41, residual norm     58.4757, e_new=e_old^c with c=0.986294
  Nonlinear solver: iteration 42, residual norm     55.1738, e_new=e_old^c with c=0.985714
  Nonlinear solver: iteration 43, residual norm     51.9465, e_new=e_old^c with c=0.984971
  Nonlinear solver: iteration 44, residual norm     48.7676, e_new=e_old^c with c=0.984014
  Nonlinear solver: iteration 45, residual norm     45.6074, e_new=e_old^c with c=0.982764
  Nonlinear solver: iteration 46, residual norm     42.4283, e_new=e_old^c with c=0.981086
  Nonlinear solver: iteration 47, residual norm     39.1785, e_new=e_old^c with c=0.978737
  Nonlinear solver: iteration 48, residual norm     35.7781, e_new=e_old^c with c=0.975248
  Nonlinear solver: iteration 49, residual norm     32.0901, e_new=e_old^c with c=0.96959
  Nonlinear solver: iteration 50, residual norm     27.8472, e_new=e_old^c with c=0.959114
  Nonlinear solver: iteration 51, residual norm     22.4064, e_new=e_old^c with c=0.934655
  Nonlinear solver: iteration 52, residual norm     13.7817, e_new=e_old^c with c=0.843695
  Nonlinear solver: iteration 53, residual norm     3.88614, e_new=e_old^c with c=0.517439
  Nonlinear solver: iteration 54, residual norm    0.949822, e_new=e_old^c with c=-0.0379251
  Nonlinear solver: iteration 55, residual norm     0.18313, e_new=e_old^c with c=32.975
  Nonlinear solver: iteration 56, residual norm   0.0547614, e_new=e_old^c with c=1.71115
  Nonlinear solver: iteration 57, residual norm   0.0131692, e_new=e_old^c with c=1.49061
  Nonlinear solver: iteration 58, residual norm  0.00452367, e_new=e_old^c with c=1.24679
  Nonlinear solver: iteration 59, residual norm  0.00121809, e_new=e_old^c with c=1.24304
  Nonlinear solver: iteration 60, residual norm 0.000444613, e_new=e_old^c with c=1.15019
Solution done in 60 iterations, residual norm 0.000444613: SNES_CONVERGED_SNORM_RELATIVE: Newton computed step size small; || delta x || < stol || x ||, KSP_CONVERGED_ITERATING: returned if the solver is not yet finished
File "out/solver_structure.txt" written.
Total user time: 25s
File "out/mappings_between_meshes_log.txt" written.
