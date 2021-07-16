[![CircleCI](https://circleci.com/gh/maierbn/opendihu/tree/develop.svg?style=svg)](https://circleci.com/gh/maierbn/opendihu/tree/develop)
[![Build Status](https://travis-ci.com/maierbn/opendihu.svg?branch=stable)](https://travis-ci.com/maierbn/opendihu)
[![CodeFactor](https://www.codefactor.io/repository/github/maierbn/opendihu/badge/develop)](https://www.codefactor.io/repository/github/maierbn/opendihu/overview/develop)
[![ReadTheDocs](https://readthedocs.org/projects/opendihu/badge/?version=latest)](https://opendihu.readthedocs.io/en/latest/)

Link to Documentation: https://opendihu.readthedocs.io/en/latest/

# Overview
OpenDiHu is a software framework to solve 1D, 2D, and 3D multi-physics problems in parallel with the Finite Element Method.
It is used in the domain of skeletal muscle simulations: Electrophysiology, contraction, neuro-chemo-electro-mechanics.
Design goals are usability, performance and extensibility.

The software is developed at [SGS](https://www.ipvs.uni-stuttgart.de/abteilungen/sgs/index.html?__locale=en) and [IANS](https://www.ians.uni-stuttgart.de/institute/) at the [University of Stuttgart](https://www.uni-stuttgart.de/en/index.html).

# Installation
Refer to the documentation for detailed [installation instructions](https://opendihu.readthedocs.io/en/latest/user/installation.html).

However, if you usually skip instructions, do the following:
```
git clone https://github.com/maierbn/opendihu.git && cd opendihu && make
```
If there are error messages, have a look at the log file `config.log`.

# Literature / How to cite

View the dissertation of Benjamin Maier, submitted in June 2021:
[Maier, B. (2021). Scalable Biophysical Simulations of the Neuromuscular System](https://arxiv.org/abs/2107.07104)

Currently, we have this conference paper from 2019:
[Maier, B., Emamy, N., Krämer, A., & Mehl, M. (2019). Highly parallel multi-physics simulation of muscular activation and EMG. In COUPLED VIII: proceedings of the VIII International Conference on Computational Methods for Coupled Problems in Science and Engineering (pp. 610-621). CIMNE](https://upcommons.upc.edu/handle/2117/190149)

More descriptive literature is expected to be published in late 2021.
