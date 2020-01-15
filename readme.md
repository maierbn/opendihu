[![CircleCI](https://circleci.com/gh/maierbn/opendihu/tree/develop.svg?style=svg)](https://circleci.com/gh/maierbn/opendihu/tree/develop)
[![Build Status](https://travis-ci.com/maierbn/opendihu.svg?branch=stable)](https://travis-ci.com/maierbn/opendihu)
[![CodeFactor](https://www.codefactor.io/repository/github/maierbn/opendihu/badge/develop)](https://www.codefactor.io/repository/github/maierbn/opendihu/overview/develop)
[![ReadTheDocs](https://readthedocs.org/projects/opendihu/badge/?version=latest)](https://opendihu.readthedocs.io/en/latest/)

Link to Documentation: https://opendihu.readthedocs.io/en/latest/

# Overview
Opendihu is a software framework that solves static and dynamic multi-physics problems, spatially discretized in 1D, 2D and 3D by the finite element method.
The core design goals are usability, performance and extensibility. 

It is developed at [SGS](https://www.ipvs.uni-stuttgart.de/abteilungen/sgs/index.html?__locale=en) and [IANS](https://www.ians.uni-stuttgart.de/institute/) at [University of Stuttgart](https://www.uni-stuttgart.de/en/index.html)
and funded by the [DFG](https://www.dfg.de/en/funded_projects/current_projects_programmes/list/projectdetails/index.jsp?id=277536708) as part of the IRTG ["Soft Tissue Robotics"](https://www.irtg.auckland.ac.nz/) and the [BW-Stiftung](https://www.bwstiftung.de/) under the [HPC 2](https://www.bwstiftung.de/hpcii/#c15824) call.

Opendihu is the compatible HPC version of [OpenCMISS](http://opencmiss.org/). The native format of OpenCMISS, the `Exfile` format, can be used for input and output in Opendihu.
The implemented skeletal muscle models are the same and are based on the same [CellML](https://www.cellml.org/) descriptions of subcellular models. 
For details, read the [introduction on the documentation website](https://opendihu.readthedocs.io/en/latest/introduction.html).

# Installation
Documentation including detailed [installation instructions](https://opendihu.readthedocs.io/en/latest/user/installation.html) can be found at [opendihu.readthedocs.io](https://opendihu.readthedocs.io/en/latest/).

If you usually skip instructions (and have Ubuntu), try the following
```
git clone https://github.com/maierbn/opendihu.git && cd opendihu
sudo apt-get update && \
  sudo apt-get install -y gfortran libopenmpi-dev libx11-* python2.7 git apt-utils make software-properties-common zlib1g-dev cmake libssl-dev bison flex
make
```
and see what happens. If there are error messages, look into the log file `config.log`.
