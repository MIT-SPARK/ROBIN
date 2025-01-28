# ROBIN 

ROBIN is a library for outlier rejection based on compatibility graphs.

If you find this library helpful or use it in your projects, please cite:
```bibtex
@InProceedings{Shi21icra-robin,
	title={{ROBIN:} a Graph-Theoretic Approach to Reject Outliers in Robust Estimation using Invariants},
	author={J. Shi and H. Yang and L. Carlone},
	booktitle={IEEE Intl. Conf. on Robotics and Automation (ICRA)},
	note = {arXiv preprint: 2011.03659},
	pdf={https://arxiv.org/pdf/2011.03659.pdf},
	year={2021}
}
```
and
```bibtex
@article{Shi22arxiv-PACE,
  author = {J. Shi and H. Yang and L. Carlone},
  title = {Optimal and Robust Category-level Perception: Object Pose and Shape Estimation from {2D and 3D} Semantic Keypoints},
  journal = {arXiv preprint: 2206.12498},
  pdf = {https://arxiv.org/pdf/2206.12498.pdf},
  Year = {2022}
}
```

If you are interested in more works from us, please visit our lab page [here](http://web.mit.edu/sparklab/).

--- 

# :gear: Build & Installation

## :package: Dependency 

ROBIN has the following dependencies:
1. OpenMP
2. Eigen3

Thus, follow the below commandline:

```bash
sudo apt-get install gcc g++ build-essential libeigen3-dev cmake python3-pip python3-dev git ninja-build -y
```

## C++ Installation 

Run the following to build the library using CMake (inside the repo root directory):
```bash
mkdir build && cd build
cmake .. && make
sudo make install
```
The following CMake options are provided:
```
BUILD_DOCS: Build documentation. Default: OFF
BUILD_TESTS: Enable testing with ctest. Default: ON
BUILD_MATLAB_BINDINGS: Build MATLAB bindings. Default: OFF
USE_ASAN: Enable address sanitizer. Default: OFF
ENABLE_DIAGNOSTIC_PRINT: Enable printing of diagnostic messages. Default: OFF
```

## Python Installation 

It's simple! To install Python bindings, just run:

```bash
pip3 install -e python/
```

--- 

# Third-party Data
Some of the testing data are from the [Network Repository](http://networkrepository.com/index.php). 
For more information, please refer to: 

Rossi, Ryan, and Nesreen Ahmed. "The network data repository with interactive graph analytics and visualization." Twenty-Ninth AAAI Conference on Artificial Intelligence. 2015.

# Third-party Code
- Parallel Maximum Clique (PMC) Library: 
    - License: GNU General Public License (https://github.com/ryanrossi/pmc/blob/master/LICENSE.md) 
- pybind11: https://github.com/pybind/pybind11
    - License: BSD-style (https://github.com/pybind/pybind11/blob/master/LICENSE)
- Catch2: https://github.com/catchorg/Catch2
    - License: Boost (https://github.com/catchorg/Catch2/blob/devel/LICENSE.txt)

# Known Issues
For errors in MATLAB with missing CXXABI:
```
export LD_PRELOAD=/usr/lib/gcc/x86_64-linux-gnu/7/libstdc++.so
```

# FAQ
* How to fix errors like "ModuleNotFoundError: No module named 'robin_py.robin_py'"?
  This might be caused by a mismatch between the Python interpreter versions the binding is built for and the interpreter 
  that the binding is installed on. Make sure to activate the correct virtual environment when calling `cmake ..`
