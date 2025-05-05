# Protecting Information in a Parametrically Driven Hybrid Quantum System

This repository contains the numerical codes and data used in the paper: [arXiv:2207.14354](https://arxiv.org/abs/2207.14354)
(Authors: Siddharth Tiwary, Harsh Sharma, Himadri Shekhar Dhar)

---

## Repository Structure

### 1. `vDMRG/`  
This folder contains C++ and Python source codes and simulation data used for the **variational Density Matrix Renormalization Group (vDMRG)** calculations. These simulations correspond to **Figure 3** and **Figure 4** in the paper.

#### Contents:
- `src/`: C++ source code files for the vDMRG simulation  
- `data.zip`: Output data files generated from simulations  
- `post-processing/`: Python scripts for post-processing and plotting  
- `README.txt`: Specific instructions and notes about vDMRG codes and parameters

#### Dependencies:
- **C++ Code**: Runs vDMRG simulation
  - Requires the **[Armadillo](https://arma.sourceforge.net/)** C++ library for linear algebra & scientific computing.
  - A C++11 compatible compiler (e.g., `g++`)

> Armadillo Citation:  
> C. Sanderson and R. Curtin, *Math. Comput. Appl.* **24(3)**, 70 (2019).

- **Python Code**: Required for data post-processing
  - Uses **[QuTiP](http://qutip.org/)** (Quantum Toolbox in Python)

> QuTiP Citation:  
> J. R. Johansson, P. D. Nation, and F. Nori, *Comput. Phys. Commun.* **183**, 1760 (2012);  
> *Comput. Phys. Commun.* **184**, 1234 (2013).

#### How to Compile and Run the C++ Code:
```bash
# Navigate to the vDMRG folder
cd vDMRG/src

# Compile the code (example: main.cpp as the entry file)
g++ --std=c++11 Squeezed_DMRG.cpp FF.cpp -o vdmrg_exec -O3 -larmadillo

# Run the executable
./vdmrg_exec < ParameterFile.txt
```

> Note: Ensure that the Armadillo library is properly installed and accessible on your system. You may need to add `-I` and `-L` flags to specify include/library paths if Armadillo is installed in a custom location. Data files generated from C++ code is saved in the current directory, ParameterFile.txt and the executable can be placed and run in proper folders accordingly.

---

### 2. `Semiclassical/`  
This folder contains **Mathematica notebooks and data** related to the **semiclassical mean-field analysis** of the Lindblad master equation dynamics. These simulations correspond to **Figure 2** in the paper.

#### Contents:
- `Semiclassical_Dynamics.nb`: Main notebook to reproduce semiclassical data.  
- `data/`: Contains pre-generated data for analysis

---

## Notes

- If you use this code or data in your work, please cite the original paper.
- For issues or contributions, feel free to open a GitHub issue or pull request.
