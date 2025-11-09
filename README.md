# Beyond Phasors: Solving Non-Sinusoidal Electrical Circuits using Geometry

This repository contains the supplementary material, validation scripts, and MATLAB benchmarks for the paper:

**"Beyond Phasors: Solving Non-Sinusoidal Electrical Circuits using Geometry"**
*Javier Castillo-Martínez, Raul Baños and Francisco G. Montoya*

*(Link to article/DOI - [ADD LINK/DOI WHEN PUBLISHED])*

## 📂 Repository Contents

1.  **`Material_Suplementary___Beyond_Phasors.pdf`**:
    The complete supplementary material for the article. It includes:
    * Detailed derivations for the single-frequency (N=1) validation.
    * Derivations for canonical (R, L, C) and combined (RL, RC) loads.
    * Full derivations of the N=2 and N=3 case studies presented in the main paper.
    * Extended case studies (N=5) for series and parallel circuits.
    * **Appendix H**: A detailed description of the script implementations, required libraries, and the computational performance benchmark results.

2.  **`/scripts`** (Suggested directory for the .m files):
    * **`Comparative_times_GAful.m`**:
        The *benchmarking* script used to generate the computational performance comparisons shown in Appendix H. This script measures the execution time of the "GA-Real (Native)" algorithm against the classic "Phasor (Native Complex)" method for all four RLC circuit cases.
    * **`Resolution_GA_clifford.m`**:
        A *validation* script that verifies the numerical consistency between three methods: (1) GA-Real (Native), (2) Phasors (Native), and (3) the high-level Roto-flex method implemented using the **`Clifford`** library for MATLAB.
    * **`Resolution_GA_sugar.m`**:
        A *validation* script that performs the same numerical check as the one above but uses the **`SUGAR`** library for MATLAB.

## ⚙️ Dependencies

To run the scripts in this repository, you will need:

* **MATLAB**
* The following Geometric Algebra libraries for MATLAB:
    * **`Clifford`** (Required by `Resolution_GA_clifford.m`)
    * **`SUGAR`** (Required by `Resolution_GA_sugar.m`)
    * **`GA-FuL`** (Required by `Comparative_times_GAful.m` for vector parsing)

## 🚀 How to Use

1.  Clone or download this repository.
2.  Ensure you have MATLAB installed.
3.  Download the `Clifford`, `SUGAR`, and `GA-FuL` libraries and add them to your MATLAB path.
4.  To **validate the numerical results**, run `Resolution_GA_clifford.m` or `Resolution_GA_sugar.m` and follow the on-screen prompts to select the desired case study.
5.  To **reproduce the performance benchmarks**, run `Comparative_times_GAful.m`.

1.  Clone or download this repository.
2.  Ensure you have MATLAB installed.
3.  Download the `Clifford`, `SUGAR`, and `GA-FuL` libraries and add them to your MATLAB path.
4.  To **validate the numerical results**, run `Resolution_GA_clifford_rev15.m` or `Resolution_GA_sugar_rev15.m` and follow the on-screen prompts to select the desired case study.
5.  To **reproduce the performance benchmarks**, run `comparativa_tiempos_egipcio_rev08.m`.
