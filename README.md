# Beyond-Phasors-Solving-Non-Sinusoidal-Electrical-Circuits-using-Geometry
Direct solution of RLC circuits with GA
**"Beyond Phasors: Solving Non-Sinusoidal Electrical Circuits using Geometry"**
*Javier Castillo-Martínez, Raul Baños and Francisco G. Montoya*

*(Link to article/DOI - [ADD LINK/DOI WHEN PUBLISHED])*

## 📂 Repository Contents

1.  **`Material_Suplementary___Beyond_Phasors.pdf`**:
    [cite_start]The complete supplementary material for the article[cite: 381]. It includes:
    * [cite_start]Detailed derivations for the single-frequency (N=1) validation[cite: 390].
    * [cite_start]Derivations for canonical (R, L, C) and combined (RL, RC) loads[cite: 443, 506].
    * [cite_start]Full derivations of the N=2 and N=3 case studies presented in the main paper[cite: 552, 608].
    * [cite_start]Extended case studies (N=5) for series and parallel circuits[cite: 650, 692].
    * [cite_start]**Appendix H**: A detailed description of the script implementations, required libraries, and the computational performance benchmark results[cite: 723].

2.  **`/scripts`** (Suggested directory for the .m files):
    * **`comparativa_tiempos_egipcio_rev08.m`**:
        [cite_start]The *benchmarking* script (v24.4) used to generate the computational performance comparisons shown in Appendix H[cite: 765, 787]. This script measures the execution time of the "GA-Real (Native)" algorithm against the classic "Phasor (Native Complex)" method for all four RLC circuit cases.
    * **`Resolution_GA_clifford_rev15.m`**:
        A *validation* script (v17.1) that verifies the numerical consistency between three methods: (1) GA-Real (Native), (2) Phasors (Native), and (3) the high-level Roto-flex method implemented using the **`Clifford`** library for MATLAB.
    * **`Resolution_GA_sugar_rev15.m`**:
        A *validation* script (v19.1) that performs the same numerical check as the one above but uses the **`SUGAR`** library for MATLAB.

## ⚙️ Dependencies

To run the scripts in this repository, you will need:

* **MATLAB**
* [cite_start]The following Geometric Algebra libraries for MATLAB (as mentioned in the scripts and Appendix H [cite: 724, 731, 748]):
    * **`Clifford`** (Used by `Resolution_GA_clifford_rev15.m`)
    * **`SUGAR`** (Used by `Resolution_GA_sugar_rev15.m`)
    * **`GA-FuL`** (Used by `comparativa_tiempos_egipcio_rev08.m` for vector parsing)

## 🚀 How to Use

1.  Clone or download this repository.
2.  Ensure you have MATLAB installed.
3.  Download the `Clifford`, `SUGAR`, and `GA-FuL` libraries and add them to your MATLAB path.
4.  To **validate the numerical results**, run `Resolution_GA_clifford_rev15.m` or `Resolution_GA_sugar_rev15.m` and follow the on-screen prompts to select the desired case study.
5.  To **reproduce the performance benchmarks**, run `comparativa_tiempos_egipcio_rev08.m`.
