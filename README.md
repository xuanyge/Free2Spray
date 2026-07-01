This automated script is developed for large-scale cold spray deposition simulations based on the Eulerian approach in ABAQUS.

For the detailed implementation procedure, please refer to the attached paper. Before running the script, carefully verify the file paths and parameter settings specified in the script. The default configuration corresponds to the representative case study presented in Section 4 of the paper.

For other simulation scenarios, please modify the following settings accordingly:

Eulerian domain and substrate geometry
Particle and substrate velocity and temperature
Mesh size
Analysis step time
Boundary conditions
Number of execution loops

The current version of the script uses the GE model and Cu as the material by default. If other material models or materials are required, additional modifications to the script and/or calibration of material parameters will be necessary.

In addition, when the model contains more than 10 million mesh elements, a computer with at least 32 GB of RAM is required.

If you have any questions, please feel free to contact: xuanyu.ge@polimi.it.

If you use this program, please cite the following paper,

Free2Spray: X. Ge, L. Zhou, A.A. Lordejani, A.H. Astaraee, S. Bagherifard, M. Guagliano, Large-Scale Multi-Particle Cold Spray Simulation Framework for Deposit Morphology and Deformation Analysis, Additive Manufacturing  (2026) 105285.

GE material model: X. Ge, L. Zhou, S. Bagherifard, M. Guagliano, High-Fidelity modeling of cold spray: Improved constitutive material model and nonlocal mesh sensitivity mitigation, International Journal of Mechanical Sciences 321 (2026) 111654.
