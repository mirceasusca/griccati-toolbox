# griccati-toolbox

The current project focuses on the implementation and validation of a
MATLAB toolbox to determine the existance and then, obtaining the
**stabilizable solution (X)**  and corresponding
**stabilizing positive feedback (F)** for matrix algebraic Riccati equations
(AREs) in both the continuous-time case and discrete-time case alike,
using the deflating subspaces of matrix pencils.

* For a continuous-time system and performance cost function
    *     dx(t)/dt = Ax(t) + Bu(t), with J = Integral{x'Qx + u'Ru + 2*x'Lu}dt, such that Re \Lambda{A+B*F} < 0;
* For a discrete-time system and performance cost function
    *     x[k+1] = Ax[k] + B[k], with J = Sum{x'Qx + u'Ru + 2*x'Lu}, such that |\Lambda{A+B*F}| < 1.


Besides the general-purpose ARE solver, a set of functionalities are also provided in the fields of
System Theory, along with Optimal and Robust Control, which are implemented in their
corresponding subfolders, described below.

Project structure:
-
* **control**: contains applications of the ARE solver for the continuous-time and
iscrete-time cases in Optimal and Robust Control, such as:
    * LQR, LQE, LQG with the fully-coupled structures and noise covariance specification;
    * H2, Hinf norm minimization problems;
    * it also has the LQG problem seen as a H2 minimization problem, instead of encompassing
    the individually obtained linear quadratic regulator and estimator.

* **dstools**: set of functionalities for descriptor systems, taken from [2];
of use in this project are *sl_gzero*, *sl_klf* (Fortran), *gklf*, *gsfstab* (MATLAB), for determining the numerical
rank of a matrix pencil, obtaining the Kronecker-like form of a possibly singular matrix pencil
(using an optimized staircase algorithm), and for generalized eigenvalue assignment.

* **pencils**: main functionalities of the toolbox, providing functions for matrix pencils and the
ARE solver.
    * computing the maximal proper stable deflating subspace of a matrix pencil in *mpdefsub.m*;
    * solves the generalized algebraic Riccati system (GCTARS/GDTARS) in *gricsv.m*;
    * creates a continuous or discrete-time Popov quadruplet, which encompasses an ARE problem;
    in [1] the authors use the term Popov triplet to encompass the A, B and P = [Q,L; L^T, R] matrices
    of an ARE, but the thesis uses a Popov "quadruplet", which also stores if the system is
    continuous or discrete;
    * computes the corresponding Hamiltonian and symplectic pencils for the
    linear quadratic optimization problem; their stable subspaces are used for solving the AREs
    in *create_hamiltonian_pencil.m*;
    * computes iterative refinement steps for the ARE solutions in *gricir.m*.

* **system_theory**: provides functions for checking system theoretical properties, such as:
    * verifying the bounded real lemma for a given system and gamma value;
    * checking if an ARE system has a positive (semi)definite stabilizing solution X;
    * computing normed coprime factorizations for a system;
    * scaling a Hinf problem to solve for gamma < 1 and backward rescaling;
    * formulating and solving a Kalman-Yakubovich-Popov (KPY) system in continuous-time,
    along with a Kalman-Szego-Popov-Yakubovich (KSPY) system in discrete-time;
    * verifying the small gain theorem for a positive feedback connection of two compatible systems;
    * inverting a system in state-space form.

* **tests**: various tests for AREs using benchmark problems, randomly-generated AREs with specified
size and comparison with MATLAB functions such as *lqr*, *care*, *dare*. Additionally, benchmark control
applications are tested for the control functions provided, such as the two tank control problem, active suspension
system, flexible beam, boost converter, inverted pendulum and others.

* **utils**: several utility functions such as checking a system for stabilizability,
detectability, positive definiteness (with specified tolerance), generalized plant
partitioning, extracting the Kronecker structure of an arbitrary matrix pencil etc.

Installation:
-
In order to run the toolbox functions, a Fortran compiler must be installed
and configured in MATLAB. See [2] and MATLAB *mex* function documentation. The implemented functions use the
routines *sl_gzero* and *sl_klf* from **dstools**. If these requirements are assured, the MATLAB functions
provided should be directly executable.

Observations:
-
The ARE solver and Optimal and Robust Control functions based on it have the possibility of using balancing,
which improve the obtained solutions.

All main functions described in the thesis are documented, are commented in relevant points and have commented
sequences of code which can be used for debugging.

The archive **pencils_old_geigassign.zip** containts an own implementation of the generalized
eigenvalue assignment algorithm using generalized Sylvester equations. In this state of the implementation, it does not
search for the numerical transformations with the minimum condition number.
The actual method used for the eigenvalue assignment problem for a controllable
matrix pencil is the routine *gsfstab* from **dstools**. Kept for educational purposes.

Bibliography:
-
1. V. Ionescu, C. Oară, M. Weiss, *Generalized Riccati Theory and Robust Control: A
Popov Function Approach*, John Wiley & Sons Ltd, England, 1999.

2. A. Varga, *DSTOOLS – The Descriptor System Tools for MATLAB*, 2018.
https://sites.google.com/site/andreasvargacontact/home/software/dstools.

An extended bibliography of the literature, encompassing also relevant implementation aspects,
is presented in the thesis, found at https://www.researchgate.net/profile/Mircea_Susca.
