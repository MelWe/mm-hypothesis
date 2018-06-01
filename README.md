# mm-hypothesis
Code repository for "Heuristic Framework for Testing the Multi-Manifold Hypothesis".

-----------

This directory is intended to be a reproducibility repository for Use Cases Described in the WiSDM paper: Heuristic Framework for Multi-Scale Testing of the Multi-Manifold Hypothesis. There are 3 subdirectories: datasets, code, and results. Datasets contains the the xyz coordinate file for the LIDAR image of the Golden Gate bridge, the MATLAB program that generates samples of the Sphere-Line configuration, and a spreadsheet of Word2Vec data (common_preprocessed.xlsx).
Code contains three sub-directories: vidimcode, vidimunpackedcode and MMMethodcode.
Vidimcode contains code for variance-based intrinsic dimension. There are two MATLAB programs: testvariancebaseddimLIDAR published on the Google Drive on August22 and computevidim. The latter program reorganizes the main calling program so it can be called more clearly by the MMMethodcode
MMMethodcode contains a main calling program: testMMmethodNov1.m which specifies and executes the uses cases in the paper (Sphere-Line and LIDAR). The other programs are invoked to construct the dyadic linear multi-manifold model for the data and the hypothesis testing distribution.
Vidimimpacked code contains the individual procedures used to compute intrinsic dimension. It also contains the idimtests program which is used to compute idim related statistics in a structure idimstats.
The results directory contains the files output by the execution of the method for the two use cases and the Word2Vec example. The tables, figures and description in the Use Case section of the paper reference these.
The sections of the paper describing the implementation summarize the functions implemented by the code in the code sub-directory.
