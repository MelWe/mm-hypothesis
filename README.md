# mm-hypothesis
This is the code repository for "Heuristic Framework for Testing the Multi-Manifold Hypothesis" by F. Patricia Medina, Linda Ness, Melanie Webr and Karamatou Yacoubou Djima.

-----------

### Abstract
When analyzing empirical data, we often find that global linear models overestimate the number of parameters required. In such cases, we may ask whether the data lies on or near a manifold or a set of manifolds (a so-called multi-manifold) of lower dimension than the ambient space. This question can be phrased as a (multi-) manifold hypothesis. The identification of such intrinsic multiscale features
is a cornerstone of data analysis and representation, and has given rise to a large body of work on manifold learning. In this work, we review key results on multi- scale data analysis and intrinsic dimension followed by the introduction of a heuris- tic, multiscale framework for testing the multi-manifold hypothesis. Our method implements a hypothesis test on a set of spline-interpolated manifolds constructed from variance-based intrinsic dimensions. The workflow is suitable for empirical data analysis as we demonstrate on two use cases.

*This directory contains an implementation of the multi-manifold hypothesis testing framework as well as the data sets and results for both use cases.*

### Code 
The folder *codes* contains implementations for the following subroutines, described in detail in the paper:
1) _MMMethodcode_ contains a main calling program: _testMMmethodNov1.m_ which specifies and executes the use cases in the paper (Sphere-Line and LiDAR). The other programs are invoked to construct the dyadic linear multi-manifold model for the data and the hypothesis testing distribution.
2) _Vidimcode_ contains implementations for variance-based intrinsic dimension. There are two MATLAB programs: _testvariancebaseddimLIDAR_ and _computevidim_. The latter program reorganizes the main calling program so it can be called more clearly by the _MMMethodcode_.
3) _Vidimimunpackedcode_ contains the individual procedures used to compute intrinsic dimension. It also contains the idimtests program which is used to compute idim related statistics in a structure idimstats.
4) _GMSTidim.m_ estimates the intrinsic dimension at each point using a variant of the geodesic minimal spanning tree (GMST) to estimate the intrinsic dimension and entropy of the manifold on which the data lie. For global GMST we use a toolbox for dimensionality reduction by Laurens van der Maaten (https://lvdmaaten.github.io/drtoolbox/). The scripts _testvariancebaseddimandTV_ and _testvariancebaseddimandTVLIDAR_ are slight modifications of original code to work on neighborhood created using k nearest neighbors instead of epsilon balls.

### Datasets
The folder *usecases* contains subfolders for three use cases with the raw data files in _/data_:
1) LiDAR: The dataset contains (x,y,z)-coordinates for LiDAR imaging data of the Golden Gate bridge. The files contains a preprocessed version of USGS LiDAR data, the preprocessing is described in the paper.
2) Sphere-Line: We folder contains MATLAB code for generating samples from a sphere-line configuration.
3) W2V: Pre-processed Word2Vec data in MATLAB format.

### Results 
The results directory contains the file outputs by the execution of all codes. This includes the results underlying the figures and tables in the paper.






