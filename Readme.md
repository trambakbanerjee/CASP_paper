What is CASP?
======

CASP (Coordinate-wise Adaptive Shrinkage Prediction) is a novel shrinkage rule for prediction in a high-dimensional non-exchangeable hierarchical Gaussian model with an unknown spiked covariance structure. CASP uses results on the behavior of eigenvalues and eigenvectors of high-dimensional sample covariance matrix ([2], [3], [4]) to develop a bias-correction principle that leads to an efficient approach for evaluating the Bayes predictors corresponding to popular loss functions such as
quadratic, generalized absolute, and linex losses.

How to use this repository?
----------

This repository holds the scripts that reproduce the analysis in the paper [1]. Send me an email if anything does not work as expected. If you are looking for the associated R-package for CASP, please [visit this page](https://github.com/trambakbanerjee/casp#casp).

References
=======
[1.] Improved Shrinkage Prediction under a Spiked Covariance Structure _(under review)_     
Banerjee, T., Mukherjee, G. and Paul, D.

[2.] Baik, J. and J. W. Silverstein (2006). Eigenvalues of large sample covariance matrices of spiked population models. Journal of Multivariate Analysis 97(6), 1382–1408.

[3.] Onatski, A. (2012). Asymptotics of the principal components estimator of large factor models with weakly influential factors. Journal of Econometrics 168(2), 244–258.

[4.] Paul, D. (2007). Asymptotics of sample eigenstructure for a large dimensional spiked covariance model. Statistica Sinica, 1617–1642.
