"""
    Module for Principal Component Analysis.
    
    Features:
    * pca and kernel pca
    * pca through singular value decomposition (SVD)

    Author: Alexis Mignon (copyright)
    Date: 10/01/2012
    e-mail: alexis.mignon@gmail.com
    
    Modified by: Sibt ul Hussain (sibt.ul.hussain dot gmail.com), date=01-01-2013
"""

import numpy as np
from scipy.sparse.linalg.eigen.arpack import eigsh
from scipy.linalg import eigh

def full_pca(data):
    """
        Performs the complete eigen decomposition of
        the covariance matrix.
        
        arguments:
        * data: 2D numpy array where each row is a sample and
            each column a feature.
        
        return:
        * w: the eigen values of the covariance matrix sorted in from 
              highest to lowest.
        * u: the corresponding eigen vectors. u[:,i] is the vector
             corresponding to w[i]
             
        Notes: If you want to compute only a few number of principal
               components, you should consider using 'pca' or 'svd_pca'.
    """
    cov = np.cov(data.T)
    w, u = eigh(cov, overwrite_a=True)
    return w[::-1], u[:, ::-1]

def pca(data, k):
    """
        Performs the eigen decomposition of the covariance matrix.
        
        arguments:
        * data: 2D numpy array where each row is a sample and
                each column a feature.
        * k: number of principal components to keep.
        
        return:
        * w: the eigen values of the covariance matrix sorted in from 
              highest to lowest.
        * u: the corresponding eigen vectors. u[:,i] is the vector
             corresponding to w[i]
             
        Notes: If the number of samples is much smaller than the number
               of features, you should consider the use of 'svd_pca'.
    """
    cov = np.cov(data.T)
    w, u = eigsh(cov, k=k, which='LA')
    return w[::-1], u[:, ::-1]

def extern_pca(data, k):
    """
        Performs the eigen decomposition of the covariance matrix based
        on the eigen decomposition of the exterior product matrix.
        
        
        arguments:
        * data: 2D numpy array where each row is a sample and
                each column a feature.
        * k: number of principal components to keep.
        
        return:
        * w: the eigen values of the covariance matrix sorted in from 
              highest to lowest.
        * u: the corresponding eigen vectors. u[:,i] is the vector
             corresponding to w[i]
             
        Notes: This function computes PCA, based on the exterior product
               matrix (C = X*X.T/(n-1)) instead of the covariance matrix
               (C = X.T*X) and uses relations based of the singular
               value decomposition to compute the corresponding the
               final eigen vectors. While this can be much faster when 
               the number of samples is much smaller than the number
               of features, it can lead to loss of precisions.
               
               The (centered) data matrix X can be decomposed as:
                    X.T = U * S * v.T
               On computes the eigen decomposition of :
                    X * X.T = v*S^2*v.T
               and the eigen vectors of the covariance matrix are
               computed as :
                    U = X.T * v * S^(-1)
    """
    data_m = data - data.mean(0)
    K = np.dot(data_m, data_m.T)
    w, v = eigsh(K, k=k, which='LA')
    U = np.dot(data.T, v / np.sqrt(w))
    return w[::-1] / (len(data) - 1), U[:, ::-1]

def full_kpca(data):
    """
        Performs the complete eigen decomposition of a kernel matrix.
        
        arguments:
        * data: 2D numpy array representing the symmetric kernel matrix.
        
        return:
        * w: the eigen values of the covariance matrix sorted in from 
              highest to lowest.
        * u: the corresponding eigen vectors. u[:,i] is the vector
             corresponding to w[i]
             
        Notes: If you want to compute only a few number of principal
               components, you should consider using 'kpca'.
    """
    w, u = eigh(data, overwrite_a=True)
    return w[::-1], u[:, ::-1]


def kpca(data, k):
    """
        Performs the eigen decomposition of the kernel matrix.
        
        arguments:
        * data: 2D numpy array representing the symmetric kernel matrix.
        * k: number of principal components to keep.
        
        return:
        * w: the eigen values of the covariance matrix sorted in from 
              highest to lowest.
        * u: the corresponding eigen vectors. u[:,i] is the vector
             corresponding to w[i]
             
        Notes: If you want to perform the full decomposition, consider 
               using 'full_kpca' instead.
    """
    w, u = eigsh(data, k=k, which='LA')
    return w[::-1], u[:, ::-1]

class PCA(object):
    """
        PCA object to perform Principal Component Analysis.
    """
    def __init__(self, k=None, kernel=False, extern=False):
        """
            Constructor.
            
            arguments:
            * k: number of principal components to compute. 'None'
                 (default) means that all components are computed.
            * kernel: perform PCA on kernel matrices (default is False)
            * extern: use extern product to perform PCA (default is 
                   False). Use this option when the number of samples
                   is much smaller than the number of features.
        """
    
        self.k = k
        self.kernel = kernel
        self.extern = extern
    
    def fit(self, X):
        """
            Performs PCA on the data array X.
            arguments:
            * X: 2D numpy array. In case the array represents a kernel
                 matrix, X should be symmetric. Otherwise each row
                 represents a sample and each column represents a
                 feature.
        """
        if self.k is None :
            if self.kernel :
                pca_func = full_kpca
            else :
                pca_func = full_pca
            self.eigen_values_, self.eigen_vectors_ = pca_func(X)
        else :
            if self.kernel :
                pca_func = kpca
            elif self.extern :
                pca_func = extern_pca
            else :
                pca_func = pca
            self.eigen_values_, self.eigen_vectors_ = pca_func(X, self.k)

        if not self.kernel :
            self.mean = X.mean(0)
        return self
        
    def transform(self, X, whiten=False, ncomp=None):
        """
            Project data on the principal components. If the whitening
            option is used, components will be normalized so that they
            have the same contribution.
            
            arguments:
            * X: 2D numpy array of data to project.
            * whiten: (default is False) all components are normalized
                so that they have the same contribution.
            * ncomp : number of components to use during projection... (default use K components)
                works only for non-kernel case
                
            returns:
            * prX : projection of X on the principal components.
            
            Notes: In the case of Kernel PCA, X[i] represents the value
               of the kernel between sample i and the j-th sample used
               at train time. Thus, if fit was called with a NxN kernel
               matrix, X should be a MxN matrix.
               
               The projection in the kernel case is made to be equivalent
               to the projection in the linear case.
               
                   X.T = U * S * v.T
                   C = 1/(N-1) * X.T * X
                   X.T * X = U*S^2*U.T
                   K = X * X.T = v*S^2*v.T
                   
                   U = X.T * v * S^(-1)
               
               The projection with PCA is :
                   X' = X * U
                   X' = X * X.T * v * S^(-1)
                   X' = K * v * S^(-1)
                   
               For whiten PCA :
                   X' = X * U * S^(-1) * sqrt(N-1)
                   X' = X * X.T * v * S^(-1) * S^(-1) * sqrt(N-1)
                   X' = K * S^(-2) * sqrt(N-1)
        """
        if ncomp is None:
            ncomp = self.k  # use all components...
        
        
        if self.kernel :
            pr = np.dot(X, self.eigen_vectors_)
            if whiten :
                pr /= self.eigen_values_ / np.sqrt(X.shape[0] - 1)
            else :
                pr /= np.sqrt(self.eigen_values_)
        else :
            pr = np.dot(X - self.mean, self.eigen_vectors_[:, :ncomp])
            if whiten:
                pr /= np.sqrt(self.eigen_values_[:ncomp])
        return pr
