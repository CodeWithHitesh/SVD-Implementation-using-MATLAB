# SVD-Implementation-using-MATLAB
Singular Value Decomposition (SVD) is used widely in signal processing. Noise reduction and image compression are some of the applications of SVD. I’ll be using it for reducing noise in an image.  

The singular value decomposition takes an m \ x \ n matrix A and decomposes it into A = U \Sigma V^{T}. \Sigma is a diagonal matrix that contains the singular values of A. U and V are orthogonal matrices where U is an m \ x \ m matrix and V is an n \ x \ n matrix.  

Though we have inbuilt functions in MATLAB for finding the SVD of an image but I’ve implemented my own function whose output is same as that of the inbuilt function. For comparing our output, we will be using Peak signal-to-noise ratio, often abbreviated PSNR, which is a term for the ratio between the maximum possible power of a signal and the power of corrupting noise that affects the fidelity of its representation.


## Output
We are using PSNR for comparing the output.

#### For original and noisy image: 

PSNR value - 23.2502

MSE - 307.6494



#### For original and output image: 

PSNR value - 30.3783

MSE - 59.6004

