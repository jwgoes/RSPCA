% Initializes the PCA eigenvectors and eigenvalues to a random 
% kk-dimensional subspace.
%
% kk - the dimension of the subspace which we seek
% dd - the dimension of the overall space
%
function [ vectors, values ] = pca_initialize_random_orthogonal( kk, dd )

vectors = random( 'normal', 0, 1, dd, kk );
vectors = vectors * pinv( sqrtm( vectors' * vectors ) );
vectors = real( vectors );

values = ones( 1, kk ) / ( ( dd - kk ) * ( kk) );
