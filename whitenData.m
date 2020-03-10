
function out = whitenData( data, covMatrix, varargin )
  % out = whitenData( data, covMatrix [, type ] )
  %
  % Assumes data is zero mean
  %
  % Inputs:
  % data - a multidimensional array; it is assumed that the last dimension
  %   is the random vector dimension
  % covMatrix - a covariance matrix
  %
  % Optional Inputs:
  % type - the type of whitening to perform.
  %   'chol' for Cholesky or 'pca' for Principal Component Analysis
  %
  % Outputs:
  % out - an array of whitened vectors of size data
  %
  % Written by Nicholas Dwork - Copyright 2020
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 2
    disp( 'Usage:  out = whitenData( data, covMatrix [, type ] )' );
    return
  end

  p = inputParser;
  p.addOptional( 'type', 'chol', @(x) true );
  p.parse( varargin{:} );
  type = p.Results.type;

  if size( data, ndims( data ) ) ~= size( covMatrix, 1 ) || ...
     size( data, ndims( data ) ) ~= size( covMatrix, 2 )
     error( 'Dimension mismatch' );
  end

  sData = size( data );
  reshaped = transpose( reshape( data, [ prod( sData(1:end-1) ) sData(end) ] ) );
  
  if strcmp( type, 'chol' )
    L = chol( covMatrix, 'lower' );
    whitened = L \ reshaped;

  elseif strcmp( type, 'pca' )
    % PCA whitening
    %  https://stats.stackexchange.com/questions/95806/how-to-whiten-the-data-using-principal-component-analysis
    [u,s,~] = svd( covMatrix, 'econ' );
    whitened = sqrt(s) \ u' * reshaped;

  else
    error( 'Unknown whitening type' );

  end

  out = reshape( whitened, sData );
end
