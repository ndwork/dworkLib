
function out = smoothData( in, varargin )
  % out = smoothData( in [, N, 'gaussian', sigma, 'op', op ] );
  %
  % Inputs:
  % in - either a 2D array (representing an image) or a 3D array (representing a color image)
  %   if op is 'notransp' then in should be the image
  %   if op is 'transp' then adjoint is applied to in
  %
  % Optional Inputs:
  % N - either a scalar or a 1D array specifying the size of each dimension of the
  %     smoothing kernel
  % sigma - by default, the smoothing operator is a box car filter.
  %         if sigma is provided, then the smoothing operator is a Gaussian filter
  % op - either 'notransp' (default, meaning filter), or 'transp' (meaning return the
  %      adjoint of the filter)
  %
  % Written by Nicholas - Copyright 2021
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 1
    disp( 'Usage:  out = smoothData( in [, N, ''gaussian'', sigma, ''op'', op ] );' );
    return
  end

  defaultSigma = 0;
  p = inputParser;
  p.addOptional( 'N', [] );
  p.addParameter( 'gaussian', defaultSigma, @isnumeric );
  p.addParameter( 'op', 'notransp', @(x) true );
  p.parse( varargin{:} );
  N = p.Results.N;
  sigma = p.Results.gaussian;
  op = p.Results.op;

  if numel( N ) == 0
    N = ceil( max( 5*sigma, 3 ) );
    N( mod( N, 2 ) == 0 ) = N( mod( N, 2 ) == 0 ) + 1;
  end

  numDims = ndims( in );
  if numel( N ) == 1, N = N * ones( 1, numDims ); end
  if numel( N ) ~= numDims
    error( 'N must have one element or the same number of elements as dimensions of in' );
  end

  if min( mod( N, 2 ) ) == 0, error( 'All values of N must be odd' ); end

  if sigma > 0
    % Gaussian filter with standard deviation sigma
    h = makeGaussFilter( N, sigma );
  else
    % Box car average (or mean) filter
    h = ones( N ) / sum( N(:) );
  end

  if strcmp( op, 'transp' )

    nH = numel( h );   %#ok<NASGU>
    halfN = ceil( N / 2 );   %#ok<NASGU>

    cmd = '[ ';
    for indx = 1 : numel(N)-1
      cmd = [ cmd, 'I', num2str(indx), ', ' ];   %#ok<AGROW>
    end
    cmd = [ cmd, 'I', num2str(numel(N)), ' ] = ind2sub( N, (1:nH)'' );' ];
    eval( cmd );

    hShifts = [];
    cmd = 'hShifts = [ ';
    for indx = 1 : numel(N)-1
      cmd = [ cmd, 'I', num2str(indx), ' - halfN(', num2str(indx), '), ' ];   %#ok<AGROW>
    end
    cmd = [ cmd, 'I', num2str(numel(N)), ' - halfN(', num2str(numel(N)), ') ];' ];
    eval( cmd );

    nH = numel( h );
    segLength = 100;
    nSegs = ceil( nH / segLength );
    cmd = [ 'segCells = cell( ', repmat( '1,', [ 1 numel(N) ] ), 'nSegs );' ];
    eval( cmd );
    parfor segIndx = 1 : nSegs
      tmp = zeros( size(in) );
      for hIndx = (segIndx-1) * segLength + 1 : min( segIndx * segLength, nH )
        shifted = shiftImg( in, hShifts(hIndx,:) );   %#ok<PFBNS>
        tmp = tmp + h( hIndx ) * shifted; %#ok<PFBNS>
      end
      segCells{ segIndx } = tmp;
    end
    out = cell2mat( segCells );
    out = sum( out, numel(N) + 1 );

  else

    out = convn( in, h, 'same' );  % ok to do convolution since h is symmetric

  end
end
