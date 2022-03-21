
function out = mrs_reconSense( kData, sMaps, varargin )
  % out = mrs_reconSense( kData, sMaps [, 'kTraj', 'kTraj', 'mu', mu, 'sImg', sImg ] )
  %
  % Performs model based reconstruction for MR spectroscopy
  %
  % Inputs:
  % kData - array of Fourier values.
  %   For Cartesian trajectory, it is of size M x N x nCoils x N1 x N2 x ... x NF
  % sMaps - sensitivity maps of size M x N x nCoils
  %
  % Optional Inputs:
  % kTraj - nTraj x 2 array of (ky,kx)
  % mu - tikhonov regularization parameter
  %
  % Outputs:
  % out - reconstructed spectrums of size M x N x nF
  %
  % Written by Nicholas Dwork - Copyright 2021
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular purpose.

  if nargin < 1
    disp([ 'Usage:  out = mrs_reconSense( kData, sMaps [, ', ...
      ' ''kTraj'', kTraj, ''mu'', mu, ''sImg'', sImg ] )' ]);
    if nargout > 0, out = []; end
    return
  end

  p = inputParser;
  p.addParameter( 'kTraj', [], @isnumeric );
  p.addParameter( 'mu', 0, @ispositive );
  p.addParameter( 'sImg', [], @ispositive );
  p.parse( varargin{:} );
  kTraj = p.Results.kTraj;
  mu = p.Results.mu;
  sImg = p.Results.sImg;


  if numel( kTraj ) == 0
    [ M, N, nCoils, nF ] = size( kData );

    sKData = size( kData );  Ns = sKData(4:end);
    kData = reshape( kData, [ sKData(1:3) prod( Ns ) ] );

  else
    M = sImg(1);
    N = sImg(2);
    sKData = size( kData );
    nCoils = sKData( 2 );
    Ns = sKData(3:end);
    nF = prod( Ns );
    kData = reshape( kData, [ sKData(1:2) nF ] );

  end


  function out = applyS( x, op )
    if nargin < 2  ||  strcmp( op, 'notransp' )
      % x is size M x N x nSlices x nF
      x = reshape( x, [ M N 1 nF ] );
      Sx = bsxfun( @times, x, sMaps );

      if numel( kTraj ) == 0
        out = fftshift2( ufft2( ifftshift2( Sx ) ) );
      else
        out = iGrid_2D( Sx, kTraj );
      end

    elseif strcmp( op, 'transp' )
      % x is size M x N x nSlices x nCoils x nF;
      if numel( kTraj ) == 0
        tmp = fftshift2( uifft2( ifftshift2( x ) ) );
      else
        tmp = iGridT_2D( x, kTraj, sImg );
      end

      tmp = reshape( tmp, [ M N nCoils nF ] );
      out = bsxfun( @times, tmp, conj( sMaps ) );
      out = permute( sum( out, 3 ), [ 1 2 4 3 ] );
    else

      error( 'Unrecognized operation for applyS' );
    end
  end

  doCheckAdjoint = false;
  if doCheckAdjoint == true
    x0 = rand( M, N, nSlices, nF );
    [ checkS, errCheckS ] = checkAdjoint( x0, @applyS );
    if checkS == false
      error([ 'applyS adjoint failed with error: ', num2str( errCheckS ) ]);
    else
      disp([ 'applyS adjoint passed with error: ', num2str( errCheckS ) ]);
    end
  end

  nData = numel( kData(:) );
  function out = g( in )
    Sin = applyS( in );
    out = 0.5 * norm( Sin(:) - kData(:) )^2 + 0.5 * mu / nData * norm( in(:) )^2;
  end

  STb = applyS( kData, 'transp' );
  function out = gGrad( in )
    out = applyS( applyS( in ), 'transp' ) + mu / nData * in - STb;
  end

  x0 = mrs_reconRefPeak( kData, 'kTraj', kTraj, 'sImg', sImg );

  %t = 1d-2;
  t = 1d-3;
  tol = 1d-3;
  [out,objectiveValues,relDiffs] = gradDescent( x0, @gGrad, 'g', @g, 't', t, ...
    'tol', tol, 'N', 10000, 'useLineSearch', true, 'beta', 0.5, ...
    'printEvery', 1, 'verbose', true );   %#ok<ASGLU>

  out = reshape( out, [ M N Ns ] );
end

