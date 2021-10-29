
function F = iGrid_2D( data, traj, varargin )
  % F = iGrid_2D( data, traj, [ 'alpha', alpha, 'W', W, 'nC', nC ] )
  %
  % Inverse Gridding based on EE369C notes by John Pauly and Beatty et. al., IEEE TMI, 2005
  % Definitions and details according to http://nicholasdwork.com/tutorials/dworkGridding.pdf
  %
  % Inputs
  %   data is a 2D array specifying the volume to be encoded
  %   traj is a Mx2 array specifying the k-space trajectory.
  %     The first/second column is ky/kx
  %     The units are normalized to [-0.5,0.5).
  %
  % Optional Inputs:
  %   alpha is the oversampling factor >= 1  (default of 1.5)
  %   W is the window width in pixels  (default of 8)
  %   nC is the number of points to sample the convolution kernel (default of 500)
  %
  % Output:
  %   F the estimates of the Fourier coefficients along the trajectory
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 1
    disp( 'Usage:  F = iGrid_2D( data, traj, [ ''alpha'', alpha, ''W'', W, ''nC'', nC ] ) ');
    return;
  end

  checknum = @(x) numel(x) == 0 || ( isnumeric(x) && isscalar(x) && (x >= 1) );
  p = inputParser;
  p.addParameter( 'alpha', [], checknum );
  p.addParameter( 'W', [], checknum );
  p.addParameter( 'nC', [], checknum );
  p.parse( varargin{:} );
  alpha = p.Results.alpha;
  W = p.Results.W;
  nC = p.Results.nC;

  [ Ny, Nx, ~ ] = size( data );

  if numel( alpha ) == 0, alphaY = 1.5; alphaX = 1.5;
  elseif numel( alpha ) == 1, alphaY = alpha; alphaX = alpha;
  elseif numel( alpha ) == 2, alphaY = alpha(1);  alphaX = alpha(2);
  end

  if ( alphaY ~= 1 ) || ( alphaX ~= 1 )
    sData = size( data );
    sData( 1 : 2 ) = ceil( [ alphaY alphaX ] .* sData( 1 : 2 ) );
    data = padData( data, sData );
  end

  if ~isreal( traj ), traj = [ real(traj(:))  imag(traj(:)) ]; end

  % Make the Kaiser Bessel convolution kernel
  Gy = Ny;
  [kCy,Cy,cY] = makeKbKernel( Gy, Ny, 'alpha', alphaY, 'W', W, 'nC', nC );
  Gx = Nx;
  [kCx,Cx,cX] = makeKbKernel( Gx, Nx, 'alpha', alphaX, 'W', W, 'nC', nC );

  % Pre-emphasize the image
  preEmphasized = bsxfun( @rdivide, data, cY * transpose( cX ) );

  % Perform an fft
  fftData = fftshift2( fft2( ifftshift2( preEmphasized ) ) ) / ...
    ( size( preEmphasized, 1 ) * size( preEmphasized, 2 ) );

  % Perform a circular convolution
  sData = size( fftData );
  F = applyC_2D( fftData, sData(1:2), traj, kCy, kCx, Cy, Cx );
end

