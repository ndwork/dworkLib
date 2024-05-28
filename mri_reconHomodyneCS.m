
function recon =  mri_reconHomodyneCS( kData, sFSR, varargin )
  % recon =  mri_reconHomodyneCS( kData, sFSR [, 'lambda', lambda, 'wavSplit', wavSplit ] )
  %
  % Inputs
  % kData - a 2D array with value equal to 0 wherever data was not collected or
  %         a 3D array where the third dimension is the coils
  %
  % Outputs
  % recon - a 2D array representing the reconstructed image
  %
  % Written by Nicholas Dwork - Copyright 2024
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular purpose.

  p = inputParser;
  p.addParameter( 'doCheckAdjoints', false );
  p.addParameter( 'epsilon', [], @isnonnegative );
  p.addParameter( 'verbose', true );
  p.addParameter( 'wavSplit', [], @isnumeric );
  p.parse( varargin{:} );
  doCheckAdjoints = p.Results.doCheckAdjoints;
  epsilon = p.Results.epsilon;
  verbose = p.Results.verbose;
  wavSplit = p.Results.wavSplit;

  sImg = size( kData, [ 1 2 ] );
  if numel( wavSplit ) == 0
    wavSplit = makeWavSplit( sImg );
  end

  sampleMask = ( kData ~= 0 );
  nu = ceil( ( sImg(1) + 1 ) / 2 ) + round( sFSR(1)/2 ) - 1;  % the last index of the partial Fourier data
  sTopPortion = [ nu sImg(2) ];

  [~,phaseImg] = mri_reconPFHomodyne( kData, sFSR );
  phases = angle( phaseImg );

  function out = applyA(in, op)
    if nargin < 2 || strcmp( op, 'notransp' )
      tmp = zeros( sImg );
      tmp( 1 : nu, : ) = reshape( in, sTopPortion );
      Pin = mri_reconPFHomodyne( tmp, sFSR, 'phases', phases );
      out = wtDaubechies2( Pin, wavSplit );
    else
      in = reshape( in, sImg );
      WHin = iwtDaubechies2( in, wavSplit );
      out = mri_reconPFHomodyne( WHin, sFSR, 'phases', phases, 'op', 'transp' );
      out = out( 1 : nu, : );
    end
    out = out(:);
  end

  f0 = kData( 1 : nu, : );

  if doCheckAdjoints == true
    innerProd = @(x,y) real( dotP( x, y ) );
    [chkA,errA] = checkAdjoint( f0(:), @applyA, 'innerProd', innerProd );
    if chkA == false
      error([ 'Adjoint check of A failed with error ', num2str(errA) ]);
    end
  end

  function out = g( in )
    out = norm( in(:), 1 );
  end

  b = kData( sampleMask == 1 );
  function out = h( in )  % Indicator function for data consistency
    normDiff = norm( in( sampleMask( 1 : nu, : ) == 1 ) - b, 2 );
    if normDiff > epsilon, out = inf; else, out = 0; end
  end

  function out = proxh( in, t )   %#ok<INUSD>
    if epsilon == 0
      out = in;
      out( sampleMask( 1 : nu, : ) == 1 ) = b;
    else
      normDiff = norm( in( sampleMask( 1 : nu, : ) == 1 ) - b, 2 );
      out = in;
      if normDiff > epsilon
        out( sampleMask( 1 : nu, : ) == 1 ) = ( in( sampleMask( 1 : nu, : ) == 1 ) - b ) * ( epsilon / normDiff ) + b;
      end
    end
  end

  function out = proxg( in, t )
    out = softThresh( in, t );
  end

  function out = proxgConj( in, sigma )
    out = proxConj( @proxg, in, sigma );
  end

  % Use PDHG to solve minimize || A f ||_1 subject to || D f - b ||_2 <= epsilon

  f0 = kData( 1 : nu, : );
  normA = powerIteration( @applyA, rand( size( f0 ) ) );

  tau = 1d-2 / normA;

  [fStar,objValues,relDiffs] = pdhg( f0(:), @proxh, @proxgConj, tau, 'A', @applyA, 'f', @h, 'g', @g, ...
    'N', 10000, 'normA', normA, 'printEvery', 100, 'verbose', verbose );   %#ok<ASGLU>

  kRecon = zeros( sImg );
  kRecon( 1 : nu, : ) = reshape( fStar, sTopPortion );
  recon = mri_reconPFHomodyne( kRecon, sFSR, 'phases', phases );
end
