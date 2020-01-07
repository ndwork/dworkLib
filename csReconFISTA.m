
function [recon,oValues] = csReconFISTA( samples, lambda, varargin )
  % recon = csReconFISTA( samples, lambda [, 'debug', debug, 'nIter', nIter, ...
  %   'polish', polish, 'printEvery', printEvery, 'verbose', verbose, ...
  %   'waveletType', waveletType ] )
  %
  % This routine minimizes 0.5 * || Ax - b ||_2^2 + lambda || W x ||_1
  %   where A is sampleMask * Fourier Transform * real part, and
  %   W is the wavelet transform.
  %
  % Inputs:
  % samples - a 2D array that is zero wherever a sample wasn't acquired
  % lambda - regularization parameter
  %
  % Optional Inputs:
  % debug - if true, reduces the default number of iterations to 30 and forces verbose
  %         statements during optimization
  % nIter - the number of iterations that FISTA will perform (default is 100)
  % polish - if set to true, adds a polishing step (default is false)
  % printEvery - FISTA prints a verbose statement every printEvery iterations
  % verbose - if true, prints informative statements
  % waveletType - either 'Deaubechies' for Deaubechies-4 (default) or 'Haar'
  %
  % Written by Nicholas Dwork - Copyright 2017
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular purpose.

  wavSplit = zeros(4);  wavSplit(1,1) = 1;
  %wavSplit = [1 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0;];
  %wavSplit = [ 1 0; 0 0; ];

  p = inputParser;
  p.addParameter( 'checkAdjoints', false, @islogical );
  p.addParameter( 'debug', false, @(x) isnumeric(x) || islogical(x) );
  p.addParameter( 'nIter', [], @ispositive );
  p.addParameter( 'polish', false, @(x) isnumeric(x) || islogical(x) );
  p.addParameter( 'printEvery', 1, @ispositive );
  p.addParameter( 'wavSplit', wavSplit, @isnumeric );
  p.addParameter( 'verbose', false, @(x) isnumeric(x) || islogical(x) );
  p.addParameter( 'waveletType', 'Deaubechies', @(x) true );
  p.parse( varargin{:} );
  checkAdjoints = p.Results.checkAdjoints;
  debug = p.Results.debug;
  nIter = p.Results.nIter;
  polish = p.Results.polish;
  printEvery = p.Results.printEvery;
  waveletType = p.Results.waveletType;
  verbose = p.Results.verbose;

  if numel( nIter ) == 0
    if debug == true
      nIter = 30;
    else
      nIter = 100;
    end
  end

  nSamples = numel( samples );
  M = ( samples ~= 0 );

  % RI = [ Re; Im; ]
  % A = M F RI , A' = (RI)' * F' * M
  % A' * A = (RI)' * F' * M * F * RI
  % gGrad = A'*A*x - A'*b;

  function Fx = F( x )
    Fx = 1/sqrt(nSamples) .* fftc( x );
  end

  function Fadjy = Fadj( y )
    Fadjy = sqrt(nSamples) * ifftc( y );
  end

  function out = A( x )
    Fx = F( x );
    out = Fx( M == 1 );
  end

  function out = Aadj( y )
    MTy = zeros( size( samples ) );
    MTy( M==1 ) = y;
    out = Fadj( MTy );
  end

  if checkAdjoints == true
    % Variable used during debugging of this routine
    checkResult = csReconFISTA_checkAdjoints( samples, @F, @Fadj, @A, @Aadj );
    if checkResult ~= false, disp( 'Adjoints test passed' ); end
  end

  b = samples( M == 1 );

  function out = g( x )
    diff = A( x ) - b;
    out = 0.5 * norm( diff(:), 2 ).^2;
  end

  Aadjb = Aadj( b );
  function out = gGrad( x )
    out = Aadj( A( x ) ) - Aadjb;
  end

  if strcmp( waveletType, 'Deaubechies' )
    W = @(x) wtDeaubechies2( x, wavSplit );
    WT = @(y) iwtDeaubechies2( y, wavSplit );

  elseif strcmp( waveletType, 'Haar' )
    W = @(x) wtHaar2( x, wavSplit );
    WT = @(y) iwtHaar2( y, wavSplit );

  else
    error( 'Unrecognized wavelet type' );
  end

  proxth = @(x,t) WT( proxL1Complex( W(x), t*lambda ) );
    % The proximal operator of || W x ||_1 was determined using
    % Vandenberghe's notes from EE 236C; slide of "The proximal mapping" entitled
    % "Composition with Affine Mapping"

  function out = h( x )
    Wx = W(x);
    out = sum( abs( Wx(:) ) );
  end

  x0 = ifftc( samples );  % Could possibly make this better with inverse gridding

  t = 1 / numel( samples );
  if debug
    %[recon,oValues] = fista( x0, @g, @gGrad, proxth, 'h', @h, 'verbose', verbose );   %#ok<ASGLU>
    [recon,oValues] = fista_wLS( x0, @g, @gGrad, proxth, 'h', @h, ...
      't0', t, 'N', nIter, 'verbose', true, 'printEvery', printEvery );                                                                     %#ok<ASGLU>
  else
    %recon = fista( x0, @g, @gGrad, proxth );   %#ok<UNRCH>
    recon = fista_wLS( x0, @g, @gGrad, proxth, 't0', t, 'N', nIter, ...
      'verbose', verbose, 'printEvery', printEvery );
  end

  if polish
    maskW = proxL1Complex( W(recon), t*lambda ) ~= 0;
    proxth = @(x,t) WT( maskW .* W(x) );
    if debug
      [recon,oValues2] = fista_wLS( recon, @g, gGrad, proxth, 'h', @h, ...
        'N', nIter, 't0', t, 'verbose', verbose, 'printEvery', printEvery );   %#ok<ASGLU>
    else
      recon = fista_wLS( recon, @g, gGrad, proxth, 'h', @h, 'N', nIter, 't0', t, ...
        'verbose', verbose, 'printEvery', printEvery );
    end
  end
end


function out = csReconFISTA_checkAdjoints( samples, F, Fadj, A, Aadj )
  % Check to make sure that Fadj is the adjoint of F
  randSamples = rand( size(samples) ) + 1i * rand( size(samples) );

  if checkAdjoint( randSamples, F, Fadj ) ~= true
    error( 'Fadj is not the transpose of F' );
  end

  if checkAdjoint( randSamples, A, Aadj ) ~= true
    error( 'Aadj is not the transpose of A' );
  end

  out = true;
end

