
function out = diChromInterp( highResImg, lowResImg, sigma, varargin )
  % out = diChromInterp( highResImg, lowResImg, sigma [ 'contrasts', true/false, ...
  %   'lambda', labmda, 'N', N, 'optAlg', optAlg, 'showScale', showScale, ...
  %   'verbose', verbose ] )
  %
  % This algorithm was published as "Di-chromatic Interpolation of Magnetic Resonance Metabolic
  % Imagery" in the Magnetic Resonance Materials in Physics, Biology, and Medicine Journal.
  % http://link.springer.com/article/10.1007/s10334-020-00903-y
  %
  % Inputs:
  % highResImg - a 2D array representing a grayscale of the high resolution image
  % lowResImg - a 2D array representing a grayscale of the low resolution image
  % sigma - the standard deviation of the Gaussian blur kernel in pixels
  %
  % Optional Inputs:
  % contrasts - if true (default) then use both highResImg and -highResImg for
  %   the reference image and determine the best result.  If false, then only use
  %   highResImg.
  % lambda - tikhonov regularization parameter
  % N - number of iterations for fista
  % optAlg - either 'projSubgrad', 'fista', or 'fista_wLS'
  % showScale - the scale factor for displayed images if verbose is true
  % verbose - show images
  %
  % Outputs:
  % srPyr - the super resolved pyruvate image
  %
  % Written by Nicholas Dwork, Copyright 2019

  p = inputParser;
  p.addRequired( 'highResImg', @isnumeric );
  p.addRequired( 'lowResImg', @isnumeric );
  p.addParameter( 'contrasts', true, @(x) islogical(x) || isnumeric(x) );
  p.addParameter( 'lambda', 1d3, @isnumeric );  % Tikhonov regularization parameter
  p.addParameter( 'N', 60, @isnumeric );
  p.addParameter( 'optAlg', 'fista_wLS', @(x) true );
  p.addParameter( 'showScale', 10, @ispositive );
  p.addParameter( 'verbose', false, @(x) isnumeric(x) || islogical(x) );
  p.parse( highResImg, lowResImg, varargin{:} );
  contrasts = p.Results.contrasts;
  lambda = p.Results.lambda;
  N = p.Results.N;
  optAlg = p.Results.optAlg;
  showScale = p.Results.showScale;
  verbose = p.Results.verbose;

  maxH = max( highResImg(:) );     highResImg = highResImg / maxH;
  maxL = max( lowResImg(:) );      lowResImg = lowResImg / maxL;

  if contrasts == true

    diChromInterps = cell(2,1);  
    parfor i = 1 : 2
      if i == 1
        thisInterp = diChromInterp_contrast( highResImg, lowResImg, sigma, ...
          lambda, N, optAlg, verbose );
      else
        thisInterp = diChromInterp_contrast( 1-highResImg, lowResImg, sigma, ...
          lambda, N, optAlg, verbose );
      end
      diChromInterps{i} = thisInterp;
    end
    diChromPos = diChromInterps{1};
    diChromNeg = diChromInterps{2};

    sHighRes = size( highResImg );
    lowResInterp = imresize( lowResImg, sHighRes, 'bilinear' );
    posDiff = norm( lowResInterp(:) - diChromPos(:) );
    negDiff = norm( lowResInterp(:) - diChromNeg(:) );

    if negDiff < posDiff
      out = diChromNeg;
    else
      out = diChromPos;
    end

    if verbose == true
      figure;  imshowscale( [ diChromPos diChromNeg out ], showScale );
    end
  else

    out = diChromInterp_contrast( highResImg, lowResImg, sigma, ...
      lambda, N, optAlg, verbose );

    if verbose == true
      figure;  imshowscale( out, showScale );
    end

  end

  out = out * maxL;
end


function srImg = diChromInterp_contrast( highResImg, lowResImg, sigma, ...
	lambda, N, optAlg, verbose )

  lambda = lambda / numel( lowResImg );

  sHighRes = size( highResImg );
  sLowRes = size( lowResImg );
  if ismatrix( highResImg )
    srImg0 = imresize( lowResImg, sHighRes, 'bilinear' );
  else
    srImg0 = imresize3( lowResImg, sHighRes, 'linear' );
  end
  ws = srImg0; 

  function out = f( x )
    smoothX = smoothData( x, 'gaussian', sigma );
    if ismatrix( smoothX )
      out = imresize( smoothX, sLowRes, 'nearest' );
    else
      out = imresize3( smoothX, sLowRes, 'nearest' );
    end
  end

  D = @(x) computeGradient( x );

  function [fx, wDx] = applyA( x )
    fx = f( x );
    Dx = D( x );
    if numel( ws ) > 0 
      wDx = bsxfun( @times, Dx, lambda * ws );
    else
      wDx = lambda * Dx;
    end
  end

  upSampFactor = round( sHighRes ./ sLowRes );
  upSampShifts = round( upSampFactor / 2 ) - 1;
  function out = fAdj( x )
    upSamp = upsampleData( x, upSampFactor, 'sOut', sHighRes, 'S', upSampShifts );
    out = smoothData( upSamp, 'gaussian', sigma, 'op', 'transp' );
  end

  DT = @(x) computeGradient( x, 'transp' );

  function out = applyAdjA( x1, x2 )
    out1 = fAdj( x1 );
    wx2 = x2;
    if numel( ws ) > 0, wx2 = bsxfun( @times, wx2, ws ); end
    out2 = lambda * DT( wx2 );
    out = out1 + out2;
  end

  nLowRes = numel( lowResImg );
  function out = A( x, op )
    if nargin < 2 || strcmp( op, 'notransp' )
      x = reshape( x, sHighRes );
      [ fx, wDx ] = applyA( x );
      out = [ fx(:); wDx(:); ];
    else
      x1 = reshape( x(1:nLowRes), sLowRes );
      x2 = reshape( x(nLowRes+1:end), [ sHighRes ndims( highResImg ) ] );
      out = applyAdjA( x1, x2 );
      out = out(:);
    end
  end

  Dw = D( highResImg );
  if numel( ws ) > 0, Dw = bsxfun( @times, Dw, ws ); end
  b = [ lowResImg(:); lambda * Dw(:); ];

  function out = g( x )
    out = 0.5 * norm( A(x) - b, 2 ).^2;
  end

  ATb = A( b, 'transp' );
  function out = gGrad( x )
    out = A( A( x ), 'transp' ) - ATb;
  end

  function out = h( x )
    out = 0;
    if min( x ) < 0, out=Inf; end
  end

  %proxth = @(x,t) min( max( x, 0 ), 1 );
  proxth = @(x,t) max( x, 0 );

  %[check,err] = checkAdjoint( srImg0, @f, 'fAdj', @fAdj );
  %[check,err] = checkAdjoint( srImg0, D, 'fAdj', DT );
  %[check,err] = checkAdjoint( srImg0, @A );


  if strcmp( optAlg, 'projSubgrad' )
    t = 1.0;
    proj = @(x) min( max( x, 0 ), 1 );
    srImg = projSubgrad( srImg0(:), @gGrad, proj, 't', t, 'verbose', verbose );
  elseif strcmp( optAlg, 'fista' )
    [srImg,objectiveValues] = fista( srImg0(:), @g, @gGrad, proxth, ...
      'h', @h, 'N', N, 'verbose', verbose, 'printEvery', 10 );
  else
    [srImg,objectiveValues] = fista_wLS( srImg0(:), @g, @gGrad, proxth, ...
      'h', @h, 'N', N, 'verbose', verbose, 'printEvery', 10 );
  end
  srImg = reshape( srImg, sHighRes );


  function out = objectiveValue( x )
    [fx, Dx] = applyA( x );
    if numel( ws ) > 0, Dx = bsxfun( @times, Dx, ws ); end
    out1 = norm( fx(:) - lowResImg(:), 2 ).^2;
    out2 = lambda * norm( Dx(:) - Dw(:), 2 ).^2;
    out = 0.5 * ( out1 + out2 );
  end

  if verbose ~= 0 && numel( gcp( 'nocreate' ) ) == 0
    disp([ 'Final objectiveValue: ', num2str( objectiveValue( srImg ) ) ]);
    if exist( 'lsqrFlag', 'var' )
      disp([ 'lsqr Flag: ', num2str(lsqrFlag) ]);
    end
    srPyrNearest = imresize( lowResImg, sHighRes, 'nearest' );
    figure; imshowscale( srPyrNearest, 5 );  titlenice( 'Initial Guess' );
    figure; imshowscale( srImg, 5 );  titlenice( 'Super Res Pyruvate' );
    figure; imshowscale( highResImg, 5 );  titlenice( 'Proton Image' );

    if exist( 'objectiveValues', 'var' )
      figure; plotnice( objectiveValues ); titlenice( 'FISTA Objective Values' );
    end
  end
end


