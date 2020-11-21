
function out = bilateralFilter( img, varargin )
  % out = bilateralFilter( img, [ 'S', S, 'sigmaD', sigmaD, ...
  %   'sigmaR', sigmaR, 'verbose', verbose ] );
  %
  % Inputs
  %   img: the input image to be denoised
  %     either a MxN 2D grayscale image or a MxNx3 color image
  %
  % Optional Inputs:
  %   S: the length of each side of the kernel (must be odd)
  %   sigmaD: euclidean distance spread parameter
  %   sigmaR: photometric distance spread parameter
  %   verbose: displays status
  %
  % Outputs
  %   out: the denoised image

  if nargin < 1
    disp( 'Usage:  out = bilateralFilter( img, [ ''S'', S, ''sigmaD'', sigmaD, ... ' );
    disp( '  ''sigmaR'', sigmaR, ''verbose'', verbose ] );' );
    return
  end
  
  isPositiveAndOdd = @(x) (x>0) && (mod(x,2) == 1);
  isNumericAndPositive = @(x) isnumeric(x) && (x>0);

  defaultSigmaD = 3;
  defaultSigmaR = 0.3;
  p = inputParser;
  p.addRequired( 'img', @isnumeric );
  p.addParameter( 'S', [], isPositiveAndOdd );
  p.addParameter( 'sigmaD', defaultSigmaD, isNumericAndPositive );
  p.addParameter( 'sigmaR', defaultSigmaR, isNumericAndPositive );
  p.addParameter( 'verbose', false, @(x) islogical(x) || isnumeric(x) );
  p.parse( img, varargin{:} );
  s = p.Results.S;
  sigmaD = p.Results.sigmaD;
  sigmaR = p.Results.sigmaR;
  verbose = p.Results.verbose;
  
  if numel( s ) == 0, s = 3 * sigmaD; end

  if ismatrix(img)
    out = bilateralFilter_2D( img, s, sigmaD, sigmaR, verbose );
  elseif ndims(img)==3
    out = bilateralFilter_3D( img, s, sigmaD, sigmaR, verbose );
  end
end


function out = bilateralFilter_2D( img, s, sigmaD, sigmaR, verbose )

  dKernel = fspecial( 'gaussian', s, sigmaD );
  halfS = floor(s/2);
  varR = sigmaR * sigmaR;

  sImg = size( img );

  firstJ = ceil(s/2);   lastJ = sImg(1)-floor(s/2);
  firstI = ceil(s/2);   lastI = sImg(2)-floor(s/2);
  jCells = cell( lastJ - firstJ + 1, 1 );

  p = parforProgress( lastJ-firstJ+1 );
  parfor j=1:lastJ-firstJ+1
    if verbose == true, p.progress( j, 50 ); end

    tmp = zeros( 1, lastI-firstI+1 );
    imgLine = img( j+firstJ-1-halfS:j+firstJ-1+halfS, : );   %#ok<PFBNS>

    for i=1:lastI-firstI+1
      %subImg = img( j-halfS:j+halfS, i-halfS:i+halfS );
      subImg = imgLine( :, i+firstI-1-halfS : i+firstI-1+halfS );

      pKernel = exp( -( img(j+firstJ-1,i+firstI-1) - subImg ).^2 / ( 2 * varR ) );

      weights = pKernel .* dKernel;
      weights = weights / sum( weights(:) );

      tmp(i) = sum( subImg(:) .* weights(:) );
    end
    jCells{j} = tmp;
  end
  p.clean;
  out = img;
  out(firstJ:lastJ,firstI:lastI) = cell2mat( jCells );
end


function out = bilateralFilter_3D( img, s, sigmaD, sigmaR, verbose )

  dKernel = fspecial3d( 'gaussian', s, sigmaD );
  halfS = floor(s/2);
  varR = sigmaR * sigmaR;

  sImg = size( img );
  firstJ = ceil(s/2);   lastJ = sImg(1)-floor(s/2);
  firstI = ceil(s/2);   lastI = sImg(2)-floor(s/2);
  firstK = ceil(s/2);   lastK = sImg(3)-floor(s/2);
  sliceCells = cell(1,1,lastK-firstK+1);

  p = parforProgress( lastK-firstK+1 );
  parfor k=1:lastK-firstK+1
    if verbose == true, p.progress( k ); end   %#ok<PFBNS>

    tmp = zeros( sImg(1:2) );   %#ok<PFBNS>
    kImg = img(:,:,k+firstK-1-halfS:k+firstK-1+halfS);   %#ok<PFBNS>

    for j=1:lastJ-firstJ+1
      for i=1:lastI-firstI+1

        subImg = kImg( j+firstJ-1-halfS:j+firstJ-1+halfS, ...
                       i+firstI-1-halfS:i+firstI-1+halfS, : );

        pKernel = exp( -( img(j,i,k) - subImg ).^2 / ( 2 * varR ) );

        weights = pKernel .* dKernel;
        weights = weights / sum( weights(:) );

        tmp(j,i) = sum( subImg(:) .* weights(:) );
      end
    end
    sliceCells{1,1,k} = tmp;
  end
  p.clean;
  out = img;
  out(firstJ:lastJ,firstI:lastI,firstK:lastK) = cell2mat( sliceCells );
end


