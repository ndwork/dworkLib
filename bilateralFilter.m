
function out = bilateralFilter( img, varargin )
  % out = bilateralFilter( img, [ 'S', S, 'sigmaD', sigmaD, 'sigmaR', sigmaR ] );
  %
  % Inputs
  %   img: the input image to be denoised
  %     either a MxN 2D grayscale image or a MxNx3 color image
  %   S: (optional) the length of each side of the kernel (must be odd)
  %   sigmaD: (optional) euclidean distance spread parameter
  %   sigmaR: (optional) photometric distance spread parameter
  %
  % Outputs
  %   out: the denoised image

  if ismatrix(img)
    out = bilateralFilter_2D( img, varargin{:} );
  elseif ndims(img)==3
    out = bilateralFilter_3D( img, varargin{:} );
  end
end


function out = bilateralFilter_3D( img, varargin )
  % out = bilateralFilter( img, [ 'S', S, 'sigmaD', sigmaD, 'sigmaR', sigmaR ] );
  %
  % Inputs
  %   img: the input image to be denoised (a 2D array)
  %   S: (optional) the length of each side of the kernel (must be odd)
  %   sigmaD: (optional) euclidean distance spread parameter
  %   sigmaR: (optional) photometric distance spread parameter
  %
  % Outputs
  %   out: the denoised image


  isIntAndPositive = @(x) isInteger(x) && (x>0);
  isIntAndPositiveAndOdd = @(x) isIntAndPositive(x) && (mod(x,2)~=0);
  isNumericAndPositive = @(x) isnumeric(x) && (x>0);

  defaultS = 31;
  defaultSigmaD = 3;
  defaultSigmaR = 0.3;
  p = inputParser;
  p.addRequired( 'img', @isnumeric );
  p.addParamValue( 'S', defaultS, isIntAndPositiveAndOdd );
  p.addParamValue( 'sigmaD', defaultSigmaD, isNumericAndPositive );
  p.addParamValue( 'sigmaR', defaultSigmaR, isNumericAndPositive );
  p.parse( img, varargin{:} );
  s = p.Results.S;
  sigmaD = p.Results.sigmaD;
  sigmaR = p.Results.sigmaR;

  dKernel = fspecial3d_gaussian( 'gaussian', s, sigmaD );
  halfS = floor(s/2);
  varR = sigmaR * sigmaR;

  sImg = size( img );
  out = zeros( size(img) );
  for j=ceil(s/2):sImg(1)-floor(s/2)
    for i=ceil(s/2):sImg(2)-floor(s/2)
      for k=ceil(s/2):sImg(3)-floor(s/2)

        pKernel = exp( -( img(j,i,k) - subImg ).^2 / ( 2 * varR ) );

        weights = pKernel .* dKernel;
        weights = weights / sum( weights(:) );

        subImg = img( j-halfS:j+halfS, i-halfS:i+halfS, k-halfS:k+halfS );
        out(j,i,k) = sum( subImg(:) .* weights(:) );
      end
    end
  end

end


function out = bilateralFilter_2D( img, varargin )
  % out = bilateralFilter( img, [ 'S', S, 'sigmaD', sigmaD, 'sigmaR', sigmaR ] );
  %
  % Inputs
  %   img: the input image to be denoised (a 2D array)
  %   S: (optional) the length of each side of the kernel (must be odd)
  %   sigmaD: (optional) euclidean distance spread parameter
  %   sigmaR: (optional) photometric distance spread parameter
  %
  % Outputs
  %   out: the denoised image


  isIntAndPositive = @(x) isInteger(x) && (x>0);
  isIntAndPositiveAndOdd = @(x) isIntAndPositive(x) && (mod(x,2)~=0);
  isNumericAndPositive = @(x) isnumeric(x) && (x>0);

  defaultS = 31;
  defaultSigmaD = 3;
  defaultSigmaR = 0.3;
  p = inputParser;
  p.addRequired( 'img', @isnumeric );
  p.addParamValue( 'S', defaultS, isIntAndPositiveAndOdd );
  p.addParamValue( 'sigmaD', defaultSigmaD, isNumericAndPositive );
  p.addParamValue( 'sigmaR', defaultSigmaR, isNumericAndPositive );
  p.parse( img, varargin{:} );
  s = p.Results.S;
  sigmaD = p.Results.sigmaD;
  sigmaR = p.Results.sigmaR;

  dKernel = fspecial( 'gaussian', s, sigmaD );
  halfS = floor(s/2);
  varR = sigmaR * sigmaR;

  sImg = size( img );
  out = zeros( size(img) );
  for j=ceil(s/2):sImg(1)-floor(s/2)
    for i=ceil(s/2):sImg(2)-floor(s/2)

      subImg = img( j-halfS:j+halfS, i-halfS:i+halfS );

      pKernel = exp( -( img(j,i) - subImg ).^2 / ( 2 * varR ) );

      weights = pKernel .* dKernel;
      weights = weights / sum( weights(:) );

      out(j,i) = sum( subImg(:) .* weights(:) );
    end
  end

end
