function out = structuralBilateralFilter( img, varargin )
%out = structuralBilateralFilter( img, [ 'kSize', kSize, ...
  %  'searchSize', searchSize, 'sigmaS', sigmaS, 'sigmaD', sigmaD ] );
  %
  % Inputs
  %   img: the input image to be denoised (a 2D array)
  %   kSize: (optional) the length of each side of the kernel
  %   searchSize: (optional) the length of each size of the square
  %     neighborhood
  %   sigmaS: (optional) structural spread parameter
  %   sigmaD: (optional) euclidean distance spread parameter
  %
  % Outputs
  %   out: the denoised image
  
  if ismatrix(img)
    out = structuralBilateralFilter_2D( img, varargin{:} );
  elseif ndims(img)==3
    out = structuralBilateralFilter_3D( img, varargin{:} );
  end
end


function out = structuralBilateralFilter_3D( img, varargin )
  %out = structuralBilateralFilter( img, [ 'kSize', kSize, ...
  %  'searchSize', searchSize,'sigmaS', sigmaS, 'sigmaD', sigmaD ] );
  %
  % Inputs
  %   img: the input image to be denoised (a 3D array)
  %   kSize: (optional) the length of each side of the kernel
  %   searchSize: (optional) the length of each size of the square
  %     neighborhood
  %   sigmaS: (optional) structural spread parameter
  %   sigmaD: (optional) euclidean distance spread parameter
  %
  % Outputs
  %   out: the denoised image

  isIntAndPositive = @(x) isInteger(x) && (x>0);
  isIntAndPositiveAndOdd = @(x) isIntAndPositive(x) && (mod(x,2)~=0);
  isNumericAndPositive = @(x) isnumeric(x) && (x>0);

  defaultKSize = 7;
  defaultSSize = 31;
  defaultSigmaS = 1.0;   % assumes images scaled from 0 to 1
  defaultSigmaD = 3;
  p = inputParser;
  p.addRequired( 'img', @isnumeric );
  p.addParamValue( 'kSize', defaultKSize, isIntAndPositiveAndOdd );
  p.addParamValue( 'searchSize', defaultSSize, isIntAndPositiveAndOdd );
  p.addParamValue( 'sigmaS', defaultSigmaS, isNumericAndPositive );
  p.addParamValue( 'sigmaD', defaultSigmaD, isNumericAndPositive );
  p.parse( img, varargin{:} );
  kSize = p.Results.kSize;
  sSize = p.Results.searchSize;
  sigmaS = p.Results.sigmaS;
  sigmaD = p.Results.sigmaD;


  dKernel = fspecial( 'gaussian', sSize, sigmaD );
  halfSearchSize = floor( sSize/2 );
  halfKSize = floor( kSize/2 );
  sVar = sigmaS*sigmaS;

  [M, N] = size( img );
  borderSize = halfKSize + halfSearchSize;

  out = zeros( size(img) );
  
  parfor j=borderSize+1:M-borderSize
    if mod(j,5)==0
      disp(['SBF: Working on ', num2str(j), ' of ', num2str(M)]);
    end

    localWeights = zeros( sSize, sSize );
    for i=borderSize+1:N-borderSize

      kernel = img( j-halfKSize : j+halfKSize, ...
                    i-halfKSize : i+halfKSize );

      for jP=0:sSize-1
        vJ = j-halfSearchSize+jP;

        for iP=0:sSize-1          
          vI = i-halfSearchSize+iP;

          v = img( vJ-halfKSize : vJ+halfKSize, ...
                   vI-halfKSize : vI+halfKSize  );

          distSq = ( kernel - v ) .* ( kernel - v );
          distSq = sum( distSq(:) ); %L2 norm squared

          localWeights( jP+1, iP+1 ) = exp( - distSq / sVar );
        end
      end

      localWeights = localWeights .* dKernel;
      localWeights = localWeights / sum( localWeights(:) );

      subImg = img( j-halfSearchSize : j+halfSearchSize, ...
                    i-halfSearchSize : i+halfSearchSize  );

      tmp = sum( localWeights(:) .* subImg(:) );
      out(j,i) = tmp;

    end
  end

end


function out = structuralBilateralFilter_2D( img, varargin )
  %out = structuralBilateralFilter( img, [ 'kSize', kSize, ...
  %  'searchSize', searchSize, 'sigmaS', sigmaS, 'sigmaD', sigmaD ] );
  %
  % Inputs
  %   img: the input image to be denoised (a 2D array)
  %   kSize: (optional) the length of each side of the kernel
  %   searchSize: (optional) the length of each size of the square
  %     neighborhood
  %   sigmaS: (optional) structural spread parameter
  %   sigmaD: (optional) euclidean distance spread parameter
  %
  % Outputs
  %   out: the denoised image

  isIntAndPositive = @(x) isInteger(x) && (x>0);
  isIntAndPositiveAndOdd = @(x) isIntAndPositive(x) && (mod(x,2)~=0);
  isNumericAndPositive = @(x) isnumeric(x) && (x>0);

  defaultKSize = 7;
  defaultSSize = 31;
  defaultSigmaS = 1.0;   % assumes images scaled from 0 to 1
  defaultSigmaD = 3;
  p = inputParser;
  p.addRequired( 'img', @isnumeric );
  p.addParamValue( 'kSize', defaultKSize, isIntAndPositiveAndOdd );
  p.addParamValue( 'searchSize', defaultSSize, isIntAndPositiveAndOdd );
  p.addParamValue( 'sigmaS', defaultSigmaS, isNumericAndPositive );
  p.addParamValue( 'sigmaD', defaultSigmaD, isNumericAndPositive );
  p.parse( img, varargin{:} );
  kSize = p.Results.kSize;
  sSize = p.Results.searchSize;
  sigmaS = p.Results.sigmaS;
  sigmaD = p.Results.sigmaD;


  dKernel = fspecial( 'gaussian', sSize, sigmaD );
  halfSearchSize = floor( sSize/2 );
  halfKSize = floor( kSize/2 );
  sVar = sigmaS*sigmaS;

  [M, N] = size( img );
  borderSize = halfKSize + halfSearchSize;

  out = zeros( size(img) );
  
  parfor j=borderSize+1:M-borderSize
    if mod(j,5)==0
      disp(['SBF: Working on ', num2str(j), ' of ', num2str(M)]);
    end

    localWeights = zeros( sSize, sSize );
    for i=borderSize+1:N-borderSize

      kernel = img( j-halfKSize : j+halfKSize, ...
                    i-halfKSize : i+halfKSize );

      for jP=0:sSize-1
        vJ = j-halfSearchSize+jP;

        for iP=0:sSize-1          
          vI = i-halfSearchSize+iP;

          v = img( vJ-halfKSize : vJ+halfKSize, ...
                   vI-halfKSize : vI+halfKSize  );

          distSq = ( kernel - v ) .* ( kernel - v );
          distSq = sum( distSq(:) ); %L2 norm squared

          localWeights( jP+1, iP+1 ) = exp( - distSq / sVar );
        end
      end

      localWeights = localWeights .* dKernel;
      localWeights = localWeights / sum( localWeights(:) );

      subImg = img( j-halfSearchSize : j+halfSearchSize, ...
                    i-halfSearchSize : i+halfSearchSize  );

      tmp = sum( localWeights(:) .* subImg(:) );
      out(j,i) = tmp;

    end
  end

end
