function out = nonlocalMeans( img, varargin )
  %out = nonlocalMeans( img, [ 'kSize', kSize, 'searchSize', searchSize, ...
  %  'sigmaS', sigmaS ] );
  %
  % Inputs
  %   img: the input image to be denoised (a 2D array)
  %   kSize: (optional) the length of each side of the kernel
  %   searchSize: (optional) the length of each size of the square
  %     neighborhood
  %   sigmaS: (optional) structural spread parameter
  %
  % Outputs
  %   out: the denoised image
  
  if ismatrix(img)
    out = nonlocalMeans_2D( img, varargin{:} );
  elseif ndims(img)==3
    out = nonlocalMeans_3D( img, varargin{:} );
  end
end

function out = nonlocalMeans_3D( img, varargin )
  %out = nonlocalMeans( img, [ 'kSize', kSize, 'searchSize', searchSize, ...
  %  'sigmaS', sigmaS ] );
  %
  % Inputs
  %   img: the input image to be denoised (a 3D array)
  %   kSize: (optional) the length of each side of the kernel
  %   searchSize: (optional) the length of each size of the square
  %     neighborhood
  %   sigmaS: (optional) structural spread parameter
  %
  % Outputs
  %   out: the denoised image

  isIntAndPositive = @(x) isInteger(x) && (x>0);
  isIntAndPositiveAndOdd = @(x) isIntAndPositive(x) && (mod(x,2)~=0);
  isNumericAndPositive = @(x) isnumeric(x) && (x>0);

  defaultKSize = 7;
  defaultSSize = 31;
  defaultSigmaS = 1.0;   % assumes images scaled from 0 to 1
  p = inputParser;
  p.addRequired( 'img', @isnumeric );
  p.addParameter( 'kSize', defaultKSize, isIntAndPositiveAndOdd );
  p.addParameter( 'searchSize', defaultSSize, isIntAndPositiveAndOdd );
  p.addParameter( 'sigmaS', defaultSigmaS, isNumericAndPositive );
  p.parse( img, varargin{:} );
  kSize = p.Results.kSize;
  sSize = p.Results.searchSize;
  sigmaS = p.Results.sigmaS;


  halfSearchSize = floor( sSize/2 );
  halfKSize = floor( kSize/2 );
  sVar = sigmaS*sigmaS;

  [M, N, K] = size( img );
  borderSize = halfKSize + halfSearchSize;

  out = zeros( size(img) );
  
  parfor j=borderSize+1:M-borderSize
    if mod(j,5)==0
      disp(['NLM: Working on ', num2str(j), ' of ', num2str(M)]);
    end

    localWeights = zeros( sSize, sSize );
    for i=borderSize+1:N-borderSize
      
      for k=borderSize+1:K-borderSize

        kernel = img( j-halfKSize : j+halfKSize, ...
                      i-halfKSize : i+halfKSize, ...
                      k-halfKSize : k+halfKSize  );

        for jP=0:sSize-1
          vJ = j-halfSearchSize+jP;

          for iP=0:sSize-1          
            vI = i-halfSearchSize+iP;

            for kP=0:sSize-1
              vK = k-halfSearchSize+kP;

              v = img( vJ-halfKSize : vJ+halfKSize, ...
                       vI-halfKSize : vI+halfKSize, ...
                       vK-halfKSize : vK+halfKSize  );

              distSq = ( kernel - v ) .* ( kernel - v );
              distSq = sum( distSq(:) ); %L2 norm squared

              localWeights( jP+1, iP+1, kP+1 ) = exp( - distSq / sVar );
            end
          end
        end

        localWeights = localWeights / sum( localWeights(:) );

        subImg = img( j-halfSearchSize : j+halfSearchSize, ...
                      i-halfSearchSize : i+halfSearchSize, ...
                      k-halfSearchSize : k+halfSearchSize  );

        tmp = sum( localWeights(:) .* subImg(:) );
        out(j,i,k) = tmp;

      end
    end
  end
end


function out = nonlocalMeans_2D( img, varargin )
  %out = nonlocalMeans( img, [ 'kSize', kSize, 'searchSize', searchSize, ...
  %  'sigmaS', sigmaS ] );
  %
  % Inputs
  %   img: the input image to be denoised (a 2D array)
  %   kSize: (optional) the length of each side of the kernel
  %   searchSize: (optional) the length of each size of the square
  %     neighborhood
  %   sigmaS: (optional) structural spread parameter
  %
  % Outputs
  %   out: the denoised image

  isPositive = @(x) x>0;
  isIntAndPositiveAndOdd = @(x) isPositive(x) && (mod(x,2)~=0);
  isNumericAndPositive = @(x) isnumeric(x) && (x>0);

  defaultKSize = 7;
  defaultSSize = 31;
  defaultSigmaS = 1.0;   % assumes images scaled from 0 to 1
  p = inputParser;
  p.addRequired( 'img', @isnumeric );
  p.addParameter( 'kSize', defaultKSize, isIntAndPositiveAndOdd );
  p.addParameter( 'searchSize', defaultSSize, isIntAndPositiveAndOdd );
  p.addParameter( 'sigmaS', defaultSigmaS, isNumericAndPositive );
  p.parse( img, varargin{:} );
  kSize = p.Results.kSize;
  sSize = p.Results.searchSize;
  sigmaS = p.Results.sigmaS;


  halfSearchSize = floor( sSize/2 );
  halfKSize = floor( kSize/2 );
  sVar = sigmaS*sigmaS;

  [M, N] = size( img );
  borderSize = halfKSize + halfSearchSize;

  out = zeros( size(img) );
  
  parfor j=borderSize+1:M-borderSize
    if mod(j,5)==0
      disp(['NLM: Working on ', num2str(j), ' of ', num2str(M)]);
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

      localWeights = localWeights / sum( localWeights(:) );

      subImg = img( j-halfSearchSize : j+halfSearchSize, ...
                    i-halfSearchSize : i+halfSearchSize  );

      tmp = sum( localWeights(:) .* subImg(:) );
      out(j,i) = tmp;

    end
  end

end
