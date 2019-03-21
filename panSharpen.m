
function out = panSharpen( M, C, varargin )
  % out = panSharpen( M, C/R [, G, B, IR, 'chWs', chWs, 'lambda', lambda, 'sigma', sigma ] )
  %
  % The algorithm implemented is described in "Fusion of Multispectral and Panchromatic
  % Images Using a Restoration-Based Method" by Zhenhua Li et al.
  %
  % Inputs:
  % Either
  %   C - a 3D array representing the red, green, and blue channels.
  %   R - a 2D array representing the red channel
  %
  % Optional Inputs:
  % G,B - green and blue images (2D arrays supplied if R is supplied instead of C)
  % IR - optional ifrared data
  % chWs - the weights used to combine each channel into the monochrome image (1D array)
  % sigma - the standard deviation of the Gaussian blur kernel
  %
  % Written by Nicholas Dwork, Copyright 2019
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  ds = size( M, 1 ) / size( C, 1 );

  p = inputParser;
  p.addOptional( 'G', [], @isnumeric );
  p.addOptional( 'B', [], @isnumeric );
  p.addOptional( 'IR', [], @isnumeric );
  p.addParameter( 'chWs', [], @isnumeric );
  p.addParameter( 'lambda', 1, @(x) isnumeric(x) && x>0 && numel(x)==1 );
  p.addParameter( 'sigma', ds/2, @isnumeric );
  p.parse( varargin{:} );
  G = p.Results.G;
  B = p.Results.B;
  IR = p.Results.IR;
  chWs = p.Results.chWs;
  lambda = p.Results.lambda;
  sigma = p.Results.sigma;

  if ismatrix( C )
    C = cat( 3, R, G, B );
    if numel( IR ) > 0, C = cat( 3, C, IR ); end
  end
  sC = size( C );
  nChannels = size( C, 3 );

  if numel( chWs ) == 0
    if nChannels == 3
      chWs = [ 0.2989, 0.5870, 0.1140 ];  % luminance weights
    elseif nChannels ==4
      chWs = [ 0.23, 0.24, 0.11, 0.42 ];  % quickbird weights
    end
  end

  x0 = repmat( M, [1 1 nChannels] );
  kSize = ceil( max( 5*sigma, 3 ) );
  for kIndx=1:numel(kSize)
    if mod( kSize(kIndx), 2 ) == 0, kSize(kIndx) = kSize(kIndx) + 1; end
  end

  sM = size(M);
  sCh = [ size(C,1), size(C,2) ];

  lum = zeros( sCh );  % luminance image
  function out = mA( in, arg )
    if nargin < 2 || strcmp( arg, 'notransp' )
      in = reshape( in, size(x0) );
      lum = sum( volVectorProd( in, chWs, 3 ), 3 );
      out = lum(:);
    else
      % Return the adjoint
      in = reshape( in, sM );
      tmp = repmat( in, [1 1 nChannels] );
      out = volVectorProd( tmp, chWs, 3 );
      out = out(:);
    end
  end
  %[out,err] = checkAdjoint( x0(:), @mA );

  xs = ones(sM(1),1) * (1:sM(2));
  ys = (1:sM(1))' * ones(1,sM(2));
  xqs = imresize( xs, sCh, 'nearest' );
  yqs = imresize( ys, sCh, 'nearest' );
  function out = cA( in, arg )
    if nargin < 2 || strcmp( arg, 'notransp' )
      in = reshape( in, size(x0) );
      out = zeros( sC );
      for ch=1:nChannels
        smoothed = smoothImg( in(:,:,ch), kSize, 'gaussian', sigma );
        out( :, :, ch ) = imresize( smoothed, sCh, 'nearest' );
      end
      out = out(:);

    else
      % Return the adjoint
      in = reshape( in, sC );
      out = zeros( size(x0) );
      for ch = 1:nChannels
        tmp = zeros( sM );

        for i=1:size(C,2)
          for j=1:size(C,1)
            tmp( yqs(j,i), xqs(j,i) ) = in( j, i, ch );
          end
        end

        out(:,:,ch) = smoothImg( tmp, kSize, 'gaussian', sigma, 'op', 'transp' );
      end

      out = out(:);
    end
  end
  %[out,err] = checkAdjoint( x0(:), @cA );

  lambdaSq = lambda * lambda;
  function out = A( in, arg )
    if nargin < 2 || strcmp( arg, 'notransp' )
      mOut = mA( in );
      cOut = lambdaSq * cA( in );
      out = [ mOut; cOut; ];
    else
      % Return the adjoint of the A operator
      mIn = in(1:numel(M));
      cIn = in(numel(M)+1:end);
      out = mA( mIn, arg ) + lambdaSq * cA( cIn, arg );
    end
  end
  %[out,err] = checkAdjoint( x0(:), @A );

  b = [ M(:); lambdaSq * C(:); ];
  [x,flag] = lsqr( @A, b, [], [], [], [], x0(:) );    %#ok<ASGLU>
  out = reshape( x, [ size(M) nChannels ] );
end

