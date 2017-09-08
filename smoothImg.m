
function smoothed = smoothImg( img, varargin )
  % out = smoothImg( img [, N ] );
  % out = smoothImg( img [, N, 'gaussian', sigma ] );

  defaultN = 5;
  defaultSigma = 0;
  p = inputParser;
  p.addOptional( 'N', defaultN );
  p.addParameter( 'gaussian', defaultSigma, @isnumeric );
  p.parse( varargin{:} );
  N = p.Results.N;
  sigma = p.Results.gaussian;

  if sigma > 0
    %smoothed = imgaussfilt( img, sigma );
    h = fspecial( 'gaussian', N, sigma );
    smoothed = imfilter( img, h );
  else
    h = fspecial( 'average', N );
    smoothed = imfilter( img, h );
  end
end
