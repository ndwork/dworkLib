
function out = uifft( in, varargin )
  % out = uifft( in [, dim ] )
  %   or
  % out = uifft( in [, n, dim ] )
  %
  % Compute the unitary inverse fft
  %
  % Written by Nicholas Dwork, Copyright 2019
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addOptional( 'n', [], @(x) numel(x) == 0 || ispositive(x) );
  p.addOptional( 'dim', [], @isnumeric );
  p.parse( varargin{:} );
  dim = p.Results.dim;
  n = p.Results.n;

  if numel( dim ) == 0 && numel( n ) ~= 0
    dim = n;  n = [];
  end

  if numel( dim ) == 0
    tmp = squeeze( in );
    if size( tmp, 1 ) > 1
      dim = size( tmp, 1 );
    else
      dim = size( tmp, 2 );
    end
  end

  if numel( n ) == 0
    n = size( in, dim );
  end

  out =  ifft( in, n, dim ) * sqrt( n );
end
