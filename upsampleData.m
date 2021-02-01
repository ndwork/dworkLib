
function out = upsampleData( data, U, varargin )
  % out = upsampleData( data, U [, 'S', S, 'sOut', sOut ] )
  %
  % Inputs:
  % data - array
  % U - upsample factor.  Either an scalar integer (assuming same factor in all dimensions)
  %     or an array of integers (one for each dimension)
  %
  % Optional Inputs:
  % S - amount to shift input, either a scalar integer or an array
  % sOut - the size of the output image.  By default, equals size(data) .* U;
  %
  % Outputs:
  % out - the upsampled image of size sOut
  %
  % Written by Nicholas Dwork - Copyright 2021
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addRequired( 'U', @isnumeric );
  p.addParameter( 'S', zeros( ndims(data), 1 ), @isnumeric );
  p.addParameter( 'sOut', [] );
  p.parse( U, varargin{:} );
  S = p.Results.S;
  sOut = p.Results.sOut;

  sData = size( data );
  nDims = ndims( data );

  if numel(U) == 1, U = U * ones( nDims, 1 ); end
  if numel(S) == 1, S = S * ones( nDims, 1 ); end
  if numel( sOut ) == 0, sOut = sData(:) .* U(:); end

  out = zeros( sOut(:)' );
  cmd = 'out( ';
  for dimIndx = 1 : nDims
    theseIndxs = [ 'S(', num2str(dimIndx), ') + (0:sData(', num2str(dimIndx), ...
      ')-1) * U(', num2str(dimIndx), ') + 1' ];
    cmd = [ cmd, theseIndxs ];  %#ok<AGROW>
    if dimIndx < nDims
      cmd = [ cmd, ', ' ];  %#ok<AGROW>
    else
      cmd = [ cmd, ' ) = data;' ];  %#ok<AGROW>
    end
  end

  eval( cmd );
end
