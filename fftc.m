
function out = fftc( in, varargin )
  % out = fftc( in, varargin )
  %
  % Compute the cenetered fft
  %
  % Inputs:
  % in - an array (of any number of dimensions)
  %
  % Written by Nicholas Dwork, Copyright 2019
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParers;
  p.addOptional( 'type', [], @(x) true )
  p.parse( varargin{:} );
  type = p.Results.type;

  switch ndims( in )
    case 1
      out = fftshift( fft( in ) );
    case 2
      out = fftshift( fft2( in ) );
    otherwise
      out = fftshift( fftn( in ) );
  end

  if strcmp( type, 'unitary' )
    out = out / sqrt( numel(in) );
  end

end
