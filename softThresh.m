
function out = softThresh( in, thresh )
  % out = softThresh( in, thresh )
  %
  % Applies the proximal operator of the L1 norm (proxL1)
  %
  % Note that this functionality for real inputs is implemented by Matlab
  % with the wthresh function, but it is part of the Wavelet Toolbox.  If
  % you don't have that toolbox, then this function will work fine.  :)
  %
  % Inputs:
  % in - scalar or array
  % thresh - the soft thresholding value
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 2
    disp( 'Usage: out = softThresh( in, thresh )' );
    return
  end

  if thresh == 0

    out = in;

  elseif isreal( in )

    out = sign(in) .* max( ( abs(in) - thresh ), 0 );

  else

    magIn = abs( in );
    magOut = max( ( magIn - thresh ), 0 );
    out = in ./ magIn .* magOut;
    out( magIn == 0 ) = 0;

  end

end
