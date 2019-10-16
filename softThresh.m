
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

  if isreal( in )

    out = sign(in) .* max( ( abs(in) - thresh ), 0 );
    
  else

    magIn = abs( in );

    scalingFactors = thresh ./ magIn;
    scalingFactors( magIn <= thresh ) = 1;

    projsOntoUnitMag = in .* scalingFactors;

    out = in - projsOntoUnitMag;

  end

end
