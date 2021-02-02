
function sig = iwtDaubechies( wt, varargin )
  % sig = iwtDaubechies( wt [, split] );
  % Performs an inverse Daubechies wavelet transform of a signal with 
  %   circular boundary conditions
  % Note that this is a unitary transform
  %
  % Inputs:
  % wt - 1D array representing the wavelet transform of the signal
  %
  % Optional Inputs:
  % split - array specifying the number of levels of the wavelet transform.
  %   by default, split is 1 (indicating only one level).
  %   Example: [1 1 1 0] will have 3 levels.  The size of the last portion
  %   in the final level will be double that of the other portions since it
  %   wasn't split.
  %
  % Written by Nicholas Dwork - Copyright 2019
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  defaultSplit = 1;
  p = inputParser;
  p.addOptional( 'split', defaultSplit );
  p.parse( varargin{:} );
  split = p.Results.split;

  nwt = numel(wt);
  wt1 = wt(1:nwt/2);
  wt2 = wt(nwt/2+1:end);

  nSplit = numel(split);
  if nSplit > 1
    split1 = split(1:nSplit/2);
    split2 = split(nSplit/2+1:end);

    if sum(split1)>0
      wt1 = iwtDaubechies( wt1, split1 );
    end
    if sum(split2)>0
      wt2 = iwtDaubechies( wt2, split2 );
    end
  end

  tmp = upsample( wt1, 2 );
  tmp3 = 3 * tmp;
  tmpSqrt3 = tmp * sqrt(3);
  sig1 = tmp - tmpSqrt3 + ...
    circshift( tmp3 - tmpSqrt3, 1 ) + ...
    circshift( tmp3 + tmpSqrt3, 2 ) + ...
    circshift( tmp + tmpSqrt3, 3 );
  
  tmp = upsample( wt2, 2 );
  tmp3 = 3 * tmp;
  tmpSqrt3 = tmp * sqrt(3);
  sig2 = -( tmp + tmpSqrt3 ) + ...
    circshift( tmp3 + tmpSqrt3, 1 ) + ...
    circshift( -( tmp3 - tmpSqrt3 ), 2 ) + ...
    circshift( tmp - tmpSqrt3, 3 );
  
  sig = ( sig1 + sig2 ) / ( 4 * sqrt(2) );
end
