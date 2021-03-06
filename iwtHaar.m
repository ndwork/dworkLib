
function sig = iwtHaar( wt, varargin )
  % sig = iwtHaar( wt [, split] );
  % Performs an inverse Haar wavelet transform of a signal
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
  % Written by Nicholas Dwork - Copyright 2016
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
    s1 = split(1:nSplit/2);
    s2 = split(nSplit/2+1:end);

    if sum(s1)>0
      wt1 = iwtHaar( wt1, s1 );
    end
    if sum(s2)>0
      wt2 = iwtHaar( wt2, s2 );
    end
  end

  sig1 = upsample( wt1, 2 );
  sig1(2:2:end) = wt1;

  sig2 = upsample( wt2, 2 );
  sig2 = [ sig2(1); sig2(2:end) - sig2(1:end-1); ];

  sig = 1/sqrt(2) * ( sig1 + sig2 );
end
