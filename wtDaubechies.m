
function wt = wtDaubechies( sig, varargin )
  % wt = wtDaubechies( sig [, split] );
  %
  % Performs a Daubechies wavelet transform of a signal with circular boundary conditions
  % Based on the Wikipedia page on the Daubechies wavelet transform and
  % (http://wavelets.pybytes.com/wavelet/db2/)
  % Note that this is a unitary transform
  %
  % Inputs:
  % sig - 1D array representing the signal to be transformed
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
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  defaultSplit = 1;
  p = inputParser;
  p.addOptional( 'split', defaultSplit );
  p.parse( varargin{:} );
  split = p.Results.split;

  sigSqrt3 = sig * sqrt(3);
  sig3 = 3 * sig;

  sigPsigSqrt3 = sig + sigSqrt3;
  sig3PsigSqrt3 = sig3 + sigSqrt3;
  sigMsigSqrt3 = sig - sigSqrt3;
  sig3MsigSqrt3 = sig3 - sigSqrt3;
  
  wt1 = ( sigMsigSqrt3 + ...
    circshift( sig3MsigSqrt3, -1 ) + ...
    circshift( sig3PsigSqrt3, -2 ) + ...
    circshift( sigPsigSqrt3, -3 ) );
  wt1 = wt1(1:2:end);

  wt2 = ( -sigPsigSqrt3 + ...
    circshift( sig3PsigSqrt3, -1 ) + ...
    circshift( -sig3MsigSqrt3, -2 ) + ...
    circshift( sigMsigSqrt3, -3 ) );
  wt2 = wt2(1:2:end);

  nSplit = numel(split);
  if nSplit > 1
    split1 = split(1:nSplit/2);
    split2 = split(nSplit/2+1:end);

    if sum(split1)>0
      wt1 = wtDaubechies( wt1, split1 );
    end    
    if sum(split2)>0
      wt2 = wtDaubechies( wt2, split2 );
    end
  end

  wt = [wt1; wt2;] / ( 4 * sqrt(2) );
end

