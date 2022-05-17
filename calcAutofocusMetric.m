
function out = calcAutofocusMetric( in )
  % Calculates the autofocus metric according to "Automatic compensation of motion
  % artifacts in MRI" by Loktyushin et al.
  %
  % Input:
  % in - an array of any number of dimensions
  %
  % Output:
  % out - the value of the autofocus metric
  %
  % Written by Nicholas Dwork, Copyright 2022
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular purpose.

  gradIn = computeGradient( in );

  nDimsIn = ndims( in );
  gradIn = reshape( gradIn, [ numel( in ) nDimsIn ] );

  out = 0;
  for dimIndx = 1 : nDimsIn
    v = sqrt( ...
      gradIn(:,dimIndx) .* conj( gradIn(:,dimIndx) ) ./ ...
      gradIn(:,dimIndx)' * gradIn(:,dimIndx) );

    out = out + -v' * log( v );
  end

end

