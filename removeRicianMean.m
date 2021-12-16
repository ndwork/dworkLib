
function out = removeRicianMean( data, noise )
  % out = removeRicianMean( data, noise )
  %
  % Removes the mean of the Rician noise according to "The Rician Distribution of
  % Noisy MRI Data" by Gudbjartsson and Patz.
  %
  % Inputs:
  % data - a real array of data with Rician noise
  % noise - either a complex array of the underlying Gaussian noise or a portion of the
  %         real data without any signal
  %
  % Outputs:
  % out - an array of the size of data with an estimate of the mean of the noise removed
  %
  % Written by Nicholas Dwork, Copyright 2021
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular purpose.

  if isreal( noise )
    gStd = mean( noise(:) ) * sqrt( 0.5 * pi );  % Gausian standard deviation
  else
    gStd = mean([ std( real(noise) ), std( imag(noise) ) ]);  % Gausian standard deviation
  end

  out = sqrt( abs( abs( data ).^2 - gStd^2 ) );

  lowIndxs = abs( data ) / gStd < 2;
  out( lowIndxs ) = data( lowIndxs ) - ( gStd / sqrt( 0.5 * pi ) );  % Subtract noise mean
end
