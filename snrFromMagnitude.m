
function snr = snrFromMagnitude( in, noiseCoords )
  % snr = snrFromMagnitude( in, noiseCoords )
  %
  % Calculates the SNR of the underlying complex data with Guassian additive noise
  % from the magnitude data according to "The Rician Distribution of Noisy MRI Data"
  % by Gudbjartsson and Patz
  %
  % Inputs:
  % in - array of magnitude data
  % noiseCoords - indices of noise region within input array with values
  %   [ minNoise_dim1 maxNoise_dim1 ... minNoise_dimN maxNoise_dimN ]
  %   Any dimension where noise coordinates are not supplied is assumed to be
  %   all noise.
  %
  % Outputs:
  % snr - scalar representing the signal to noise ratio
  %
  % Written by Nicholas Dwork - Copyright 2021
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  subIndxs = repmat( {' : '}, 1, ndims(in) );
  for dimIndx = 1 : numel( noiseCoords ) / 2
    subIndxs{ dimIndx } = [ ...
      num2str( noiseCoords( 2 * ( dimIndx - 1 ) + 1 ) ), ' : ', ...
      num2str( noiseCoords( 2 * ( dimIndx - 1 ) + 2 ) ) ];
  end
  noise = [];
  subIndxs = join( subIndxs, ', ' );
  str2eval = [ 'noise = in( ', subIndxs{1}, ' );' ];
  eval( str2eval );

  varM = var( noise(:) );  % variance of magnitude values
  varNoise = varM / ( 2 - 0.5 * pi );

  signalPower = abs( in(:).^2 - varNoise );  % Unbiased signal power
  signalPower = reshape( signalPower, size( in ) );

  snr = signalPower / varNoise;
end
