
function mdm = calcMetricMDM( img, varargin )
  % mdm = calcMetricMDM( img [, 'p', p, 'q' q ] )
  %
  % Calculates the MDM metric based on the Minkowski distance as described
  % in "Efficient No-Reference Quality Assessment and Classification Model
  % for Contrast Distorted Images" by Nafchi and Cheriet, 2018
  %
  % Inputs:
  % img - array
  %
  % Optional Inputs:
  % p - Mdm parameter (default is 64)
  % q - Mdm parameter (default is 8)
  %
  % Written by Nicholas Dwork, Copyright 2019
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  ip = inputParser;
  ip.addParameter( 'p', 64, @ispositive );
  ip.addParameter( 'q', 8, @ispositive );
  ip.parse( varargin{:} );
  p = ip.Results.p;
  q = ip.Results.q;

  imgq = img .^ q;
  meanq = mean( imgq(:) );
  nImg = numel( img );
  
  mdm = ( sum( abs( imgq(:) - meanq ).^p ) / nImg ).^(1/(4*p));
end
