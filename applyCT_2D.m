
function out = applyCT_2D( f, domTraj, rangeTraj, kCy, kCx, Cy, Cx )
  % out = applyCT_2D( f, domTraj, N, kCy, kCx, Cy, Cx )
  % or
  % out = applyCT_2D( f, N, rangeTraj, kCy, kCx, Cy, Cx )
  % or
  % out = applyCT_2D( f, domTraj, rangeTraj, kCy, kCx, Cy, Cx )
  %
  % Applies a continuous circular convolution of a kernel as detailed in
  % http://nicholasdwork.com/tutorials/dworkGridding.pdf
  %
  % Inputs:
  %   f - An nTraj array representing the values of the function evaluated at each point in domTraj
  %   domTraj - An nTraj x 2 array specifying the ky / kx (first / second column) coordinates
  %             of the domain trajectory points
  %             OR
  %             a two element array specifying the size of the grid in the Fourier domain
  %   rangeTraj - An nNew x 2 array specifying the ky / kx (first / second column) coordinates
  %             of the new points
  %             OR
  %             a two element array specifying the size of the grid in the Fourier domain
  %   kCy - array of convolution kernel domain values in y dimension
  %   kCx - array of convolution kernel domain values in x dimension
  %   Cy - array of convolution kernel values in y dimension
  %   Cx - array of convolution kernel values in x dimension
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  % Adjoint is just applyC with the trajectories reversed
  out = applyC_2D( f, rangeTraj, domTraj, kCy, kCx, Cy, Cx );
end

