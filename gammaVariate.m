function out = gammaVariate( alpha1, beta1, A0, t0, dValues )
  % Standard parameterization of the gamma variate function
  %
  % Inputs:
  % dValues - domain values specifying where to evaluate the function
  %
  % Ouputs:
  % out - values of gamma variate function at sample times
  %
  % Written by Nicholas Dwork, 2019
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  dts = dValues - t0;

  gvs = A0 * dts.^alpha1 .* exp( -dts / beta1 );

  gvs( dts <= 0 ) = 0;
  out = gvs(:);
end

