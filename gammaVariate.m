function out = gammaVariate( dValues, alpha, beta, A0, t0, varargin )
  % out = gammaVariate( dValues, alpha, beta, A0, t0 [, deriv ] )
  %
  % Standard parameterization of the gamma variate function
  %
  % Inputs:
  % dValues - domain values specifying where to evaluate the function
  %
  % Optional Inputs:
  % deriv - if true, then evaluate the derivative
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

  p = inputParser;
  p.addOptional( 'deriv', 0 );
  p.parse( varargin{:} );
  deriv = p.Results.deriv;

  dts = dValues - t0;

  if deriv == 0
    out = A0 * dts.^alpha .* exp( -dts / beta );
  else
    % calculate the first derivative
    out = A0 * dts.^(alpha-1) .* exp( -dts / beta ) .* ( alpha - dts / beta );
  end

  out( dts <= 0 ) = 0;
  out = out(:);
end

