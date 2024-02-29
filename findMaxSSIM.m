
function [out,k] = findMaxSSIM( est, truth, varargin )
  % Result is max_k ssim( k * est - truth )
  %
  % out = findMaxSSIM( est, truth [, 'tol', tol, 'verbose', true/false ] )
  %
  % Inputs:
  % est - the estimate array
  % truth - the truth array
  %
  % Optional Inputs:
  % tol - the tolerance of k
  % verbose - whether or not to display verbosity messages
  %
  % Written by Nicholas Dwork, Copyright 2024
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 2
    disp( 'Usage:  out = findMaxSSIM( est, truth [, ''tol'', tol, ''verbose'', true/false ] )' );
    if nargout > 0, out = []; end
    return;
  end

  p = inputParser;
  p.addParameter( 'tol', [] );
  p.addParameter( 'verbose', false );
  p.parse( varargin{:} );
  tol = p.Results.tol;
  verbose = p.Results.verbose;

  f = @(k) -ssim( k * est, truth );

  LB = 0;
  UB = max( est(:) ) / mean( truth(:) ) * 10;

  k = goldenSectionSearch( f, LB, UB, 'tol', tol, 'verbose', verbose );
  out = f( k );
end
