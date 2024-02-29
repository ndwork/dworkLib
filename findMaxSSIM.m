
function [out,k] = findMaxSSIM( recon, trueRecon, varargin )
  % Result is max_k ssim( k * recon - trueRecon )
  %
  % out = findMaxSSIM( recon, trueRecon [, 'tol', tol, 'verbose', true/false ] )
  %
  % Written by Nicholas Dwork, Copyright 2024
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 2
    disp( 'Usage:  out = findMaxSSIM( recon, trueRecon [, ''verbose'', true/false ] )' );
    if nargout > 0, out = []; end
    return;
  end

  p = inputParser;
  p.addParameter( 'tol', [] );
  p.addParameter( 'verbose', true );
  p.parse( varargin{:} );
  tol = p.Results.tol;
  verbose = p.Results.verbose;

  f = @(k) -ssim( k*recon, trueRecon );

  LB = 0;
  UB = max( recon(:) ) / mean( trueRecon(:) ) * 10;

  k = goldenSectionSearch( f, LB, UB, 'tol', 1d-4, 'verbose', verbose );
  out = f( k );
end
