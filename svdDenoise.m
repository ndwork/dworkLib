
function out = svdDenoise( data, newRank )
  % out = svdDenoise( data, newRank )
  %
  % Implements the denoising method of "Application of low-rank approximation using 
  % truncated singular value decomposition for noise reduction in hyperpolarized 13C
  % NMR spectroscopy" by Francischello et al.
  %
  % Inputs:
  % data - an array of dimensions greater than 1
  % newRank - the new Rank of the denoised data.  1 <= newRank <= size( data, last )
  %
  % Written by Nicholas Dwork, Copyright 2020
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 1
    disp( 'Usage:  out = svdDenoise( data, newRank )' );
    return;
  end

  sData = size( data );
  if newRank > sData( end )
    error( 'newRank must be less than the size of the last dimension of data.' );
  end

  data = reshape( data, [ prod( sData(1:end-1) ) sData(end) ] );
  [u,s,v] = svd( data, 'econ' );
  s = diag( s );

  if newRank < numel(s)
    s( newRank + 1 : end ) = 0;
  end
  out = u * diag(s )* v';

  out = reshape( out, sData );
end
