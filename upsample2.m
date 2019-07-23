
function out = upsample2( img, n )
  % out = upsample2( img, n )
  % Peforms two dimensional upsampling
  %
  % Inputs:
  % n - if n is a scalar, upsamples by n in both dimensions
  % n - if n is a two element array, upsamples by n(i) in the ith dimension
  %
  % Outputs:
  % out - upsampled image
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if numel(n)==1, n = n * ones(ndims(img),1); end

  tmp = upsample( img, n(1) );
  out = upsample( tmp', n(2) )';
end
