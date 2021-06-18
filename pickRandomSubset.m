

function [ subset, complement ] = pickRandomSubset( set, nSubset )
  % [ subset, complement ] = pickRandomSubset( set, nSubset )
  %
  % Picks a random subset from a set (and possibly identifies its complement)
  %
  % Inputs:
  % set - an array that is the set to select from
  % nSubset - the number of elements in the subset
  %
  % Outputs:
  % subset - an array that is the randomly selected subset
  % complement - an array that is set - subset
  %
  % Written by Nicholas Dwork, Copyright 2021
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 2
    subset = [];  complement = [];
    disp( 'Usage:  [ subset, complement ] = pickRandomSubset( set, nSubset )' );
    return
  end

  nSet = numel( set );
  subsetIndxs = randsample( 1:nSet, nSubset );

  subset = set( subsetIndxs );

  if nargout > 1
    complement = setdiff( set, subset );
  end
end