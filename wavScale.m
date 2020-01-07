
function out = wavScale( wt, wavSplit, varargin )
  % out = wavScale( wt, wavSplit [, 'range', range ] );
  %
  % Written by Nicholas Dwork - Copyright 2019
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 1
    disp( 'Usage:  out = wavScale( wt, wavSplit [, ''range'', range ] )' );
    return
  end

  p = inputParser;
  p.addParameter( 'range', [0 1] );
  p.parse( varargin{:} );
  range = p.Results.range;

  sWT = size( wt );
  wt11 = wt( 1 : sWT(1)/2, 1 : sWT(2)/2 );
  wt12 = wt( 1 : sWT(1)/2, sWT(2)/2+1 : end );
  wt21 = wt( sWT(1)/2+1 : end, 1 : sWT(2)/2 );
  wt22 = wt( sWT(1)/2+1 : end, sWT(2)/2+1 : end );

  nSplit = numel(wavSplit);
  if nSplit > 1
    sSplit = size(wavSplit);
    s11 = wavSplit( 1:sSplit(1)/2, 1:sSplit(2)/2 );
    s12 = wavSplit( 1:sSplit(1)/2, sSplit(2)/2+1:end );
    s21 = wavSplit( sSplit(2)/2+1:end, 1:sSplit(1)/2 );
    s22 = wavSplit( sSplit(2)/2+1:end, sSplit(2)/2+1:end );

    if sum( s11(:) ) > 0
      if max( mod(size(wt11),2) ) > 0
        error('wavShow: improper dimensions of image');
      end
      wt11 = wavScale( wt11, s11 );
    end
    if sum( s12(:) ) > 0
      if max( mod(size(wt12),2) ) > 0
        error('wavShow: improper dimensions of image');
      end
      wt12 = wavScale( wt12, s12 );
    end
    if sum( s21(:) ) > 0
      if max( mod(size(wt21),2) ) > 0
        error('wavShow: improper dimensions of image');
      end
      wt21 = wavScale( wt21, s21 );
    end
    if sum( s22(:) ) > 0
      if max( mod(size(wt22),2) ) > 0
        error('wavShow: improper dimensions of image');
      end
      wt22 = wavScale( wt22, s22 );
    end
  end

  scaled11 = scaleImg( wt11, range );
  scaled12 = scaleImg( wt12, range );
  scaled21 = scaleImg( wt21, range );
  scaled22 = scaleImg( wt22, range );

  out = [ scaled11 scaled12; scaled21 scaled22 ];
end
