
function out = makeLoraksMatrix( in, varargin )
  % out = makeLoraksMatrix( in [, 'op', op, 'R', R, 'sKData', sKData ] );
  %
  % Makes the C LORAKS matrix from the k-space img in.
  % The C LORAKS matrix is defined in "Low-Rank Modeling of Local k-Space
  % Neighborhoods (LORAKS) for Constrained MRI" by Justin Haldar
  %
  % Inputs:
  % in - if 'notransp' then a 2D complex array specifying the k-space values
  %    - if 'transp' then a LORAKS matrix that must be converted to a 2D array
  %
  % Optional Inputs:
  % op - either 'notransp' (the default) or 'transp'
  %
  % Written by Nicholas Dwork, Copyright 2020

  p = inputParser;
  p.addParameter( 'op', 'notransp', @(x) true );
  p.addParameter( 'R', 3, @ispositive );
  p.addParameter( 'sKData', [], @(x) isnumeric(x) && numel(x) == 2 );
  p.parse( varargin{:} );
  op = p.Results.op;
  R = p.Results.R;
  sKData = p.Results.sKData;

  rImg = makeRadialImg( [ 2*R+1 2*R+1 ] ) <= R;
  K = sum( rImg(:) );

  if strcmp( op, 'notransp' )

    sKData = size( in );
    NR = ( sKData(1) - 2*R ) * ( sKData(2) - 2*R );
    out = zeros( K, NR );

    cIndx = 1;
    for i = R+1 : sKData(2)-R
      for j = R+1 : sKData(1)-R
        subKData = in( i-R:i+R, j-R:j+R );
        out(:,cIndx) = subKData( rImg == 1 );
        cIndx = cIndx + 1;
      end
    end

  elseif strcmp( op, 'transp' )
    if numel( sKData ) == 0
      error( 'Must supply sKData for transp operation' );
    end
    [ kIn, NR ] = size( in );   %#ok<ASGLU>
    if kIn ~= K, error( 'Incorrect input dimensions' ); end

    out = zeros( sKData );
    cIndx = 1;
    tmp = zeros( 2*R+1 );
    for i = R+1 : sKData(2)-R
      for j = R+1 : sKData(1)-R
        tmp( rImg == 1 ) = in(:,cIndx);
        out(i-R:i+R,j-R:j+R) = out(i-R:i+R,j-R:j+R) + tmp;
        cIndx = cIndx + 1;
      end
    end

  else
    error( 'Unrecognized value of op' );
  end

end
