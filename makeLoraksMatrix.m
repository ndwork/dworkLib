
function out = makeLoraksMatrix( in, varargin )
  % out = makeLoraksMatrix( in [, 'op', op, 'R', R, 'sKData', sKData ] );
  %
  % Makes the C LORAKS matrix from the k-space array in.
  % The C LORAKS matrix is defined in "Low-Rank Modeling of Local k-Space
  % Neighborhoods (LORAKS) for Constrained MRI" by Justin Haldar
  %
  % Inputs:
  % in - if 'notransp' then a 2D x nCoils complex array specifying the k-space values
  %    - if 'transp' then a LORAKS matrix that must be converted to a 2D array
  %
  % Optional Inputs:
  % op - either 'notransp' (the default) or 'transp'
  %
  % Written by Nicholas Dwork, Copyright 2020

  p = inputParser;
  p.addParameter( 'op', 'notransp', @(x) true );
  p.addParameter( 'R', 3, @ispositive );
  p.addParameter( 'sKData', [], @(x) isnumeric(x) );
  p.parse( varargin{:} );
  op = p.Results.op;
  R = p.Results.R;
  sKData = p.Results.sKData;

  if strcmp( op, 'notransp' )
    sKData = size( in );
  else
    if numel( sKData ) == 0, error( 'Must supply sKData for transp operation' ); end
  end
  if numel( sKData ) > 2
    nCoils = sKData(3);
  else
    nCoils = 1;
  end

  rImg = makeRadialImg( [ 2*R+1 2*R+1 ] ) <= R;
  K = sum( rImg(:) );
  NR = ( sKData(1) - 2*R ) * ( sKData(2) - 2*R );

  L = @(x) makeLoraksCoilMatrix( x, 'notransp', R, sKData(1:2), rImg, K, NR );
  Lh = @(x) makeLoraksCoilMatrix( x, 'transp', R, sKData(1:2), rImg, K, NR );

  doCheckAdjoint = false;
  if doCheckAdjoint == true
    tmp = zeros( sKData(1:2) );
    [checkL,adjErrL] = checkAdjoint( tmp, L, Lh );
    if checkL ~= true
      error([ 'L and Lh failed ajoint test with error ', num2str(adjErrL) ]);
    end
  end

  if strcmp( op, 'notransp' )
    out = cell( 1, nCoils );
    %out = zeros( K, NR * nCoils );

    parfor coilIndx = 1 : nCoils
      %startIndx = (coilIndx-1) * NR + 1;
      %endIndx = coilIndx * NR;
      %out( :, startIndx : endIndx ) = L( in(:,:,coilIndx) );
      %out( :, startIndx : endIndx ) = makeLoraksCoilMatrix( in, 'notransp', ...
      %  R, sKData(1:2), rImg, K, NR );
      out{ coilIndx } = L( in(:,:,coilIndx ) );
    end
    out = cell2mat( out );

  elseif strcmp( op, 'transp' )
    out = cell(1,1,nCoils);
    %out = zeros( sKData );
    parfor coilIndx = 1 : nCoils
      startIndx = (coilIndx-1) * NR + 1;
      endIndx = coilIndx * NR;
      thisLoraksMatrix = in( :, startIndx : endIndx );
      %out( :, :, coilIndx ) = Lh( thisLoraksMatrix );
      out{1,1,coilIndx} = Lh( thisLoraksMatrix );
      %out( :, :, coilIndx ) = makeLoraksCoilMatrix( thisLoraksMatrix, 'transp', ...
      %  R, sKData(1:2), rImg, K, NR );
    end
    out = cell2mat( out );
    
  else
    error( 'Unrecognized value of op' );
  end

end

function out = makeLoraksCoilMatrix( in, op, R, sKData, rImg, K, NR )

  if strcmp( op, 'notransp' )    
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
    [ kIn, nrIn ] = size( in );
    if kIn ~= K, error( 'Incorrect input dimensions' ); end
    if nrIn ~= NR, error( 'Incorrect input dimensions' ); end

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



