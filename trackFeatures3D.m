
function pts2 = trackFeatures3D( pts1, vol1, vol2, varargin )
  % pts2 = trackFeatures( pts1, vol1, vol2 [, 'searchWidth', searchWidth, ...
  %   'kernelWidth', kernelWidth, 'offset', offset] )
  %
  % Algorithm tracks points from img1 into img2 using Normalized Cross
  % correlation according to "Fast Normalized Cross Correlation" by Lewis.
  %
  % Inputs:
  % pts1 - an Nx3 array specifying the y/x/z location of each feature in
  %   vol1.  It has size Nx2 where N is the number of features
  % vol1 - a 3D array representing the first volume
  % vol2 - a 3D array representing the second volume
  %   features will be tracked into this volume
  %
  % Optional Inputs:
  % searchWidth - a 3 element array specifying the size of the search area
  %   If a scalar is supplied, it assumed the same for all dimensions
  %   by default, the search area is 20% of the image size
  % kernelWidth - a 3 element array specifying the size of the kernel
  %   If a scalar is supplied, it assumed the same for all dimensions
  % offset - a 3 element array specifying the y/x/z shift to center the
  %   search area (default is [0 0])
  %
  % Output:
  % pts2 - an Nx3 array specifying the y/x/z location of each feature
  %   in vol2.

  sVol = size( vol1 );
  nPts = size( pts1, 1 );

  defaultSearch = ceil( 0.2 * sVol );
  defaultKernel = 11;
  p = inputParser;
  p.addParamValue( 'searchWidth', defaultSearch );
  p.addParamValue( 'kernelWidth', defaultKernel );
  p.addParamValue( 'offset', [0 0 0] );
  p.parse( varargin{:} );
  offset = p.Results.offset;
  kw = p.Results.kernelWidth;
  sw = p.Results.searchWidth;

  if numel(kw)==1, kw=[kw kw kw]; end;
  if numel(sw)==1, sw=[sw sw sw]; end;
  hkw = floor( kw/2 );  % half kernel width
  hsw = floor( sw/2 );  % half search width

  pts2 = zeros(nPts,2);
  for i=1:nPts
    y = pts1(i,1);
    x = pts1(i,2);
    z = pts1(i,3);
    if y-hkw(1) < 1 || y+hkw(1) > sVol(1) || ...
       x-hkw(2) < 1 || x+hkw(2) > sVol(2) || ...
       z-hkw(3) < 1 || z+hkw(3) > sVol(3)
      continue;
    end;
    template = vol1(y-hkw:y+hkw,x-hkw:x+hkw,z-hkw:z+hkw);

    sB = max( y-hsw(1)+offset(1), 1 );  % search bottom
    sT = min( y+hsw(1)+offset(1), sVol(1) );  % search top
    sL = max( x-hsw(2)+offset(2), 1 );  % search left
    sR = min( x+hsw(2)+offset(2), sVol(2) );  % search right
    sA = max( z-hsw(3)+offset(3), 1 );  % search anterior
    sP = min( z+hsw(3)+offset(3), sVol(3) );  %search posterior
    search = vol2(sB:sT,sL:sR,sA:sP);
    if numel(search) == 0, continue; end;
    if min(size(search) - size(template)) < 2, continue; end;

    ncc = normxcorr3(template, search, 'same');
    [maxNCC,maxNccIndx] = max( ncc(:) );  %#ok<ASGLU>
    sSearch = size(search);
    [nccY,nccX,nccZ] = ind2sub( sSearch, maxNccIndx );
    pts2(i,1) = nccY + sB - 1;
    pts2(i,2) = nccX + sL - 1;
    pts2(i,3) = nccZ + sA - 1;
  end

end
