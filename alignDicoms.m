
function [out1,out2] = alignDicoms( in1, in2, varargin )
  % [out1,out2] = alignDicoms( in1, in2 [, info1, info2 ] )
  %
  % Shifts and scales in2 to align it with in1
  %
  % Inputs:
  % in1 - either an image array or a character array specifying a dicom filename or a directory
  %       containing dicom files
  % in2 - either an image array or a character array specifying a dicom filename or a directory
  %       containing dicom files
  %       ( must be the same type as in1 )
  %
  % Optional Inputs:
  % info1 - the structure array associated with the dicom of data 1 (retrieved using dicominfo)
  %   if in1 is a numeric array, info1 must be supplied.
  % info2 - the structure array associated with the dicom of data 2 (retrieved using dicominfo)
  %   if in2 is a numeric array, info1 must be supplied.
  %
  % Outputs:
  % out2 - the image of in2 aligned with in1
  %
  % Written by Nicholas Dwork - Copyright 2019
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addOptional( 'info1', [] );
  p.addOptional( 'info2', [] );
  p.parse( varargin{:} );
  info1 = p.Results.info1;
  info2 = p.Results.info2;

  if isnumeric( in1 )
    if ~isnumeric( in2 ), error( 'in1 and in2 must be of the same type' ); end
    if numel( info1 ) == 0 || numel( info2 ) == 0
      error( 'If supplying image data, must also supply info data' );
    end
    data1 = in1;
    data2 = in2;

  elseif ischar( in1 )
    if ~ischar( in2 ), error( 'in1 and in2 must be of the same type' ); end

    if exist( in1, 'dir' )
      % in1 is the name of a directory containing dicom files
      dicomFiles1 = dir( [ in1, '/*.DCM' ] );
      nFiles1 = numel( dicomFiles1 );      
      file1 = [ in1, '/', dicomFiles1(1).name ];
      firstImg1 = double( dicomread( file1 ) );
      info1 = dicominfo( file1 );
      sImg1 = size( firstImg1 );
      data1 = zeros( sImg1(1), sImg1(2), nFiles1 );
      z1s = zeros( nFiles1, 1 );
      for fIndx = 1 : nFiles1
        data1(:,:,fIndx) = dicomread( [ in1, '/', dicomFiles1( fIndx ).name ] );
        thisInfo1 = dicominfo( [in1, '/', dicomFiles1( fIndx ).name ] );
        thisPos = thisInfo1.ImagePositionPatient;
        z1s(fIndx) = thisPos(3);
      end
      [~,sortedZIndxs] = sort( z1s );
      data1 = data1(:,:,sortedZIndxs);
      z1s = z1s( sortedZIndxs );

      dicomFiles2 = dir( [ in2, '/*.DCM' ] );
      nFiles2 = numel( dicomFiles2 );
      file2 = [ in2, '/', dicomFiles2(2).name ];
      firstImg2 = double( dicomread( file2 ) );
      info2 = dicominfo( file2 );
      sImg2 = size( firstImg2 );
      data2 = zeros( sImg2(1), sImg2(2), nFiles2 );
      z2s = zeros( nFiles2, 1 );
      for fIndx = 1 : nFiles2
        data2(:,:,fIndx) = dicomread( [in2, '/', dicomFiles2( fIndx ).name ] );
        thisInfo2 = dicominfo( [in2, '/', dicomFiles2( fIndx ).name ] );
        thisPos = thisInfo2.ImagePositionPatient;
        z2s(fIndx) = thisPos(3);
      end
      [~,sortedZIndxs] = sort( z2s );
      data2 = data2(:,:,sortedZIndxs);
      z2s = z2s( sortedZIndxs );

    else
      % in1 is the filename of a dicom file
      if ~exist( in1, 'file' ), error( 'in1 must be a filename or a directory name' ); end
      if ~exist( in2, 'file' ), error( 'in2 must be a filename or a directory name' ); end
      data1 = double( dicomread( in1 ) );
      data2 = double( dicomread( in2 ) );
      if numel( info1 ) == 0, info1 = dicominfo( in1 ); end
      if numel( info2 ) == 0, info2 = dicominfo( in2 ); end
      z1s = info1.ImagePositionPatient;  z1s=z1s(3);
      z2s = info1.ImagePositionPatient;  z2s=z2s(3);
    end

  end


  pxlSpacing1 = info1.PixelSpacing;     pos1 = info1.ImagePositionPatient;
  pxlSpacing2 = info2.PixelSpacing;     pos2 = info2.ImagePositionPatient;

  if max( pxlSpacing2 ~= pxlSpacing1 ) == 1
    % scale image 2 so that pixels have the same size as image 1
    hScale = pxlSpacing2(1) / pxlSpacing1(1);
    vScale = pxlSpacing2(2) / pxlSpacing1(2);
    newSize = round( size( data2(:,:,1) ) .* [ vScale hScale ] );

    scaled2 = zeros( newSize(1), newSize(2), size(data2,3) );
    for i = 1 : size( data2, 3 )
      scaled2(:,:,i) = imresize( data2(:,:,i), newSize, 'nearest' );
    end
    data2 = scaled2;
  end
  if size(data2,1) < size(data1,1)
    tmp = zeros( size(data1,1), size(data2,2), size(data2,3) );
    tmp(1:size(data2,1),:,:) = data2;
    data2 = tmp;
    clear tmp;
  end
  if size(data2,2) < size(data1,2)
    tmp = zeros( size(data2,1), size(data1,2), size(data2,3) );
    tmp(:,1:size(data2,2),:) = data2;
    data2 = tmp;
    clear tmp;
  end

  % Perform in-plane shifts
  dPos = pos2 - pos1;
  dShift = dPos(2:-1:1) ./ pxlSpacing1;
  shifted2 = zeros( size( data2 ) );
  for i = 1 : size( data2, 3 )
    shifted2(:,:,i) = shiftImg( data2(:,:,i), round( dShift ) );
  end
  shifted2 = shifted2( 1 : size(data1,1), 1 : size(data1,2), : );


  closest2 = zeros( size(data1) );
  newZ2s = zeros( size(data1,3) );
  sliceThickness1 = info1.SliceThickness;
  sliceThickness2 = info2.SliceThickness;
  for i = 1 : size( data1, 3 )
    % Find the closest slice in data2
    thisZ1 = z1s( i );
    [zDist,closestIndx] = min( abs( z2s - thisZ1 ) );
    if zDist < sliceThickness1 + sliceThickness2
      closest2(:,:,i) = shifted2( :, :, closestIndx );
      newZ2s(i) = z2s( closestIndx );
    end
  end


  out2 = closest2;
  out1 = data1;
end

