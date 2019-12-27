
function testImages = listTestImages( varargin )
  % testImages = listTestImages()
  %
  % List the standard image processing images located in the testImages directory
  %
  % Outputs:
  % testImages - an array of structures for the image files in testImages
  %
  % Written by Nicholas Dwork - Copyright 2019
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  thisFilename = which( 'listTestImages' );

  thisFilenameParts = split( thisFilename, filesep() );

  thisFilenameParts{end-1} = 'testImages';
  testImagesDir = join( thisFilenameParts(1:end-1), filesep() );

  testImages = dir( testImagesDir{1} );
  testImages = testImages( 3 : end );

  testImgIndx = 1;
  nTestImages = numel( testImages );
  while testImgIndx <= nTestImages
    testImgName = testImages( testImgIndx ).name;
    if numel( regexp( testImgName, '^\.' ) ) > 0
      for i = testImgIndx+1 : nTestImages
        testImages( i-1 ) = testImages( i );
      end
      nTestImages = nTestImages - 1;
    else
      testImgIndx = testImgIndx + 1;
    end
  end

  testImages = testImages( 1 : nTestImages );

  if nargout < 1
    for testImgIndx = 1 : nTestImages
      disp( testImages( testImgIndx ).name );
    end
  end
end
