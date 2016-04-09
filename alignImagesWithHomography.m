
function proj2 = alignImagesWithHomography( img1, img2, varargin )
  % proj2 = alignVolumesWithHomography( img1, img2 [, distThresh] )
  %
  % This function aligns img2 with img1 using a Homography.
  % The algorithm it employs is to:
  %   - find features in img1 with Difference of Gaussians
  %   - track features into img2 with normalized cross correlation
  %   - Use RANSAC with the Direct Linear Transformation to identify the
  %     3x3 homography H.
  %   - Project img2 into the image 1 space.
  %
  % Inputs:
  % img1 - a 2D array
  % img2 - a 2D array
  %
  % Optional Inputs:
  % distThresh - distance threshold for RANSAC  (default is 5)
  %
  % Outputs:
  % proj2 - the projected image 2
  %
  % Written by Nicholas Dwork Copyright 2016

  defaultDistThresh = 5;
  p = inputParser;
  p.addOptional( 'distThresh', defaultDistThresh );
  p.parse( varargin{:} );
  distThresh = p.Results.distThresh;

  [pts1, pts2] = findAndTrackCorners( img1, img2 );

  H21 = ransacDltHomographyFromPts2D( pts2, pts1, distThresh );

  proj2 = projectImage( img2, H21 );

end
