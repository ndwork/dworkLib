
function proj2 = alignVolumesWithHomography( vol1, vol2, varargin )
  % proj2 = alignVolumesWithHomography( vol1, vol2 [, distThresh] )
  %
  % This function aligns vol2 with vol1 using a Homography.
  % The algorithm it employs is to:
  %   - find features in vol1 with Difference of Gaussians
  %   - track features into vol2 with normalized cross correlation
  %   - Use RANSAC with the Direct Linear Transformation to identify the
  %     4x4 homography H.
  %   - Project vol2 into the volume 1 space.
  %
  % Inputs:
  % vol1 - a 3D array
  % vol2 - a 3D array
  %
  % Optional Inputs:
  % distThresh - distance threshold for RANSAC  (default is 5)
  %
  % Outputs:
  % proj2 - the projected volume 2
  %
  % Written by Nicholas Dwork Copyright 2016

  defaultDistThresh = 5;
  p = inputParser;
  p.addOptional( 'distThresh', defaultDistThresh );
  p.parse( varargin{:} );
  distThresh = p.Results.distThresh;

  [pts1, pts2] = findAndTrackFeatures3D( vol1, vol2, 'w', 21 );

  H21 = ransacDltHomographyFromPts3D( pts2, pts1, distThresh );

  proj2 = projectVolume( vol2, H21 );
  
end
