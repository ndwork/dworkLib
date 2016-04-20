
function showFeaturesOnImg( features, img )
  % showFeaturesOnImg( features, img )
  %
  % Inputs:
  % features - 2D array of size Nx2
  %   N is the number of points
  %   The first/second column is the x/y location
  %
  % Written by Nicholas Dwork - Copyright 2016

  figure;
  imshow(img, []);
  hold on
  plot( features(:,1), features(:,2), 'r*');
  drawnow;
end
