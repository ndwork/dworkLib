
function showFeaturesOnImg( features, img )

  imshow( img, [] );

  [nFeatures,~] = size(features);
  for i=1:nFeatures
    y = features(i,1);
    x = features(i,2);

    text( x-2, y-2, num2str(i), 'color', 'blue', 'FontSize', 20 );
    text( x, y, num2str(i), 'color', 'yellow', 'FontSize', 20 );
  end

end
