
function showFeaturesOnImg( features, img )

  figure;
  imshow(img, []);
  hold on
  plot( features(:,2), features(:,1), 'r*');
  drawnow;
end
