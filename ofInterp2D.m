
function interped = ofInterp2D( img, du, dv )
  [M N] = size(img);
  xIndxs = ones(M,1) * (1:N);
  yIndxs = (1:M)' * ones(1,N);
  
  interped = interp2( xIndxs, yIndxs, img, ...
    xIndxs+du, yIndxs+dv, 'linear', 0 );
  
  %nanIndxs = find( ~isfinite( interped ) );
  %if numel(nanIndxs) > 0
  %  interped(nanIndxs) = img(nanIndxs);
  %end
end
