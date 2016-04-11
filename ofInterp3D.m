function interped = ofInterp3D( data, du, dv, dw )
  [M N K] = size(data);
  [ xIndxs, yIndxs, zIndxs ] = meshgrid( 1:N, 1:M, 1:K );

  interped = interp3( xIndxs, yIndxs, zIndxs, data, ...
    xIndxs+du, yIndxs+dv, zIndxs+dw, 'linear', 0 );

  %nanIndxs = find( ~isfinite( interped ) );
  %if numel(nanIndxs) > 0
  %  interped(nanIndxs) = img(nanIndxs);
  %end
end
