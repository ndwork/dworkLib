
function sig = iwtHaar( wt, varargin )

  defaultSplit = 1;
  p = inputParser;
  p.addOptional( 'split', defaultSplit );
  p.parse( varargin{:} );
  split = p.Results.split;

  nwt = numel(wt);
  wt1 = wt(1:nwt/2);
  wt2 = wt(nwt/2+1:end);

  nSplit = numel(split);
  if nSplit > 1
    s1 = split(1:nSplit/2);
    s2 = split(nSplit/2+1:end);

    if sum(s1)>0
      wt1 = iwtHaar( wt1, s1 );
    end
    if sum(s2)>0
      wt2 = iwtHaar( wt2, s2 );
    end
  end

  sig1 = upsample( wt1, 2 );
  sig1(2:2:end) = wt1;

  sig2 = upsample( wt2, 2 );
  sig2 = [sig2(1) sig2(2:end) - sig2(1:end-1)];

  sig = 1/sqrt(2) * ( sig1 + sig2 );
end
