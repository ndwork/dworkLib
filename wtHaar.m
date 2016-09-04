
function wt = wtHaar( sig, varargin )

  defaultSplit = 1;
  p = inputParser;
  p.addOptional( 'split', defaultSplit );
  p.parse( varargin{:} );
  split = p.Results.split;

  wt1 = 1/sqrt(2) * ( sig(1:2:end) + sig(2:2:end) );
  wt2 = 1/sqrt(2) * ( sig(1:2:end) - sig(2:2:end) );

  nSplit = numel(split);
  if nSplit > 1
    s1 = split(1:nSplit/2);
    s2 = split(nSplit/2+1:end);

    if sum(s1)>0
      wt1 = wtHaar( wt1, s1 );
    end    
    if sum(s2)>0
      wt2 = wtHaar( wt2, s2 );
    end
  end

  wt = [wt1 wt2];
end
