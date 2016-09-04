
function wt = wtHaar2( img, varargin )
  % wt = wtHaar2( img [, split ] )

  defaultSplit = 1;
  p = inputParser;
  p.addOptional( 'split', defaultSplit );
  p.parse( varargin{:} );
  split = p.Results.split;

  % H/L - High / Low pass filter
  % h/v - horizontal / vertical direction

  wt1 = img(1:2:end-1,:) + img(2:2:end,:);
  wt11 = wt1(:,1:2:end-1) + wt1(:,2:2:end);
  wt12 = wt1(:,1:2:end-1) - wt1(:,2:2:end);

  wt2 = img(1:2:end-1,:) - img(2:2:end,:);
  wt21 = wt2(:,1:2:end-1) + wt2(:,2:2:end);
  wt22 = wt2(:,1:2:end-1) - wt2(:,2:2:end);

%   Lh = img(:,1:2:end-1) + img(:,2:2:end);
%   wt11 = Lh(1:2:end-1,:) + Lh(2:2:end,:);
%   wt21 = Lh(1:2:end-1,:) - Lh(2:2:end,:);
% 
%   Hh = img(:,1:2:end-1) - img(:,2:2:end);
%   wt12 = Hh(1:2:end-1,:) + Hh(2:2:end,:);
%   wt22 = Hh(1:2:end-1,:) - Hh(2:2:end,:);

  nSplit = numel(split);
  if nSplit > 1
    sSplit = size(split);
    s11 = split(1:sSplit(1)/2,1:sSplit(2)/2);
    s12 = split(1:sSplit(1)/2,sSplit(2)/2+1:end);
    s21 = split(sSplit(2)/2+1:end,1:sSplit(1)/2);
    s22 = split(sSplit(2)/2+1:end,sSplit(2)/2+1:end);

    if sum(s11)>0
      wt11 = wtHaar( wt11, s11 );
    end
    if sum(s12)>0
      wt12 = wtHaar( wt12, s12 );
    end
    if sum(s21)>0
      wt21 = wtHaar( wt21, s21 );
    end
    if sum(s22)>0
      wt22 = wtHaar( wt22, s22 );
    end
  end

  wt = 0.5 * [ wt11 wt12; wt21 wt22 ];

end
