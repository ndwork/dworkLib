
function img = iwtHaar2( wt, varargin )
  % img = iwtHaar2( wt )

  defaultSplit = 1;
  p = inputParser;
  p.addOptional( 'split', defaultSplit );
  p.parse( varargin{:} );
  split = p.Results.split;
  
  % H/L - High / Low pass filter
  % h/v - horizontal / vertical direction

  sWT = size(wt);
  wt11 = wt(1:sWT(1)/2,1:sWT(2)/2);
  wt21 = wt(sWT(1)/2+1:end,1:sWT(2)/2);
  wt12 = wt(1:sWT(1)/2,sWT(2)/2+1:end);
  wt22 = wt(sWT(1)/2+1:end,sWT(2)/2+1:end);


  nSplit = numel(split);
  if nSplit > 1
    sSplit = size(split);
    s11 = split(1:sSplit(1)/2,1:sSplit(2)/2);
    s12 = split(1:sSplit(1)/2,sSplit(2)/2+1:end);
    s21 = split(sSplit(2)/2+1:end,1:sSplit(1)/2);
    s22 = split(sSplit(2)/2+1:end,sSplit(2)/2+1:end);

    if sum(s11)>0
      wt11 = iwtHaar( wt11, s11 );
    end
    if sum(s12)>0
      wt12 = iwtHaar( wt12, s12 );
    end
    if sum(s21)>0
      wt21 = iwtHaar( wt21, s21 );
    end
    if sum(s22)>0
      wt22 = iwtHaar( wt22, s22 );
    end
  end
  
  
  sW = size(wt11);

  recon11 = imresize( wt11, 2*sW, 'nearest' );

  mid12 = imresize( wt12, [2*sW(1) sW(2)], 'nearest' );
  mid12 = upsample2( mid12, [1 2] );
  recon12 = [mid12(:,1) mid12(:,2:end)-mid12(:,1:end-1)];

  mid21 = imresize( wt21, [sW(1) 2*sW(2)], 'nearest' );
  mid21 = upsample2( mid21, [2 1] );
  recon21 = [mid21(1,:); mid21(2:end,:)-mid21(1:end-1,:); ];

  mid22 = upsample2( wt22, 2 );
  mid22 = [mid22(1,:); mid22(2:end,:)-mid22(1:end-1,:);];
  recon22 = [mid22(:,1) mid22(:,2:end)-mid22(:,1:end-1)];

  img = 0.5 * (recon11 + recon12 + recon21 + recon22);
end

