function [x,y] = poissonDisc2D(r,varargin)
  % POISSONDIC2D Generate Poisson-disc-sampled points in two dimensions
  % 
  %   [X,Y] = RANDPD2D(R) randomly generates a collection of points (X,Y) in
  %   rectangle (-0.5,-0.5),(0.5,0.5) that satisfy the minimum-distance 
  %   requirement specified by R using the Brindson method (without using the
  %   background grid).
  %   
  %   RANDPD2D uses the RAND function in its current/existing state.
  %   
  % Ethan Johnson, 2016
  %
  % Inspiration:
  %  + http://www.cs.ubc.ca/~rbridson/docs/bridson-siggraph07-poissondisk.pdf
  % 
  % Other reading:
  %  + https://www.jasondavies.com/poisson-disc/
  %  + https://bost.ocks.org/mike/algorithms/
  
  % Defaults
  k = 30; % limit of samples before rejection
  bx = [-0.5 0.5]; by = [-0.5 0.5]; % bounds
  
  % Parse
  if nargin>1 && ~isempty(varargin{1}), k = varargin{1}; end;
  if nargin>2 && ~isempty(varargin{2}), bx = sort(varargin{2}); end;
  if nargin>3 && ~isempty(varargin{3}), by = sort(varargin{3}); end;
  
  % Grid size, shortcuts
  r2 = r^2; t=3*r2; b = [by;bx].';
  
  % Make list of samples
  l = zeros(10,2); n = 0;
  
  % Make queue (kind of; really just a list b/c processed in random order)
  q = zeros(10,2); s = 0;
  
  % Initialise
  p0 = rand(1,2).*diff(b) + b(1,:);
  n = n+1; l(n,:) = p0; % save in list
  s = s+1; q(s,:) = p0; % add to queue
  
  % Add points
  while s>0
      % First check for size overrun
      if n==size(l,1), ll=zeros(2*n,2); ll(1:n,:)=l; l=ll; clear ll; end
      if s==size(q,1), qq=zeros(2*s,2); qq(1:s,:)=q; q=qq; clear qq; end
      % Then proceed to select point in the grid
      i = floor(rand*s)+1;
      % And generate candidate(s)
      c = false; % is it a candidate?
      for j=1:k % try at most k times
          aj = 2*pi*rand; rj = sqrt(rand*t + r2); % uniformly-random annulus
          xj = q(i,2) + rj*cos(aj); yj = q(i,1) + rj*sin(aj); % around point
          if ( bx(1)<xj&&xj<bx(2) && by(1)<yj&&yj<by(2) ...
               && ~tooprox( [yj xj], l(1:n,:), r ) ) % in bounds & not close?
              n = n+1; l(n,:) = [yj xj]; % save in list
              s = s+1; q(s,:) = [yj xj]; % add to queue
              c = true; % signal that candidate was found
              break; % quit making candidates
          end
      end
      if ~c, q(i,:) = q(s,:); s = s-1; end % remove point from queue
  end; clear i j aj rj xj yj c;
  
  % Return
  x = l(1:n,2); y = l(1:n,1);
end


function tp = tooprox(p,l,e)
  % TOOPROX Check for too-proximal point relative to list of points
  %  + P is the point as an n-tuple (row-vector)
  %  + L is the list of points with n-tuples across and different points down
  %  + E is the exclusion distance (2-norm of P-L)
  
  d = sum( bsxfun(@plus,-p,l).^2, 2);
  
  tp = any( d < (e^2) );
end

