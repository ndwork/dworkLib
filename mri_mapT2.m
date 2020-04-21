
function [t2Map,m0Map,RSquared] = mri_mapT2( dataCube, TEs, varargin )
  % [t2Map,m0Map] = mri_mapT2( dataCube, TEs, [, 'alg', alg, ...
  %  'b1ScaleMap', b1ScaleMap, 'mask', mask, 'verbose', verbose ] )
  %
  %  Determines T2 values by fitting data to a model.
  
  %  This objective function is taken from "Errors in the Measurements of T2 Using
  %  Multiple-Echo MRI Techniques by Majumdar et al. (1986)
  %
  % Fits T2 to an exponential decay.  If the data is complex, uses the
  % magnitude of the data.
  %
  % Inputs:
  % dataCube - a 3D array of size MxNxK.  Each k index corresponds to an
  %   image taken with a different spin echo time.
  % TEs - a 1D array of size K specifying the spin echo times of each image
  %
  % Optional Inputs:
  % alg - either 'lsqr' or 'linear' (linear fits line to log of data)
  % mask - a 2D array of size MxN.  Only map pixels with nonzero mask values.
  % verbose - scalar; info statements made if verbose is nonzero
  %
  % Outputs:
  % t2Map - a 2D array of size MxN; units are same as TEs
  % m0Map - a 2D array of size MxN; values are proportional to M(0)
  %  Note: M(0) is only related to proton density M0 if in equilibrium at time 0
  %
  % Written by Nicholas Dwork - Copyright 2018
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addParameter( 'alg', 'lsqr', @(x) true );
  p.addParameter( 'b1ScaleMap', [], @isnumeric );
  p.addParameter( 'mask', [], @(x) isnumeric(x) || islogical(x) );
  p.addParameter( 'verbose', 0, @(x) isnumeric(x) || islogical(x) );
  p.parse( varargin{:} );
  b1ScaleMap = p.Results.b1ScaleMap;
  alg = p.Results.alg;
  mask = p.Results.mask;
  verbose = p.Results.verbose;

  sData = size( dataCube );
  if numel( mask ) == 0, mask=ones( sData(1:2) ); end
  TEs = TEs(:);

  signalModel = @(x,tmp) tmp(1) .* exp(-x./tmp(2));

  if numel( b1ScaleMap ) > 0 && strcmp( alg, 'lsqr' )
    alg = 'lsqr_non180';
  end

  if ~strcmp( alg, 'linear' ) && nargout > 2
    error( 'RSquared only available with linear fit' );
  end

  fminconOptions = optimoptions('fmincon','Display','off');
  m0MapCols = cell( 1, sData(2) );
  t2MapCols = cell( 1, sData(2) );
  RSquaredCols = cell( 1, sData(2) );

  p = parforProgress( sData(2) );
  parfor i=1:sData(2)
    if verbose ~= 0, p.progress( i, 10 ); end    %#ok<PFBNS>

    m0MapCol = zeros( sData(1), 1 );  %#ok<PFBNS>
    t2MapCol = zeros( sData(1), 1 );
    RSquaredCol = zeros( sData(1), 1 );

    for j=1:sData(1)
      if mask(j,i)==0, continue; end

      thisData = abs( squeeze( dataCube(j,i,:) ) );

      if strcmp( 'linear', alg )
        [params,thisRSquared] = fitPolyToData( 1, TEs, thisData );
        %figure; scatternice( thisData );  hold on; plotnice( evaluatePoly(params,TEs), 'r' )
        thisT2 = -1 / params(2) * params(1);
        m0MapCol(j) = params(1);
        t2MapCol(j) = thisT2;
        RSquaredCol(j) = thisRSquared;

      elseif strcmp( 'lsqr', alg )
        params = fmincon( @(tmp) norm( thisData - signalModel(TEs,tmp) ), ...
          [0; 0;], [], [], [], [], [0; 0;], [], [], fminconOptions );
        m0MapCol(j) = params(1);
        t2MapCol(j) = params(2);

      elseif strcmp( 'lsqr_non180', alg )
        spinAngle = pi * b1ScaleMap(j,i);
        params = fmincon( ...
          @(tmp) norm( thisData - signalModel_non180(TEs,tmp,spinAngle) ), ...
            [0; 0;], [], [], [], [], [0; 0;], [], [], fminconOptions );
          %@(tmp) norm( thisData1 - f(TEs,tmp) ) + norm( thisData2 - f(TEs,tmp) ), ...
          %  [0; 0;], [], [], [], [], [0; 0;], [], [], fminconOptions );
        m0MapCol(j) = params(1);
        t2MapCol(j) = params(2);
        
      else
        error( 'Unrecognized algorithm specification' );

      end
    end

    m0MapCols{i} = m0MapCol;
    t2MapCols{i} = t2MapCol;
    if strcmp( 'linear', alg )
      RSquaredCols{i} = RSquaredCol;
    end
  end
  p.clean;

  m0Map = cell2mat( m0MapCols );
  t2Map = cell2mat( t2MapCols );
  if strcmp( 'linear', alg )
    RSquared = cell2mat( RSquaredCols );
  end
end


function out = signalModel_non180( TEs, params, angle )
  % This objective function is taken from "Errors in the Measurements of T2 Using
  % Multiple-Echo MRI Techniques by Majumdar et al. (1986)
  M0 = params(1);
  T2 = params(2);

  ca = cos(angle);
  sa = sin(angle);
  saSq = sa * sa;
  tmp1 = 1 - ca;
  tmp1Sq = tmp1 * tmp1;
  tmp2 = 1 + ca;

  nTEs = numel( TEs );
  out = zeros( numel(TEs), 1 );

  out(1) = 0.5 * tmp1;
  out(2) = 0.25 * tmp1*tmp1;

  if nTEs > 2
    out(3) = 0.125 * ( tmp1*tmp1Sq + tmp1*tmp2*tmp2 ) + 0.5 * ca*saSq;
  end

  if nTEs > 3
    out(4) = 0.0625 * tmp1Sq*tmp1Sq + 0.125 * 3 * saSq*saSq + 0.5*ca*saSq*tmp1;
  end

  if nTEs > 4, error( 'Not yet implemented' ); end

  out = out .* M0 .* exp(-TEs./T2);
end
