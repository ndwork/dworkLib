
function cb = colorbarnice( varargin )
  % cb = colorbarnice( [ 'FontSize', FontSize, 'Label', Label', ...
  %   'LineWidth', LineWidth ] );
  %
  % Inputs:
  % FontSize - (default is 20)
  % Label - the label of the colorbar (default is no label)
  % LineWidth - (default is 2)
  %
  % Outputs:
  % cb - the colorbar's handle
  %
  % Written by Nicholas Dwork - Copyright 2017
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  defaultFontSize = 20;
  defaultLineWidth = 2;
  p = inputParser;
  p.addParameter( 'FontSize', defaultFontSize, @isnumeric );
  p.addParameter( 'Label', [] );
  p.addParameter( 'LineWidth', defaultLineWidth, @isnumeric );
  p.parse( varargin{:} );
  fontSize = p.Results.FontSize;
  cbLabel = p.Results.Label;
  lineWidth = p.Results.LineWidth;

  % Make space for the colorbar
  cf = gcf;  beforeFigUnits = cf.Units;  set( cf, 'units', 'pixels' );
  ca = gca;  beforeAxesUnits = ca.Units;  set( ca, 'units', 'pixels' );

  cfPos = get( cf, 'position' );  % Position is [ left bottom width height ]
  caPos = get( ca, 'position' );

  cb = colorbar( 'eastoutside' );
  cb.FontSize = fontSize;
  cb.LineWidth = lineWidth;
  beforeCbUnits = cb.Units;  set( cb, 'units', 'pixels' );

  set( cf, 'position', cfPos );
  set( ca, 'position', caPos );

  cbPos = get( cb, 'position' );

  beforeLabelUnits = cb.Label.Units;  set( cb.Label, 'units', 'pixels' );
  cbLabelPos = get( cb.Label, 'position' );
  cbLabelBuffer = 0;
  if numel(cbLabel) > 0
    cb.Label.FontSize = fontSize;
    cb.Label.String = cbLabel;
    cbLabelBuffer = 30;
  end

  newCfPos = [  cfPos(1), cfPos(2), cbPos(1)+cbLabelPos(1)+cbLabelBuffer, cfPos(4) ];
  set( cf, 'position', newCfPos );
  set( ca, 'position', caPos );

  set( cf, 'units', beforeFigUnits );
  set( ca, 'units', beforeAxesUnits );
  set( cb, 'units', beforeCbUnits );
  set( cb.Label, 'units', beforeLabelUnits );
end
