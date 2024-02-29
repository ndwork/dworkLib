
function [out, phaseImg] = mri_reconPartialFourier( in, sFSR, varargin )
  % out = mri_reconPartialFourier( in, sFSR [, 'phases', phases, 'op', op ] )
  %
  % Written according to "Partial k-space Reconstruction" by John Pauly
  %
  % Inputs:
  % in - the input array of size Ny x Nx x nCoils representing the MRI data.
  % sFSR - Either a scalar or a two element array that specifies the size of the fully
  %   sampled region, which is used to estimate the phase of the image.
  %   If a scalar, then the size of the fully sampled region is sFSR x Nx.
  %   If a two element array, then the size of the FSR is sFSR(1) x sFSR(2).
  %
  % Outputs:
  % out - an array of size Ny x Nx x nCoils that represents the coil images
  %
  % Written by Nicholas Dwork - Copyright 2023
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addParameter( 'op', 'notransp', @(x) true );
  p.addParameter( 'phases', [], @isnumeric );
  p.addParameter( 'traj', [], @isnumeric );
  p.parse( varargin{:} );
  phases = p.Results.phases;
  op = p.Results.op;
  traj = p.Results.traj;

  if numel( traj ) == 0
    % Cartesian sampling
    [out, phaseImg] = mri_reconPartialFourier_Cartesian( in, sFSR, phases, op );

  else
    % Non-Cartesian sampling
    [out, phaseImg] = mri_reconPartialFourier_Cartesian( in, sFSR, phases, traj, op );
  end

end


function [out, phaseImg] = mri_reconPartialFourier_Cartesian( in, sFSR, phases, op )

  Ny = size( in, 1 );
  ys = size2imgCoordinates( Ny );
  centerIndx = find( ys == 0, 1 );

  if numel( sFSR ) == 1
    sFSR = [ sFSR, size( in, 2 ) ];
  end

  firstY = centerIndx - ceil( ( sFSR(1) - 1 ) / 2 );
  lastY = centerIndx + floor( ( sFSR(1) - 1 ) / 2 );
  m = 2 / ( lastY - firstY );
  ramp = m * ys + 1.0;
  ramp( 1 : firstY ) = 0;
  ramp( lastY : end ) = 2;
  ramp = 2 - ramp;

  nCoils = size( in, 3 );

  if numel( phases ) == 0
    if strcmp( op, 'transp' )
      error( 'Cannot estimate image phase during transpose operation' );
    end

    firstX = centerIndx - ceil( ( sFSR(2) - 1 ) / 2 );
    lastX = centerIndx + floor( ( sFSR(2) - 1 ) / 2 );

    inLF = in;  % Low-freq data
    if nCoils == 1
      inLF( 1:firstY-1, : ) = 0;
      inLF( :, 1:firstX-1 ) = 0;
      inLF( :, lastX+1:end ) = 0;
    else
      inLF( 1:firstY-1, :, : ) = 0;
      inLF( :, 1:firstX-1, : ) = 0;
      inLF( :, lastX+1:end, : ) = 0;
    end
    imgLF = fftshift2( ifft2( ifftshift2( inLF ) ) );
    phases = angle( imgLF );
  end
  phaseImg = exp( 1i * phases );

  if strcmp( op, 'notransp' )
    inW = bsxfun( @times, in, ramp );
    imgW = fftshift2( ifft2( ifftshift2( inW ) ) );
    out = imgW .* conj( phaseImg );
    out = real( out ) .* phaseImg;

  elseif strcmp( op, 'transp' )
    in = real( in .* conj( phaseImg ) );
    inWithPhase = in .* phaseImg;
    ifft2hIn = fftshift2( ifft2h( ifftshift2( inWithPhase ) ) );
    out = bsxfun( @times, ifft2hIn, ramp );

  else
    error( 'Unrecognized operator' );
  end
end



function [out, phaseImg] = mri_reconPartialFourier_NonCartesian( in, sFSR, phases, op, traj )
  error( 'I still need to implement this' )
end
