

function epgTest

  t = 50;
  T1 = 1000;
  T2 = 200;
  alpha = 35;

  Q = zeros(3, 10);
  Q(3,1) = 1;

  Q = epgRF( Q, alpha );
  Q = epgGrad( Q, 3 );
  Q = epgRelax( Q, t, T1, T2 );

  disp(Q);
end
