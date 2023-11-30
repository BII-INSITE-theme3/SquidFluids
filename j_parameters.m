% Parameters needed to run j_full2D

% Parameters
rho = 10; % Number of stokeslets per unit length.
Lt = 8; % Length of top segment.
Lm = 8; % Length of transition region.
Lb = 10; % Length of bottom segment.
theta = pi/4; % Angle of right transition region to horizontal.
Ltot = Lt+sin(pi/4)*Lm+Lb;
Ptx = 10; % Position of top point of right boundary.
Pty = Lt+Lm/2; % Position of top point of right boundary.
dsep = 0.5; % Appendage separation.
psi = 0; % Angle of inclination between appendage pairs.
PRAx = 5; % Position of right appendage in x.
PRAy = 4; % Position of right appendage in y.
nptx = 100;
npty = nptx;
x = linspace(-(Ptx+1),(Ptx+1),nptx);
y = linspace(-(Ptx+1),(Ptx+1),npty);
