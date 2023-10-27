%% Script to calculate analytically the stokeslet cyclinder flow

% Space parameters.
Xsi = 1; Xlen = 500; % X dimensionless length, number of points, resp.
Ysi = 1; Ylen = 500; % Y dimensionless length, number of points, resp.
Xsys = linspace(-Xsi/2,Xsi/2,Xlen);
Ysys = linspace(-Ysi/2,Ysi/2,Ylen);

% Store for final velocity field.
UxSys = zeros(Xlen,Ylen);
UySys = zeros(Xlen,Ylen);

% Cylinder parameters.
a = 0.25;
f0 = (1/(1-2*log(a)))/20 * [1,0];

% loop through the space and solve.
for xi = 1:Xlen
    for yi = 1:Ylen

        R = [Xsys(xi),Ysys(yi)];
        r = norm(R);

        utemp = -f0*(2*log(r) - (a^2/r^2)) + 2*dot(f0,R)*R/r^2*(1-(a^2/r^2));

        UxSys(xi,yi) = utemp(1);
        UySys(xi,yi) = utemp(2);

    end
end

% Optional thresholding to improve the contour plots.
 thresh1 = 0.1;
 UxSys( abs(UxSys) > thresh1) = thresh1;
 UySys( abs(UySys) > thresh1) = thresh1;

%imagesc(UxSys)
%imagesc(UySys)

% thresh2 = 0.4;
% UxSys( abs(UxSys) < thresh2) = thresh2;
% UySys( abs(UySys) < thresh2) = thresh2;

figure
contour(UxSys',100)
% contour(UySys)
% axis equal
