
v = VideoWriter("R2Cent.mp4",'MPEG-4');
v.FrameRate = 10;
open(v)
iii = 0;

for cent2 = linspace(0.11,1,50)

iii = iii+1;

% Constants.
Xma = 2;
Yma = Xma;
Npts = 100; % Number of points in the "calculated" space.
eps = 0.00001; % Epsilon of the regularization.
Nstoks = 160; % Number of stokeslets on the radius.

% Preallocation.
x = linspace(-Xma/2,Xma/2,Npts); % X-values of solution space.
y = linspace(-Yma/2,Yma/2,Npts); % Y-values of solution space.
%
Ux = zeros(Npts); % fluid velocity x-component.
Uy = zeros(Npts); % fluid velocity y-component.

% Background flow
flowang = 0;
flowstr = 0.5;
Uflow = flowstr*[sin(flowang),cos(flowang)];
%
theta = linspace(0,2*pi,(Nstoks/2)+1); % Get the angles on the surface of the stokeslets.
theta = theta(1:end-1); % Remove the repeat value.
%
cent1 = -0.3;
%cent2 = -cent1;
R01 = 0.2; % Radius of the stokelet surface.
R02 = 0.2; % Radius of the stokelet surface.
stks1 = [cent1 + R01*cos(theta'),R01*sin(theta');]; % Get the stokeslet cartesian coordinates.
stks2 = [cent2 + R02*cos(theta'),R02*sin(theta');]; % Get the stokeslet cartesian coordinates.
stks = [stks1;stks2];
%

S = zeros(2*Nstoks); % Store for the Stokeslet.
Stemp = zeros(2); % Store for the Stokeslet.
F = zeros(Nstoks,2); % Store for the stokeslet's forces.
Ftemp = zeros(2,1); % Store for the stokeslet's forces.
%
U0 = zeros(Nstoks,2); % Store the for fluid velocity at the stokeslets.
U0V = zeros(1,2*Nstoks); % Vertical store the for fluid velocity at the stokeslets.

theta2 = linspace(0,2*pi,(Nstoks/2)+1);
theta2 = theta2(1:end-1);

modi = 1;

U0(1:Nstoks/2,1) = -modi*sin(theta2);
U0(1:Nstoks/2,2) = modi*cos(theta2);

U0(1:20,1) = 0;
U0(1:20,2) = 0;

U0(21:50,:) = -U0(21:50,:);

coef = zeros(80,1);
coef2 = zeros(80,1);
coef(1:60) = sin((0:59)*2*pi/60);
coef = abs(coef);
for i = 1:80
    coef2(i) = coef(mod(i-20,80)+1);
end

U0(1:80,1) = U0(1:80,1).*coef2;
U0(1:80,2) = U0(1:80,2).*coef2;

U0(1+(Nstoks/2):end,1) = modi*sin(theta2);
U0(1+(Nstoks/2):end,2) = -modi*cos(theta2);

U0(101:120,1) = 0;
U0(101:120,2) = 0;

U0(end-9:end,:) = -U0(end-9:end,:);
U0(81:100,:) = -U0(81:100,:);

for i = 1:80
    coef2(i) = coef(mod(i-40,80)+1);
end

U0(81:end,1) = U0(81:end,1).*coef2;
U0(81:end,2) = U0(81:end,2).*coef2;

% U0(:,1) = U0(:,1) - mean(U0(:,1));
% U0(:,2) = U0(:,2) - mean(U0(:,2));

U0(:,1) = U0(:,1) - Uflow(1);
U0(:,2) = U0(:,2) - Uflow(2);
%plot(coef2)
%n = 2;
%quiver(stks(1:n:end,1),stks(1:n:end,2),U0(1:n:end,1),U0(1:n:end,2))

% Put the fluid velocity into "vertical array" form.
for i = 1:Nstoks
    U0V(2*i-1) = U0(i,2);
    U0V(2*i)   = U0(i,1);
end

% Find the forces.
for l = 1:2*Nstoks % Loop through the stokeslet-components
    for k = 1:2*Nstoks % Loop through the stokeslet-components
        %
        n1 = ceil(k/2); % Get the stokelets number of the influence stokeslet.
        n2 = ceil(l/2); % Get the stokelets number of the influenced stokeslet.
        %
        stks1 = stks(n1,:); % Get the position of the influence stokeslet.
        stks2 = stks(n2,:); % Get the position of the influenced stokeslet.
        %
        r = sqrt(norm(stks1-stks2)^2 + eps^2) + eps; % Get the reg `distance' between them.
        rho = (r+eps)/(r*(r-eps)); % Get the "rho", for convenience.
        %
        S(l,k) = -(log(r)-eps*rho)*(mod(k,2)==mod(l,2)) ... % Get the log term contribution.
        + (stks1(1+mod(k,2))-stks2(1+mod(k,2)))*(stks1(1+mod(l,2))-stks2(1+mod(l,2)))*rho/r; % Get the <....> contribution.
    end
end

FV = U0V/S; % Get the forces from the stokeslet.

% Put the forces velocity into vertical array form.
for i = 1:Nstoks
    F(i,2) = FV(2*i-1);
    F(i,1) = FV(2*i);
end

% Loop over the whole space.
for i = 1:Npts
    for j = 1:Npts

        p = [x(j),y(i)]; % Get the position of consideration, i and j done make sense to me here??

        for n = 1:Nstoks

            pN = stks(n,:); % Get the position of stokeslet N.
            Ftemp = F(n,:); % Get the forces of stokeslet N.
            r = sqrt(norm(p - pN).^2 + eps^2) + eps; % Distance, considered to stokeslet N.
            rho = (r+eps)/(r*(r-eps)); % Rho, considered to stokeslet N.

            for k = 1:2
                for l = 1:2
                    Stemp(k,l) = -(log(r)-eps*rho)*(k==l) + (p(k)-pN(k))*(p(l)-pN(l))*rho/r;
                end
            end

            U = Stemp*Ftemp';
            Ux(i,j) = Ux(i,j) + U(1);
            Uy(i,j) = Uy(i,j) + U(2);

        end

    end
end

n = 3;

UxTemp = Ux+Uflow(1);
UyTemp = Uy+Uflow(2);

% Optional thresholding to improve the contour plots.
% thresh1 = 20;
% UxTemp( abs(Ux) > thresh1) = thresh1;
% UyTemp( abs(Uy) > thresh1) = thresh1;

Umag = sqrt(UxTemp.^2 + UyTemp.^2);

imagesc(x,y,Umag)
%set(gca,'YDir','normal')
hold on
%contour(x,y,UxTemp',n,'r')
plot(stks(1:80,1),stks(1:80,2),'k','LineWidth',3)
plot(stks(81:end,1),stks(81:end,2),'k','LineWidth',3)
quiver(x(1:n:end),y(1:n:end),UxTemp(1:n:end,1:n:end),UyTemp(1:n:end,1:n:end),2)
%saveas(gcf,['frames/flowAngle/' num2str(iii) '.png'])
frame = getframe(gcf);
writeVideo(v,frame);
hold off

end

close(v)
