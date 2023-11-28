% Title: Second attempt to do the inverse problem.
% Author: Stephen Williams.
% Notes: Attempt 2 without using Poissuelle flow field as background
%--------------------------------------------%

close all

%% Constants

% Constants.
Xma = 20; Yma = Xma;
Npts = 500; % Number of points in the "calculated" space.
eps = 0.00001; % Epsilon of the regularization.

% Preallocation.
x = linspace(-Xma/2,Xma/2,Npts); % X-values of solution space.
y = linspace(-Yma/2,Yma/2,Npts); % Y-values of solution space.
%
Ux = zeros(Npts); % fluid velocity x-component.
Uy = zeros(Npts); % fluid velocity y-component.

% Background flow
flowang = pi/2;
flowstr = 0.5;
Uflow = flowstr*[cos(flowang),sin(flowang)];
%

% Set the RS collection shapes

%

stks = []; % Initialise the store of positions
NstokBound = 500;
NstokCyl = 80;
theta = linspace(0,2*pi,NstokCyl+1); % Get the angles on the surface of the stokeslets.
theta = theta(1:end-1); % Remove the repeat value.

% Boundaries
wall1 = zeros(NstokBound,1);
wall1(1:NstokBound/2) = -1;
wall1 = 2*(wall1+0.5);
windowsize = 20;
wall1 = smooth(smooth(smooth(wall1,windowsize),windowsize),windowsize);
wall2 = -wall1;
wallwid = 1;
spacing = 5;
wall1 = wall1*wallwid - spacing;
wall2 = wall2*wallwid + spacing;

maxH = 5;
xwall = linspace(-maxH,maxH,NstokBound);

stks = [wall1',wall2';xwall,xwall]';

% Left cylinder pair
cent = -2;
cent1 = -0.5;
cent2 = -cent1;
R01 = 0.3; % Radius of the stokelet surface.
R02 = 0.3; % Radius of the stokelet surface.
H0 = -2; % Initial height
stks1 = [cent + cent1 + R01*cos(theta'),H0+R01*sin(theta');]; % Get the stokeslet cartesian coordinates.
stks2 = [cent + cent2 + R02*cos(theta'),H0+R02*sin(theta');]; % Get the stokeslet cartesian coordinates.

stks = [stks;stks1;stks2];

% Right cylinder pair
cent = -cent;
cent1 = -0.5;
cent2 = -cent1;
R01 = 0.3; % Radius of the stokelet surface.
R02 = 0.3; % Radius of the stokelet surface.
H0 = -2; % Initial height
stks1 = [cent + cent1 + R01*cos(theta'),H0+R01*sin(theta');]; % Get the stokeslet cartesian coordinates.
stks2 = [cent + cent2 + R02*cos(theta'),H0+R02*sin(theta');]; % Get the stokeslet cartesian coordinates.

stks = [stks;stks1;stks2];

%%

Nstoks = 2*NstokBound + 4*NstokCyl;

S = zeros(2*Nstoks); % Store for the Stokeslet.
Stemp = zeros(2); % Store for the Stokeslet.
F = zeros(Nstoks,2); % Store for the stokeslet's forces.
Ftemp = zeros(2,1); % Store for the stokeslet's forces.
%
U0 = zeros(Nstoks,2); % Store the for fluid velocity at the stokeslets.
U0V = zeros(1,2*Nstoks); % Vertical store the for fluid velocity at the stokeslets.

%% Set the flow speeds on the cylinders

% Nothing to do for 1st wall.

UcylL = zeros(NstokCyl,2);
UcylR = zeros(NstokCyl,2);

theta2 = linspace(-pi,pi,NstokCyl);

modi=-1;

UcylL(:,1) = -sin(theta2)*modi; UcylL(:,2) =  cos(theta2)*modi;

% UcylL(1:NstokCyl/4,1:2) = 0;
% UcylL(NstokCyl/4+1:NstokCyl*5/8,1:2) = -UcylL(NstokCyl/4+1:NstokCyl*5/8,1:2);
% 
% coef = zeros(NstokCyl,1);
% coef2 = zeros(NstokCyl,1);
% coef(1:NstokCyl*3/4) = sin((0:(NstokCyl*3/4)-1)*2*pi/(NstokCyl*3/4));
% coef = abs(coef);
% for i = 1:NstokCyl
%     coef2(i) = coef(mod(i-NstokCyl/4,NstokCyl)+1);
% end
% UcylL(:,1:2) = UcylL(:,1:2).*coef2;

UcylR(:,1) =  sin(theta2)*modi; UcylR(:,2) = -cos(theta2)*modi;

% UcylR(NstokCyl/4+1:NstokCyl/2,1:2) = 0;
% UcylR(1:NstokCyl/4,1:2) = -UcylR(1:NstokCyl/4,1:2);
% UcylR(end-(NstokCyl/8)+1:end,:) = -UcylR(end-(NstokCyl/8)+1:end,:);
% coef = zeros(NstokCyl,1);
% coef2 = zeros(NstokCyl,1);
% coef(1:NstokCyl*3/4) = sin((0:(NstokCyl*3/4)-1)*2*pi/(NstokCyl*3/4));
% coef = abs(coef);
% for i = 1:NstokCyl
%     coef2(i) = coef(mod(i-NstokCyl/2,NstokCyl)+1);
% end
% UcylR(:,1:2) = UcylR(:,1:2).*coef2;

Uwall1 = zeros(NstokBound,2);
Uwall2 = zeros(NstokBound,2);

UcylL1 = zeros(NstokCyl,2);
UcylL1(:,1:2) = [-sin(theta);cos(theta)]';
UcylL1(11:50,:) = -UcylL1(11:50,:);
UcylL1(1:20,:) = 0;

UcylR1 = zeros(NstokCyl,2);
UcylR1(:,1:2) = [sin(theta);-cos(theta)]';
UcylR1(1:30,:) = -UcylR1(1:30,:);
UcylR1(71:end,:) = -UcylR1(71:end,:);
UcylR1(21:40,:) = 0;

UcylL2 = zeros(NstokCyl,2);
UcylL2(:,1:2) = [-sin(theta);cos(theta)]';
UcylL2(11:50,:) = -UcylL2(11:50,:);
UcylL2(1:20,:) = 0;

UcylR2 = zeros(NstokCyl,2);
UcylR2(:,1:2) = [sin(theta);-cos(theta)]';
UcylR2(1:30,:) = -UcylR2(1:30,:);
UcylR2(71:end,:) = -UcylR2(71:end,:);
UcylR2(21:40,:) = 0;

U0 = [Uwall1;Uwall2;UcylL1;UcylR1;UcylL2;UcylR2;];

% U0(1:NstokBound,2) = -Uflow(2);
% evalue = NstokBound;
% U0(evalue+1:evalue+NstokBound,2) = -Uflow(2);
% evalue = evalue+NstokBound;
% U0(evalue+1:evalue+NstokCyl,1:2) = UcylL(:,1:2)-Uflow(2);
% evalue = evalue+NstokCyl;
% U0(evalue+1:evalue+NstokCyl,1:2) = UcylR(:,1:2)-Uflow(2);
% evalue = evalue+NstokCyl;
% U0(evalue+1:evalue+NstokCyl,1:2) = -UcylL(:,1:2)-Uflow(2);
% evalue = evalue+NstokCyl;
% U0(evalue+1:evalue+NstokCyl,1:2) = -UcylR(:,1:2)-Uflow(2);

Uflow(1) = Uflow(1);
Uflow(2) = Uflow(2);

% U0(:,1) = U0(:,1)-Uflow(1);
% U0(:,2) = U0(:,2)-Uflow(2);

%%

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

% Run to here for the forces.

%% Solve over the space.

% Loop over the whole space.
parfor i = 1:Npts
    for j = 1:Npts

        Stemp = zeros(2,2);
        tempStks = stks;
        tempF = F;
        tempX = x;
        tempY = y;

        p = [tempX(j),tempY(i)]; % Get the position of consideration, i and j done make sense to me here??

        for n = 1:Nstoks

            pN = tempStks(n,:); % Get the position of stokeslet N.
            Ftemp = tempF(n,:); % Get the forces of stokeslet N.
            r = sqrt(norm(p - pN).^2 + eps^2) + eps; % Distance, considered to stokeslet N.
            %r = sqrt(norm(p - pN).^2 + eps^2) + eps; % Distance, considered to stokeslet N.
            rho = (r+eps)/(r*(r-eps)); % Rho, considered to stokeslet N.

            for k = 1:2
                for l = 1:2
                    Stemp(k,l) = -(log(r)-eps*rho)*(k==l) + (p(k)-stks(n,k))*(p(l)-stks(n,l))*rho/r;
                    %Stemp(k,l) = -(log(r)-eps*rho)*(k==l) + (p(k)-pN(k))*(p(l)-pN(l))*rho/r;
                end
            end

            U = Stemp*Ftemp';
            Ux(i,j) = Ux(i,j) + U(1);
            Uy(i,j) = Uy(i,j) + U(2);

        end

    end
end

%%

n = 5;

x = linspace(-Xma/2,Xma/2,Npts); % X-values of solution space.
y = linspace(-Yma/2,Yma/2,Npts); % Y-values of solution space.

Uxtemp = Ux + Uflow(1);
Uytemp = Uy + Uflow(2);

% Optional thresholding to improve the contour plots.
% thresh1 = 20;
% UxTemp( abs(Ux) > thresh1) = thresh1;
% UyTemp( abs(Uy) > thresh1) = thresh1;

Umag = sqrt(Uxtemp.^2 + Uytemp.^2);

figure
imagesc(x,y,Umag)
%set(gca,'YDir','normal')
hold on
%contour(x,y,UxTemp',n,'r')
scatter(stks(:,1),stks(:,2),2,'r')
quiver(x(1:n:end),y(1:n:end),Uxtemp(1:n:end,1:n:end),Uytemp(1:n:end,1:n:end),2)
axis equal
