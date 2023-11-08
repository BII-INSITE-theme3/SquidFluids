
close all

% Constants.
Xma = 10;
Yma = Xma;
Npts = 400; % Number of points in the "calculated" space.
eps = 0.00001; % Epsilon of the regularization.
Nstoks = 2*Npts; % Number of stokeslets on the radius.

% Preallocation.
x = linspace(-Xma/2,Xma/2,Npts); % X-values of solution space.
y = linspace(-Yma/2,Yma/2,Npts); % Y-values of solution space.
%
Ux = zeros(Npts); % fluid velocity x-component.
Uy = zeros(Npts); % fluid velocity y-component.

% Background flow
flowang = 0;
flowstr = -1;
Uflow = flowstr*[cos(flowang),sin(flowang)];
%
wall1 = zeros(length(x),1);
wall1(1:length(x)/2) = -1;
wall1 = 2*(wall1+0.5);
windowsize = 30;
wall1 = smooth(smooth(smooth(wall1,windowsize),windowsize),windowsize);
wall2 = -wall1;
wallwid = 1;
spacing = 3;
wall1 = wall1*wallwid - spacing;
wall2 = wall2*wallwid + spacing;
stks = [x,x; wall1',wall2']';

%plot(stks(:,2),stks(:,1))

%%
S = zeros(2*Nstoks); % Store for the Stokeslet.
Stemp = zeros(2); % Store for the Stokeslet.
F = zeros(Nstoks,2); % Store for the stokeslet's forces.
Ftemp = zeros(2,1); % Store for the stokeslet's forces.
%
U0 = zeros(Nstoks,2); % Store the for fluid velocity at the stokeslets.
U0V = zeros(1,2*Nstoks); % Vertical store the for fluid velocity at the stokeslets.

%% Boundary flow

U0(:,1) = -Uflow(1);
U0(:,2) = -Uflow(2);

%%

% Put the fluid velocity into "vertical array" form.
for i = 1:Nstoks
    U0V(2*i-1) = U0(i,1);
    U0V(2*i)   = U0(i,2);
end

% Find the forces.
for l = 1:Nstoks % Loop through the stokeslet-components
    for k = 1:Nstoks % Loop through the stokeslet-components
        %
        stks1 = stks(l,:); % Get the position of the influence stokeslet.
        stks2 = stks(k,:); % Get the position of the influenced stokeslet.
        %
        r = sqrt(norm(stks1-stks2)^2 + eps^2) + eps; % Get the reg `distance' between them.
        rho = (r+eps)/(r*(r-eps)); % Get the "rho", for convenience.
        %
        logterm = -(log(r)-eps*rho);
        x2term = (stks1(1)-stks2(1))*(stks1(1)-stks2(1))*rho/r;
        xyterm = (stks1(1)-stks2(1))*(stks1(2)-stks2(2))*rho/r; 
        y2term = (stks1(2)-stks2(2))*(stks1(2)-stks2(2))*rho/r;

        S(2*l-1,2*k-1) = logterm + x2term;
        S(2*l,  2*k-1)   = xyterm;
        S(2*l-1,2*k)     = xyterm;
        S(2*l,  2*k)     = logterm + y2term;
    end
end

FV = U0V/S; % Get the forces from the stokeslet.

% Put the forces velocity into vertical array form.
for i = 1:Nstoks
    F(i,1) = FV(2*i-1);
    F(i,2) = FV(2*i);
end

%% Solve over the space.

% Loop over the whole space.
for i = 1:Npts
    for j = 1:Npts

        p = [x(i),y(j)]; % Get the position of consideration, i and j done make sense to me here??

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

%%

n = 20;

UxTemp = Ux+Uflow(1);
UyTemp = Uy+Uflow(2);

nullflow = ones(Npts);
stTemp1 = stks(1:Npts,2);
stTemp2 = stks(Npts+1:end,2);

for i=1:Npts

    a = find(y(:) < stTemp1(i));
    b = find(y(:)>stTemp2(i));

    nullflow(i,a) = 0;
    nullflow(i,b) = 0;

end

UxTemp = UxTemp.*nullflow;
UyTemp = UyTemp.*nullflow;

% Optional thresholding to improve the contour plots.
% thresh1 = 20;
% UxTemp( abs(Ux) > thresh1) = thresh1;
% UyTemp( abs(Uy) > thresh1) = thresh1;

Umag = sqrt(UxTemp.^2 + UyTemp.^2);

figure
imagesc(x,y,Umag)
%set(gca,'YDir','normal')
hold on
%contour(x,y,UxTemp',n,'r')
%scatter(stks(:,2),stks(:,1),'k','LineWidth',3)
quiver(x(1:n:end),y(1:n:end),UyTemp(1:n:end,1:n:end),UxTemp(1:n:end,1:n:end),2,'r')

%%

% c = jet(Npts);
% 
% for h = 1:4:Npts
%     plot(x,UxTemp(h,:),'Color',c(h,:))
%     hold on
% end

%%

UxPois = UxTemp;
UyPois = UyTemp;

% save('outputs/800pts_Poissuelle')
save('outputs/U_800pts_Poissuelle',"UxPois","UyPois","x","y")
