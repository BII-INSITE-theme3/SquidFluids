%%

%% This code is shit, make it better.
% Boundary flow field.
%U0(:,1) = 1; % Set the values (this is where the code is modified).
Ul = zeros(Nstoks/2,2);
Urtemp = zeros(Nstoks/2,2);
Ur = zeros(Nstoks/2,2);

coef = zeros(80,1);
coef2 = zeros(80,1);
coef(1:60) = sin((0:59)*2*pi/60);
coef = abs(coef);

for i = 1:80
    coef2(i) = coef(mod(i-20,80)+1);
end

%plot(coef2)

for i = 1:80

    Ul(i,1) = -sin(i*2*pi/80)*coef2(i);
    Ul(i,2) = cos(i*2*pi/80)*coef2(i);
    Urtemp(end-i+1,1) = sin(i*2*pi/80)*coef2(end-i+1);
    Urtemp(end-i+1,2) = -cos(i*2*pi/80)*coef2(end-i+1);

end

Ur = Urtemp;

% for i = 1:80
% 
%     ii = 1 + mod(i+20,80);
%     Ur(ii,1) = Urtemp(i,1);
%     Ur(ii,2) = Urtemp(i,2);
% 
% end

U0 = [Ul;Ur];

% No net force in x on either
%U0(:,1) = U0(:,1);
%U0(:,2) = U0(:,2) - 10; % background flow

%quiver(cos(theta'),sin(theta'),Ul(:,1),Ur(:,2)); hold on
%plot(Ur)

%%

clear Ux Uy

p = stks(10,:); % Get the position of consideration
Ux = 0; Uy = 0;

for n = 1:Nstoks

    stksN = stks(n,:); % Get the position of stokeslet N.
    Ftemp = F(n,:); % Get the forces of stokeslet N.

    r = sqrt(norm(p-stksN) + eps^2) + eps; % Distance, considered to stokeslet N.
    rho = (r+eps)/(r*(r-eps)); % Rho, considered to stokeslet N.

    delta(1) = p(1)-stksN(1);
    delta(2) = p(2)-stksN(2);

    for k = 1:2
        for l = 1:2
            Stemp(k,l) = -(log(r)-eps*rho)*(k==l) - (delta(k)*delta(l)*rho/r);
        end
    end

    U = Stemp*Ftemp';
    Ux = Ux + U(1);
    Uy = Uy + U(2);

end

disp([Ux Uy]);

%%



Ux = zeros(Nstoks,1);
Uy = zeros(Nstoks,1);

eps = 0.001;

f = 0;

for j = 1:Nstoks

    p = stks(j,:);

    for n = 1:Nstoks

        pN = stks(n,:); % Get the position of stokeslet N.
        Ftemp = F(n,:); % Get the forces of stokeslet N.
        r = sqrt(norm(p' - pN').^2 + eps^2) + eps; % Distance, considered to stokeslet N.
        rho = (r+eps)/(r*(r-eps)); % Rho, considered to stokeslet N.

        for k = 1:2
            for l = 1:2
                Stemp(k,l) = -(log(r)-eps*rho)*(k==l) + (p(k)-pN(k))*(p(l)-pN(l))*rho/r;
            end
        end

        U = Stemp*Ftemp';
        Ux(j) = Ux(j) + U(1);
        Uy(j) = Uy(j) + U(2);

    end

end


n=10;
%scatter(stks(:,1),stks(:,2))
hold on
quiver(stks(:,1),stks(:,2),Ux,Uy,2)
hold on;