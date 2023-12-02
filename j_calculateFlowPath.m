function [Uflowx,Uflowy] = j_calculateFlowPath(stks1,F1,x1,y1)

    Uflowx = zeros(length(x1),1);
    Uflowy = zeros(length(x1),1);
    [nStok,~] = size(stks1);
    
    for xposition = 1:length(x1)
    
        Stemp = zeros(2,2);
        tempStks = stks1(:,1:2);
        tempF = F1;
        tempX = x1;
        tempY = y1;
    
        p = [tempX(xposition),tempY(xposition)]'; % Get the position of consideration, i and j done make sense to me here??
    
        for n = 1:nStok
    
            pN = tempStks(n,:)'; % Get the position of stokeslet N.
            Ftemp = tempF(n,:)'; % Get the forces of stokeslet N.
            r = sqrt(norm(p - pN).^2 + eps^2) + eps; % Distance, considered to stokeslet N.
            rho = (r+eps)/(r*(r-eps)); % Rho, considered to stokeslet N.
    
            for k = 1:2
                for l = 1:2
                    Stemp(k,l) = -(log(r)-eps*rho)*(k==l) + (p(k)-tempStks(n,k))*(p(l)-tempStks(n,l))*rho/r;
                end
            end
    
            U = Stemp*Ftemp;
            Uflowx(xposition) = Uflowx(xposition) + U(1);
            Uflowy(xposition) = Uflowy(xposition) + U(2);
    
        end
    
    end


end

