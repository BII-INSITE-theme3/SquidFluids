
% First, run code d_..

% Exclude the parts of the fluid domain that're invalid

cent1 = [-0.3,0];
cent2 =  [0.3,0];

ind = ones(length(x));

for ii = 1:length(x)
    for jj = 1:length(y)
        r1 = norm(cent1-[x(ii),y(jj)]);
        r2 = norm(cent2-[x(ii),y(jj)]);
        if r1 <= R01
            ind(jj,ii) = 0;
        end
        if r2 <= R02
            ind(jj,ii) = 0;
        end
    end
end

UmagTemp = Umag.*ind;
UxTemp = Ux.*ind;
UyTemp = Uy.*ind;

stagZone = zeros(length(x));
thresh = 0.33;

for ii = 1:length(x)
    for jj = 1:length(y)
        
        if abs(UmagTemp(ii,jj)) > 0 && abs(UmagTemp(ii,jj)) < thresh
            stagZone(ii,jj) = 1;
        end

    end
end

sum(sum(stagZone))

imagesc(stagZone)