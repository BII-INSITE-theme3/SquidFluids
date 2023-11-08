function [Uxout,Uyout] = U_pois(xin,yin)

    load('outputs/U_800pts_Poissuelle.mat');

    a = find( min( abs(x-xin) ) );
    b = find( min( abs(y-yin) ) );

    % Some amount of interpolation?

    Uyout = UxPois(a,b);
    Uxout = UyPois(a,b);

end