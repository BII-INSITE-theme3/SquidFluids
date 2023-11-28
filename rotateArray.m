function [array] = rotateArray(array,shift)
    n = length(array);
    shift = mod(shift,n);
    temp = zeros(n,1);
    m = n-shift;
    temp(1+shift:end) = array(1:m);
    temp(1:shift) = array(m+1:end);
    array = temp;
end

