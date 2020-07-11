function[r] = corr(A,B)
    A = A - mean(A,'all');
    B = B - mean(B,'all');
    
    num = sum( A .* B, 'all');
    den1 = sum(A.^2, 'all');
    den2 = sum(B.^2, 'all');
    r = num ./ (den1 .* den2)^0.5;
end