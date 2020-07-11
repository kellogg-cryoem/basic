function[p] = plane_normal(pts)
    %pts are row vectors
    pts = pts - mean(pts);
    [U,S,W]=svd(pts,0);
    p = W(:,end);
end