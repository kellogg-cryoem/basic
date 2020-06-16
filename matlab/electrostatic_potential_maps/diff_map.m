function[diffm] = diff_map(map2,map1)
    %map1 is the reference 
    %map2matchd = imhistmatchn(map2,map1);
    %diffm = map1-map2matchd;
    diffm = map2-map1;
end