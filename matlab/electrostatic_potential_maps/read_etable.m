function[etablemap] = read_etable(file)
    fid = fopen(file);
    ii = 1;
    k = {};
    v = {};
    
    while(~feof(fid))
        ln = fgetl(fid);
        data = split(ln);
        if( strcmp( data{1}, '#' ) ~= 1 )
            %atom name,  amplitude, sigma
            ki = data{1};
            vi = [str2double(data{2}), str2double(data{3})];
            k{end+1} = ki;
            v{end+1} = vi;
            ii = ii+1;
        end
    end
    etablemap = containers.Map(k,v);
    
end