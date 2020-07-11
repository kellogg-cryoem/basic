function[] = write_bild_vector(output_file,start,dir)
    e = start + dir;
    fid = fopen(output_file,'w');
    fprintf(fid,'.arrow %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f',start(1),start(2),start(3), ...
        e(1),e(2),e(3));
    fclose(fid)
end