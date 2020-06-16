
fid = fopen('test/tyr.lst')
tyrpdbs = {};
ii = 1;
while(~feof(fid))
    tyrpdbs{ii} = fgetl(fid);
    ii = ii+1;
end

output_dir = 'test/TYR/';
refala = pdbread(strcat('test/',tyrpdbs{1}));
tm = tyrpdbs{1};
nn = strfind(tyrpdbs{1},'.pdb');
refPREFIX = tm(1:(nn-1));

for(ii = 1:length(tyrpdbs))
  
    testala = pdbread(strcat('test/',tyrpdbs{ii}));

    %[d,rmsd,trans] = pdbsuperpose(refala,testala);
    nn = strfind(tyrpdbs{ii},'.pdb');
    tm = tyrpdbs{ii};
    PREFIX = tm(1:(nn-1));
   
    
    addpath('~/projects/electrostatic_potential_maps/src/rigid_transform_3D/')

    refcoords = [squeeze([refala.Model.Atom(:).X]); ...
                 squeeze([refala.Model.Atom(:).Y]); ...
                 squeeze([refala.Model.Atom(:).Z])];
         
    size(refcoords)

    testcoords = [squeeze([testala.Model.Atom(:).X]); ...
                  squeeze([testala.Model.Atom(:).Y]); ...
                  squeeze([testala.Model.Atom(:).Z])];
    size(testcoords)
    %transforms testcoords to match refcoords
    [R,t]=rigid_transform_3D(testcoords, refcoords)

    xformcoords = R*testcoords + t;

    %replace coordinates
    for(i = 1:length(xformcoords))
        testala.Model.Atom(i).X = xformcoords(1,i);
        testala.Model.Atom(i).Y = xformcoords(2,i);
        testala.Model.Atom(i).Z = xformcoords(3,i);
    end
    pdbwrite(strcat(output_dir,PREFIX,'-aln.pdb'),testala);

%read in maps
    addpath('~/src/Matlab_EM/MRCIO/')
    addpath('~/src/Matlab_EM/EMIODist/')
    refala_map = ReadMRC(strcat('test/',refPREFIX,'.mrc'));
    testala_map = ReadMRC(strcat('test/',PREFIX,'.mrc'));

    %change basis for pdbs, then compute rotation
    %newpdb1=Rpdb2mrc*refcoords+tpdb2mrc;
    %newpdb2=Rpdb2mrc*testcoords+tpdb2mrc;
    %[R,t]=rigid_transform_3D(newpdb2, newpdb1);
        
    
    
    eul = rad2deg(rotm2eul(R,"ZYZ"));
    
    %ZYZ convention
    tformimg = imrotate3(imrotate3(imrotate3(testala_map,-(eul(3)),[0,0,1],'cubic','crop') ...
                                                        ,-(eul(2)),[0,1,0],'cubic','crop') ...
                                                        ,-(eul(1)),[0,0,1],'cubic','crop');
    %tformimg = imtranslate(tformimg,t'); 
      

%cropped_maps
    ll = [centers-bsize/2,...
          centers+bsize/2];
    cropped_ref = refala_map;
    cropped_tform = tformimg;
%cropped_ref = refala_map(ll(1,1):ll(1,2),ll(2,1):ll(2,2),ll(3,1):ll(3,2));
%cropped_tform = tformimg(ll(1,1):ll(1,2),ll(2,1):ll(2,2),ll(3,1):ll(3,2));

subplot(2,2,1)
showImage(squeeze(sum(cropped_ref,3)))
subplot(2,2,2)
showImage(squeeze(sum(testala_map,3)))
subplot(2,2,3)
showImage(squeeze(sum(cropped_tform,3)))
   % if( ii > 1 )
   %     WriteMRC(cropped_tform,0.656,strcat(output_dir,PREFIX,'-aln.mrc'));
   % else
   %     WriteMRC(cropped_ref,0.656,strcat(output_dir,refPREFIX,'-aln.mrc'));
   % end
end