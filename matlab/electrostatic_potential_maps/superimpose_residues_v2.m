
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

for(ii = 2:length(tyrpdbs))
  
    testala = pdbread(strcat('test/',tyrpdbs{ii}));

    nn = strfind(tyrpdbs{ii},'.pdb');
    tm = tyrpdbs{ii};
    PREFIX = tm(1:(nn-1));
   
    
    addpath('~/projects/electrostatic_potential_maps/src/rigid_transform_3D/')

    refcoords = [squeeze([refala.Model.Atom(:).X]); ...
                 squeeze([refala.Model.Atom(:).Y]); ...
                 squeeze([refala.Model.Atom(:).Z])];

    testcoords = [squeeze([testala.Model.Atom(:).X]); ...
                  squeeze([testala.Model.Atom(:).Y]); ...
                  squeeze([testala.Model.Atom(:).Z])];

    %transforms testcoords to match refcoords
    [R,t]=rigid_transform_3D(testcoords, refcoords)

    xformcoords = R*testcoords + t;

%read in maps
    addpath('~/src/Matlab_EM/MRCIO/')
    addpath('~/src/Matlab_EM/EMIODist/')
    refPREFIX
    PREFIX
    refala_map = ReadMRC(strcat('test/',refPREFIX,'.mrc'));
    testala_map = ReadMRC(strcat('test/',PREFIX,'.mrc'));


    %eul = rad2deg(rotm2eul(R',"ZYX"));
    %only works for same position, different chain
    %ZYX convention
    %tformimg = imrotate3(testala_map,-(eul(3)),[0,1,0],'cubic','crop');
    %tformimg = imrotate3(tformimg,-(eul(2)),[0,0,1],'cubic','crop');
    %tformimg = imrotate3(tformimg,-(eul(1)),[1,0,0],'cubic','crop');
    
    eul = rad2deg(rotm2eul(R',"XYZ"));
    
    %ZYZ convention
    tformimg = imrotate3(testala_map,-(eul(3)),[0,0,1],'cubic','crop');
    tformimg = imrotate3(tformimg,-(eul(2)),[0,1,0],'cubic','crop');
    tformimg = imrotate3(tformimg,-(eul(1)),[1,0,0],'cubic','crop');
    
    
    %????
    %tformimg = imtranslate(tformimg,t'); 
      

subplot(3,3,1)
showImage(squeeze(sum(refala_map,1)))
subplot(3,3,2)
showImage(squeeze(sum(testala_map,1)))
subplot(3,3,3)
showImage(squeeze(sum(tformimg,1)))

subplot(3,3,4)
showImage(squeeze(sum(refala_map,2)))
subplot(3,3,5)
showImage(squeeze(sum(testala_map,2)))
subplot(3,3,6)
showImage(squeeze(sum(tformimg,2)))

subplot(3,3,7)
showImage(squeeze(sum(refala_map,3)))
subplot(3,3,8)
showImage(squeeze(sum(testala_map,3)))
subplot(3,3,9)
showImage(squeeze(sum(tformimg,3)))


WriteMRC(tformimg,1,strcat(output_dir,'test-aln.mrc'));
end