check = refala_map;

rotmap = imrotate3(check,45,[0, 0, 1],'cubic','crop');
rotmap = imrotate3(rotmap,10,[0, 1, 0],'cubic','crop');
rotmap = imrotate3(rotmap,30,[0, 0, 1],'cubic','crop');


Rz1 = [cosd(45) -sind(45) 0; sind(45) cosd(45) 0; 0 0 1];
Ry = [cosd(10) 0 sind(10); 0 1 0; -sind(10) 0 cosd(10)];
Rz2 = [cosd(30) -sind(30) 0; sind(30) cosd(30) 0; 0 0 1];

R = Rz2*Ry*Rz1;
 eul = rad2deg(rotm2eul(R,"ZYZ"));
eul

rotback = imrotate3(rotmap,-(eul(3)),[0,0,1],'cubic','crop');
rotback = imrotate3(rotback,-(eul(2)),[0,1,0],'cubic','crop');
rotback = imrotate3(rotback,-(eul(1)),[0,0,1],'cubic','crop');

clf
subplot(3,3,1)
showImage(squeeze(sum(check,3)));
subplot(3,3,2)
showImage(squeeze(sum(rotmap,3)));
subplot(3,3,3)
showImage(squeeze(sum(rotback,3)));

test2coords = R*refcoords
subplot(3,3,4)
scatter(test2coords(1,:),test2coords(2,:),'blue')
hold on
scatter(refcoords(1,:),refcoords(2,:),'red')

[R,t]=rigid_transform_3D(test2coords, refcoords)
suptestcoords = R*test2coords + t;
subplot(3,3,5)
scatter(test2coords(1,:),test2coords(2,:),'blue')
hold on
scatter(refcoords(1,:),refcoords(2,:),'red')
scatter(suptestcoords(1,:),suptestcoords(2,:),'green')


%not correct
 eul = rad2deg(rotm2eul(R',"ZYZ"));
 rotback2 = imrotate3(rotmap,-(eul(3)),[0,0,1],'cubic','crop');
rotback2 = imrotate3(rotback2,-(eul(2)),[0,1,0],'cubic','crop');
rotback2 = imrotate3(rotback2,-(eul(1)),[0,0,1],'cubic','crop');
subplot(3,3,6)
showImage(squeeze(sum(rotback2,3)))


%change basis for pdbs, then compute rotation
newpdb1=Rpdb2mrc*refcoords+tpdb2mrc;
newpdb2=Rpdb2mrc*test2coords+tpdb2mrc;
[R,t]=rigid_transform_3D(newpdb2, newpdb1);
eul = rad2deg(rotm2eul(R',"ZYZ"));
 rotback2 = imrotate3(rotmap,-(eul(3)),[0,0,1],'cubic','crop');
rotback2 = imrotate3(rotback2,-(eul(2)),[0,1,0],'cubic','crop');
rotback2 = imrotate3(rotback2,-(eul(1)),[0,0,1],'cubic','crop');
subplot(3,3,7)
showImage(squeeze(sum(rotback2,3)))
