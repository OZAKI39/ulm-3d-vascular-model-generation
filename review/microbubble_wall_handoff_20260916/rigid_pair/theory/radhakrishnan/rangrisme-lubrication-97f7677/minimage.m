function [rx,ry,rz,del]=minimage(lx,ly,lz,delx,px1,py1,pz1,px2,py2,pz2)
%minimum image convention for positions

rx = px2 - px1;
ry = py2 - py1;
rz = pz2 - pz1;

cory = round(ry / ly); % assume shear rate = dvx/dy
rx = rx - cory * delx;
rx = rx - round(rx/lx) * lx;
ry = ry - cory * ly;
rz = rz - round(rz/lz) * lz;



del = sqrt(rx.*rx + ry.*ry + rz.*rz);

end

