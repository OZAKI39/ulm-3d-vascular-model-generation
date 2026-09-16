function [dvx]=min_velx(y1,y2,vx1,vx2,ly,delv)
%minimum image for velocity x

dy = y2-y1;

cory = round(dy / ly); % assume shear rate = dvx/dy
dvx = vx2 -(vx1 + cory * delv);

end
