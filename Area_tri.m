function [Grid] = Area_tri(Grid)
for i=1:Grid.nelem
    Temp=Grid.LEM(:,i)';
    a=Temp(1);
    b=Temp(2);
    c=Temp(3);
    Ax=Grid.coord(a,1);
    Ay=Grid.coord(a,2);
    Bx=Grid.coord(b,1);
    By=Grid.coord(b,2);
    Cx=Grid.coord(c,1);
    Cy=Grid.coord(c,2);
    area_temp=(Ax*(By-Cy)+Bx*(Cy-Ay)+Cx*(Ay-By))/2;
    Grid.CellVolumes(i)=abs(area_temp);
end    

end
