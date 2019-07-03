function [dummy] = Integrate(Grid,fieldvector)
[n1,n2]=size(fieldvector);
x=Grid.coord(:,1);
y=Grid.coord(:,2);   
dummy=0;
for k=1:Grid.nelem
    Temp=Grid.LEM(:,k)';
    a=Temp(1);
    b=Temp(2);
    c=Temp(3);
    xe1=x(a);
    xe2=x(b);
    xe3=x(c);
    ye1=y(a);
    ye2=y(b);
    ye3=y(c);
    xp=(xe1+xe2+xe3)/3;
    yp=(ye1+ye1+ye1)/3;
    dummy1=fieldvector(a);
    dummy2=fieldvector(b);
    dummy3=fieldvector(c);
    [dum_interpolated] =LinearInterp2D(xp,yp,xe1,ye1,xe2,ye2,xe3,ye3,...
          dummy1,dummy2,dummy3);
    dummy=dummy+dum_interpolated*Grid.CellVolumes(k);
   end

