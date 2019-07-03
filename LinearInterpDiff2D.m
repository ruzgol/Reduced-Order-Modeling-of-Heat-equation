function [dummyDifx,dummyDifz]= LinearInterpDiff2D(xp,yp,xe1,ye1,xe2,ye2,xe3,...
          ye3,dummy1,dummy2,dummy3)
ax= ((ye2-ye3))/((ye2-ye3)*(xe1-xe3)+(xe3-xe2)*(ye1-ye3));
bx= ((ye3-ye1))/((ye2-ye3)*(xe1-xe3)+(xe3-xe2)*(ye1-ye3));
cx=-((ye2-ye1))/((ye2-ye3)*(xe1-xe3)+(xe3-xe2)*(ye1-ye3));
az= ((xe3-xe2))/((ye2-ye3)*(xe1-xe3)+(xe3-xe2)*(ye1-ye3));
bz= ((xe1-xe3))/((ye2-ye3)*(xe1-xe3)+(xe3-xe2)*(ye1-ye3));
cz=-((xe1-xe2))/((ye2-ye3)*(xe1-xe3)+(xe3-xe2)*(ye1-ye3));
dummyDifx=ax*dummy1+bx*dummy2+cx*dummy3;
dummyDifz=az*dummy1+bz*dummy2+cz*dummy3;
end           