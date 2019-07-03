function [dum_interpolated] =LinearInterp2D(xp,zp,xe1,ze1,xe2,ze2,xe3,ze3,...
          dummy1,dummy2,dummy3)
a= ((ze2-ze3)*(xp-xe3)+(xe3-xe2)*(zp-ze3))/((ze2-ze3)*(xe1-xe3)+(xe3-xe2)*(ze1-ze3));
b= ((ze3-ze1)*(xp-xe3)+(xe1-xe3)*(zp-ze3))/((ze2-ze3)*(xe1-xe3)+(xe3-xe2)*(ze1-ze3));
c= 1-a-b; 
dum_interpolated=a*dummy1+b*dummy2+c*dummy3;    
end