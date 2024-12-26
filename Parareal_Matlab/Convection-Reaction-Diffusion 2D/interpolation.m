function [uF] = interpolation(xC,yC,uC,xF,yF)

 [xCinterp,yCinterp] = meshgrid(xC,yC);
 [xFinterp,yFinterp] = meshgrid(xF,yF);
 uF = interp2(xCinterp,yCinterp,uC,xFinterp,yFinterp);