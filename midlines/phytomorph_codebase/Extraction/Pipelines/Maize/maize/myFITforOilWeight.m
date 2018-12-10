function [] = myFITforOilWeight(p,W,T,P)
   p = poly(p);
   p1 = polyder(p,1);
   pv = polyval(p1,linspace(0,10^5));
   
end