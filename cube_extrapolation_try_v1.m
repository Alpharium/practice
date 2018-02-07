for x=0:10
   a = sprintf('%d %d',x,x.^3+x.^2+x);
   disp(a)
   a = sprintf('%d %d',-x,(-x).^3+(-x).^2-x);
   disp(a)
end

