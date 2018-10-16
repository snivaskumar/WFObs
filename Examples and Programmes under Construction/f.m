function z = f(w,Pf,xf,tur);
   z = zeros(size(Pf{1}));
   for i = 1:tur
        z = z + w(i)*Pf{i};
   end
   z = trace(inv(z));
end