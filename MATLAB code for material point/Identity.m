function [fourthK, fourthI, Ja, fourthJ] = Identity()
% FourIdentity.m
% 4th-order identity tensors
%--------------------------------------------------------------------------

%  allocate J, fourthJ
   J  = zeros(3, 3, 3, 3);
   J1 = zeros(3, 3, 3, 3);
   fourthJ = zeros(3, 3, 3, 3);

%  Eqn(3) in Jiang(2014)
   for i = 1:3
       for j = 1:3
           for k = 1:3
               for l = 1:3
                   J(i, j, k, l)  = KronD(i, k)*KronD(j, l);
                   J1(i, j, k, l) = KronD(j, k)*KronD(i, l);
                   fourthJ(i, j, k, l) = KronD(i, j)*KronD(k, l)/3;
               end
           end
       end
   end
   fourthI = 0.5*(J + J1);
   Ja = 0.5*(J - J1);
   fourthK = fourthI - fourthJ;   
end