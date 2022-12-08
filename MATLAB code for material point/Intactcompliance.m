function Intact = Intactcompliance(E, nu)
Intact = zeros(3, 3, 3, 3);
    for i = 1:3
        for j = 1:3
            for k = 1:3
                for l=1 : 3
                     Intact(i, j, k, l) = Intact(i, j, k, l)...
                         + (-1)*nu /E *KronD(i,j)*KronD(k,l) + (1 +...
                         nu)/2/E*(KronD(i,k)*KronD(j,l)+KronD(i,l)*KronD(j,k));
                end
            end
        end
    end
   
end