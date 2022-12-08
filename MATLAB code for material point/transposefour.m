function [outA]=transposefour(inA)
outA(1:3,1:3,1:3,1:3) = 0;
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                outA(i,j,k,l)= inA(k,l,i,j);
            end
        end
    end
end
end