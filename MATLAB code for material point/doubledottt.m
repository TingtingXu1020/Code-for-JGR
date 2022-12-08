%A is 2th order tensor and B is 2nd order tensor
%result is a scalar

function R = doubledottt(A, B)
R = 0;
for i = 1:3
    for j = 1:3
                R = R + A(i, j)*B(i, j);
    end
end
end


            
        
