% get the inverse of 4th order tensor F (with minor symmetry)


function Finverse = inversegeneral(F)
%% first method
b(:,:,1) = 1/sqrt(6)*[-1,0,0;0,-1,0;0,0,2];
b(:,:,2) = 1/sqrt(2)*[-1,0,0;0,1,0;0,0,0];
b(:,:,3) = 1/sqrt(2)*[0,0,0,;0,0,1;0,1,0];
b(:,:,4) = 1/sqrt(2)*[0,0,1,;0,0,0;1,0,0];
b(:,:,5) = 1/sqrt(2)*[0,1,0,;1,0,0;0,0,0];
b(:,:,6) = 1/sqrt(3)*[1,0,0,;0,1,0;0,0,1];
C = zeros(6,6);
for lambda = 1:6
    for eta = 1:6
        for i = 1:3
            for j =1:3
                for k=1:3
                    for l =1:3
                        C(lambda,eta) = C(lambda,eta) + F(i,j,k,l)*b(i,j,lambda)*b(k,l,eta);
                    end
                end
            end
        end
    end
end
if rcond(C) < 1e-15
    C_inv = pinv(C);
else
C_inv = inv(C);
end
F_inv = zeros(3,3,3,3);
for lambda = 1:6
    for eta = 1:6
        for i = 1:3
            for j =1:3
                for k=1:3
                    for l =1:3
                        F_inv(i,j,k,l) = F_inv(i,j,k,l) + C_inv(lambda,eta)*b(i,j,lambda)*b(k,l,eta);
                    end
                end
            end
        end
    end
end
Finverse = F_inv;

%% third method (Mandel notation)
% F_matrix = tensor2matrix(F);
% if rcond(F_matrix) < 1e-15
%     W = pinv(F_matrix);
% else
% W = inv(F_matrix);
% end
% Finverse = matrix2tensor(W);
end




