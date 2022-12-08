function T=Transform(X ,Q)
% Transform.m
%  tensors transformation between coordinate systems.
%
% Input:  second-order or fourth-order tensors
%--------------------------------------------------------------------------

if ndims(X) == 4
    T = zeros(3,3,3,3);
    for i=1:3
        for j=1:3
            for k=1:3
                for l=1:3
                    for m=1:3
                        for n=1:3
                            for o=1:3
                                for p=1:3
                                    T(m,n,o,p)= T(m,n,o,p) + Q(m,i)*Q(n,j)*Q(o,k)*Q(p,l)*X(i,j,k,l);
                                end
                            end
                        end
                    end
                end
            end
        end
    end
else
    T = zeros(3,3);
    for i=1:3
        for j=1:3
            for m=1:3
                for n=1:3
                    T(m,n)= T(m,n) + Q(m,i)*Q(n,j)*X(i,j);
                end
            end
        end
    end
end
end