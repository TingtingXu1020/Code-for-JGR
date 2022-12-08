function Hill_P = Hill_P(a, N, M, Stiffness)
N_ELE = N*M;
[Alp, Bet, w_pq] = GaussGGLQ(N, M);
Hill_P = zeros(3, 3, 3, 3);

parfor n_ele = 1 : N_ELE
    oAux(:, :, :, :, n_ele) =  HillAux(Alp, Bet, w_pq, a, Stiffness, n_ele);
end
 
Aux = sum(oAux, 5);

for i = 1:3
    for j = 1:3
        for k = 1:3
            for l = 1:3
                Hill_P(i, j, k, l) = Hill_P(i, j, k, l) + 1/2 * (Aux(i, j, k, l) + Aux(i, j, l, k));
            end
        end
    end
end
end




