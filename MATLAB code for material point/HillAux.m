function HillAux = HillAux(Alp, Bet, w_pq, a, Stiffness, n_ele)

    omega = Alp(n_ele);
    Zeta(3) = Bet(n_ele);
    Zeta(1) =(1-Zeta(3)^2)^(1/2) * cos(omega);
    Zeta(2) =(1-Zeta(3)^2)^(1/2) * sin(omega);
    
    Xi = zeros(3, 1);
    for i = 1:3
        Xi(i) = Zeta(i)/a(i);
    end
    
    K = zeros(3, 3);
    for i = 1:3
        for j = 1:3
            for k = 1:3
                for l = 1:3
                    K (i, k) = K (i, k) + Stiffness (i, j, k, l) * Xi(j) * Xi(l);
                end
            end
        end
    end
    
    D = 0;
    for m = 1:3
        for n = 1:3
            for l = 1:3
                D =  D + LeviCivita([m, n, l]) * K(m, 1) * K(n, 2) * K(l, 3);
            end
        end
    end
    
    N_bar = zeros(3, 3);
    for i = 1:3
        for j = 1:3
            for k = 1:3
                for l = 1:3
                    for m = 1:3
                        for n =1:3
                            N_bar(i, j) = N_bar(i, j) + 1/2 * LeviCivita([i, k, l]) * LeviCivita([j, m, n]) * K(k, m) * K(l, n);
                        end
                    end
                end
            end
        end
    end
    
    Green_f = zeros(3, 3, 3, 3);
    for i = 1:3
        for j = 1:3
            for k = 1:3
                for l = 1:3
                    Green_f(i, j, k, l) = Xi(k) * Xi(l) *N_bar(i, j)/D;
                end
            end
        end
    end
    
    Result= zeros(3, 3, 3, 3);
    for i = 1:3
        for j = 1:3
            for k = 1:3
                for l = 1:3                  
                    Result(i, j, k, l) = Result(i, j, k, l) + (Green_f(i, l, j, k) + Green_f(j, l, i, k));    
                end
            end
        end
    end
    
    HillAux = 1/(8*pi) * (1/2) * (2*pi-0) * (Result* w_pq(n_ele));
end