T1 = [3 8 7; 5 2 7; 1 5 6];
T2 = [9 4 5; 7 9 6; 7 6 3];

u1 = [7 5];
u2 = [3 8];
Tclass = ones(3,3);

T = T2; % run for both classes
u = u2;
for i = 1:3
    T_c1 = abs(T-u(1)).^2;
    T_c2 = abs(T-u(2)).^2;
    for j = 1:3
        for k = 1:3
            if T_c1(j,k) < T_c2(j,k)
                Tclass(j,k) = 1;
            else
                Tclass(j,k) = 2;
            end
        end
    end
    u = [mean(T(find(Tclass ==1))) mean(T(find(Tclass ==2)))]
end