function dtxt = SimpleODEALL(t,x)

global S
global L
global A_vary
global B_fom


dtxt = [S, zeros(16,1006) zeros(16,2044); B_fom*L, A_vary(:,:,5), zeros(1006,2044);
        zeros(16,1022), S, zeros(16,1006), zeros(16,1022); zeros(1006,1022), B_fom*L, A_vary(:,:,6), zeros(1006,1022);
        zeros(16,2044), S, zeros(16,1006); zeros(1006,2044), B_fom*L, A_vary(:,:,7)
        ] * x;

end
