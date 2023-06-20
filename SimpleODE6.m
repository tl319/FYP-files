function dtxt = SimpleODE6(t,x)

global S
global L
global A_vary
global B_fom



dtxt = [S, zeros(16,1006); B_fom*L, A_vary(:,:,6)] * x;

end
