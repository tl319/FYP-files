function dtxt = SimpleRedODE2(t,x)
dtxt = zeros(32, 1); 

global S
global L
global G2
global F2


dtxt = [S, zeros(16,16); G2'*L, F2] * x;

end
