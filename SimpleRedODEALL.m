function dtxt = SimpleRedODEALL(t,x)
dtxt = zeros(32, 1); 

global S
global L
global G
global F


dtxt = [S, zeros(16,16); G'*L, F] * x;

end
