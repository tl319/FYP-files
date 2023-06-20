function dtxt = SimpleRedODE3(t,x)
dtxt = zeros(32, 1); 

global S
global L
global G3
global F3


dtxt = [S, zeros(16,16); G3'*L, F3] * x;

end
