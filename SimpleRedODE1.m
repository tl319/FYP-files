function dtxt = SimpleRedODE1(t,x)
dtxt = zeros(32, 1); 

global S
global L
global G1
global F1


dtxt = [S, zeros(16,16); G1'*L, F1] * x;

end
