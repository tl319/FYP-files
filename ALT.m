%clear 
close all

load("fom.mat","A","B","C");

global A_fom
A_fom = full(A);


global B_fom
%B_fom = ones(1006,1);
B_fom=B;
global C_fom;
%C_fom = ones(1,1006);
C_fom=C;

ft = linspace(0,10,11);


f= (10+10*(ft));

sys = ss(A_fom,B_fom,C_fom,0);
bode(sys);

figure
plot(f)


global A_vary
for idx = 1:length(f)
    A_vary(:,:,idx) = A_fom;
    A_vary(5,6) = 400-f(idx);
    A_vary(6,5) = -400+f(idx);
end    

%Select S, ensuring no shared eigen values with A
global S

S = zeros(16,16);
S(1,1) = 0; S(1,2) = 1; S(2,1) = -1; S(2,2) = 0;
S(3,3) = 0; S(3,4) = 50; S(4,3) = -50; S(4,4) = 0;
S(5,5) = 0; S(5,6) = 100; S(6,5) = -100; S(6,6) = 0;
S(7,7) = 0; S(7,8) = 200; S(8,7) = -200; S(8,8) = 0;
S(9,9) = 0; S(9,10) = 400; S(10,9) = -400; S(10,10) = 0;
S(11,11) = 0; S(11,12) = 75; S(12,11) = -75; S(12,12) = 0;
S(13,13) = 0; S(13,14) = 500; S(14,13) = -500; S(14,12) = 0;
S(15,15) = 0; S(15,16) = 700; S(16,15) = -700; S(16,12) = 0;

%Select L, ensuring (S,L) Observable
global L
L = zeros(1,16);
L(1,1) = 1; L(1,2) = 0; L(1,3) = 1; L(1,4) = 0;
L(1,5) = 1; L(1,6) = 0; L(1,7) = 1; L(1,8) = 0;
L(1,9) = 1; L(1,10) = 0; L(1,11) = 1; L(1,12) = 0;
L(1,13) = 1; L(1,14) = 0; L(1,15) = 1; L(1,16) = 0;

%Scale L
L = L;

global PI
PI_Vary1= sylvester(A_vary(:,:,5),-S,-B_fom*L);
PI_Vary6= sylvester(A_vary(:,:,6),-S,-B_fom*L);
PI_Vary11= sylvester(A_vary(:,:,7),-S,-B_fom*L);



%Solve ODE, plug in ft and f for parameter
% x0_temp = rand(1022,1)'; 
% x0 = [x0_temp x0_temp x0_temp]';
tspan = [0 10];
[t,x1] = ode45(@SimpleODEALL,tspan,x0);
xT1 =x1';


%Calculate transient response
y1 = C_fom * xT1(17:1022,:);
yss1 = C_fom * PI_Vary1 * xT1(1:16,:);
y6 = C_fom * xT1(1039:2044,:);
yss6 = C_fom * PI_Vary6 * xT1(1023:1038,:);
y11 = C_fom * xT1(2061:3066,:);
yss11 = C_fom * PI_Vary11 * xT1(2045:2060,:);


dataw1 = xT1(1:16, end-20000:end)';
datay1 = y1(:,end-20000:end)';
dataw6 = xT1(1023:1038, end-20000:end)';
datay6 = y6(:,end-20000:end)';
dataw11 = xT1(2045:2060, end-20000:end)';
datay11 = y11(:,end-20000:end)';

CPI1 = dataw1 \ datay1;
CPI6 = dataw6 \ datay6;
CPI11 = dataw11 \ datay11;

error1 = CPI1 - (C_fom*PI_Vary1)';
error6 = CPI6 - (C_fom*PI_Vary6)';
error11 = CPI11 - (C_fom*PI_Vary11)';

sum1 = sum(error1);
sum2 = sum(error6);
sum3 = sum(error11);

%Create H for each second
global H_vary1;
H_vary1 = C_fom*PI_Vary1;
global H_vary6;
H_vary6 = C_fom*PI_Vary6;
global H_vary11;
H_vary11 = C_fom*PI_Vary11;

global eig_f
eig_f = [ -1+100*i, -1-100*i, -20+200*i, -20-200*i, -40+400*i, -40-400*i, -1, -100, -250, -350, -450, -550, -650, -750, -910, -1000,];


global G
G = place(S',L',eig_f);

global F
F = S- G'*L;

x0_red_temp = x0(1:32);
x0_red = [x0_red_temp];
tspan = t;
[t,x_generic_red] = ode45(@SimpleRedODE,tspan,x0_red);
xT1_red = x_generic_red';

timeline = 5:1:7;
temp1 = [5, 7]';
temp2 = [CPI1, CPI11]';
InterpCPI = interp1(temp1,temp2,timeline);
interpCPI6 = InterpCPI(2,:);
 
y_red1 = H_vary1 * xT1_red(17:32,:);
yss_red1 = H_vary1 * xT1_red(1:16,:);

y_red1_est = CPI1' * xT1_red(17:32,:);
yss_red1_est = CPI1' * xT1_red(1:16,:);

y_red6 = H_vary6 * xT1_red(17:32,:);
yss_red6 = H_vary6 * xT1_red(1:16,:);

y_red6_est = CPI6' * xT1_red(17:32,:);
yss_red6_est = CPI6' * xT1_red(1:16,:);

interp_y_red6 = interpCPI6 * xT1_red(17:32,:);
interp_yss_red6 = interpCPI6 * xT1_red(1:16,:);


y_red11 = H_vary11 * xT1_red(17:32,:);
yss_red11 = H_vary11 * xT1_red(1:16,:);

y_red11_est = CPI11' * xT1_red(17:32,:);
yss_red11_est = CPI11' * xT1_red(1:16,:);


figure
plot(t,y1);
hold on
plot(t,y_red1_est);
title('Comparing output of Original and Reduced model for p = 20')
legend('Y','Y Reduced')

figure
plot(t,y6);
hold on
plot(t,interp_y_red6)
title('Comparing output of Original and Reduced model for p = 70')
legend('Y','Y Reduced')

figure
plot(t,y11);
hold on
plot(t,y_red11_est);
title('Comparing output of Original and Reduced model for p = 120')
legend('Y','Y Reduced')

errorfinal = interpCPI6 - (C_fom*PI_Vary6);
sumf = sum(errorfinal);
