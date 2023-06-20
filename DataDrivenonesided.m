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

for idx = 1:length(f)
    PI_Vary(:,:,idx) = sylvester(A_vary(:,:,idx),-S,-B_fom*L);
end


%Solve ODE, plug in ft and f for parameter
% x0 = rand(1022,1); 
% tspan = [0 10];
% [t,x1] = ode45(@SimpleODE1,tspan,x0);
% [t,x2] = ode45(@SimpleODE6,tspan,x0);
% [t,x3] = ode45(@SimpleODE11,tspan,x0);

xT1 =x1';
xT2 =x2';
xT3 =x3';

%Calculate transient response
y1 = C_fom * xT1(17:1022,:);
y2 = C_fom * xT2(17:1022,:);
y3 = C_fom * xT3(17:1022,:);
yss1 = C_fom * PI_Vary(:,:,1) * xT1(1:16,:);
yss2 = C_fom * PI_Vary(:,:,6) * xT2(1:16,:);
yss3 = C_fom * PI_Vary(:,:,11) * xT3(1:16,:);

dataw1 = xT1(1:16, end-20000:end)';
dataw2 = xT2(1:16, end-20000:end)';
dataw3 = xT3(1:16, end-20000:end)';
datay1 = y1(:,end-20000:end)';
datay2 = y2(:,end-20000:end)';
datay3 = y3(:,end-20000:end)';

CPI1 = dataw1 \ datay1;
CPI2 = dataw2 \ datay2;
CPI3 = dataw3 \ datay3;

error1 = CPI1 - (C_fom*PI_Vary(:,:,1))';
error2 = CPI2 - (C_fom*PI_Vary(:,:,6))';
error3 = CPI3 - (C_fom*PI_Vary(:,:,11))';

%Create H for each second
global H_vary;
for idx = 1:11;
     H_vary(:,idx)= C_fom*PI_Vary(:,:,idx);
end

global eig_f
% %alt_eig_f_real = [ -1+100*i, -1-100*i, -1+200*i, -1-200*i, -1+400*i, -1-400*i, -33, -150, -250, -350, -450, -550, -650, -750, -910, -1000,];
eig_f = [ -1+100*i, -1-100*i, -20+200*i, -20-200*i, -40+400*i, -40-400*i, -1, -100, -250, -350, -450, -550, -650, -750, -910, -1000,];
% %alt_eig_f_real = [ -1+100*i, -1-100*i, -1+200*i, -1-200*i, -1+400*i, -1-400*i, -33, -150, -250, -350, -450, -550, -650, -750, -999, -1000,];
% %alt_eig_f_real = [ 0, -40, -80, -150, -200, -230, -300, -340, -410, -480, -520, -590, -650, -750, -999, -1000,];
% 

global G
G = place(S',L',eig_f);

global F
F = S- G'*L;

global eig_f1
eig_f1 = [ -1+100*i, -1-100*i, -20+200*i, -20-200*i, -40+390*i, -40-390*i, 0, -100, -250, -350, -450, -550, -650, -750, -910, -1000,];
global G1
G1 = place(S',L',eig_f1);
global F1
F1 = S- G1'*L;

global eig_f2
eig_f2 = [ -1+100*i, -1-100*i, -20+200*i, -20-200*i, -40+340*i, -40-340*i, 0, -100, -250, -350, -450, -550, -650, -750, -910, -1000,];
global G2
G2 = place(S',L',eig_f2);
global F2
F2 = S- G2'*L;

global eig_f3
eig_f3 = [ -1+100*i, -1-100*i, -20+200*i, -20-200*i, -40+290*i, -40-290*i, 0, -100, -250, -350, -450, -550, -650, -750, -910, -1000,];
global G3
G3 = place(S',L',eig_f3);
global F3
F3 = S- G3'*L;

x0_red = x0(1:32);
tspan = t;
% tspan = linspace(0,100,35693);

[t,x_generic_red] = ode45(@SimpleRedODE,tspan,x0_red);
[t,x1_red] = ode45(@SimpleRedODE1,tspan,x0_red);
[t,x2_red] = ode45(@SimpleRedODE2,tspan,x0_red);
[t,x3_red] = ode45(@SimpleRedODE3,tspan,x0_red);

xTgen_red = x_generic_red';
xT1_red = x1_red';
xT2_red = x2_red';
xT3_red = x3_red';

ygen_red1 = H_vary(:,1)' * xTgen_red(17:32,:);
ygen_red2 = H_vary(:,6)' * xTgen_red(17:32,:);
ygen_red3 = H_vary(:,11)' * xTgen_red(17:32,:);
yssgen_red1 = H_vary(:,1)' * xTgen_red(1:16,:);
yssgen_red2 = H_vary(:,6)' * xTgen_red(1:16,:);
yssgen_red3 = H_vary(:,11)' * xTgen_red(1:16,:);

ygen_red1_est = CPI1' * xTgen_red(17:32,:);
ygen_red2_est = CPI2' * xTgen_red(17:32,:);
ygen_red3_est = CPI3' * xTgen_red(17:32,:);
yssgen_red1_est = CPI1' * xTgen_red(1:16,:);
yssgen_red2_est = CPI2' * xTgen_red(1:16,:);
yssgen_red3_est = CPI3' * xTgen_red(1:16,:);
 
y_red1 = H_vary(:,1)' * xT1_red(17:32,:);
y_red2 = H_vary(:,6)' * xT2_red(17:32,:);
y_red3 = H_vary(:,11)' * xT3_red(17:32,:);
yss_red1 = H_vary(:,1)' * xT1_red(1:16,:);
yss_red2 = H_vary(:,6)' * xT2_red(1:16,:);
yss_red3 = H_vary(:,11)' * xT3_red(1:16,:);

y_red1_est = CPI1' * xT1_red(17:32,:);
y_red2_est = CPI2' * xT2_red(17:32,:);
y_red3_est = CPI3' * xT3_red(17:32,:);
yss_red1_est = CPI1' * xT1_red(1:16,:);
yss_red2_est = CPI2' * xT2_red(1:16,:);
yss_red3_est = CPI3' * xT3_red(1:16,:);

timeline = 0:1:10;
temp1 = [0, 10]';
temp2 = [CPI1, CPI3]';
InterpCPI = interp1(temp1,temp2,timeline);
interpCPI6 = InterpCPI(6,:);

errorfinal = interpCPI6' -  (C_fom*PI_Vary(:,:,6))';

estygen_red2 = interpCPI6 * xTgen_red(17:32,:);
estyssgen_red2 = interpCPI6 * xTgen_red(1:16,:);
estygen_red2_est = interpCPI6 * xTgen_red(17:32,:);
estyssgen_red2_est = interpCPI6 * xTgen_red(1:16,:);
esty_red2 = interpCPI6 * xT2_red(17:32,:);
estyss_red2 = interpCPI6 * xT2_red(1:16,:);
esty_red2_est = interpCPI6 * xT2_red(17:32,:);
estyss_red2_est = interpCPI6 * xT2_red(1:16,:);

figure
plot(t,y3);
hold on
plot(t,y_red3);

figure
plot(t,y3);
hold on
plot(t,y_red3_est)


figure
plot(t,y2);
hold on
plot(t,y_red2);

figure
plot(t,y2);
hold on
plot(t,estygen_red2);

figure
plot(t,estygen_red2)
hold on
plot(t,estyssgen_red2)
% hold on
% plot(t,y_red1_est)
% plot(t,ygen_red1_est);
% hold on
% title('Comparison of system response with estimated PI for p = 0')
% legend('y','y red with correct PI','y red with estimated PI')

