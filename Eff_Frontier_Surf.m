clear all;
clc;

load('Index20052016');

%% Set Variables

a = alpha; % row vector
S = Sigma;
% inverse of the covariance
ind = 3:10;
S = S(ind,ind); % covariance
a = transpose(a(ind)); % column vector
v = ones(size(a)); % row vector
I = eye(size(S)); % identity
[ExpStd,ExpCorr] = cov2corr(S);
iS{1} = inv(S);
iS{2} = I;
iS{3} = inv(corr2cov(ones(length(S),1)*mean(diag(S)),ExpCorr)); % average variances + corr
iS{4} = inv(corr2cov(ones(length(S),1)*mean(diag(S)),I)); % average variances + I

%% Use Covariance Matrix

i =1;
A = v'*iS{i}*v;
B = a'*iS{i}*v;
C = a'*iS{i}*a;
D = A*C-B^2;
% 1.1. Define the inline functions for the different optimisations
ta = iS{i}*a/B; % optimal risky
t0 = iS{i}*v/A; % global minimum variance
ra = ta' * a;
sa = ta' * S * ta;
st0 = t0' * S * t0;
rt0 = t0' * a;
%I II III IV V VII
x = [0:0.0025:0.5]; % risk
tI = (B/sqrt(C)) * x .* ta;
y = a'*tI; % returns
%VI
g0 = 1;
a0 = x;
omega = (B/D)*(a0*A-g0*B);
tVI = (g0 - omega) .* t0 + omega .* ta;
rVI = tVI'*a;
sVI = diag(tVI' * S * tVI);
ind2= find(sqrt(sVI)>max(x),1);
%VII
gamma = [50:-1:2];
mu = (1./gamma)*B;
g0 = 1;
tVII = (g0 - mu) .* t0 + mu .* ta;
rVII = tVII'*a;
sVII = diag(tVII' * S * tVII);
% 1.2 Plot for g0=1
figure;
plot(12*x,12*y,'Color','black')
ylabel('\alpha_p');
xlabel('\sigma_p');
hold on;
scatter(12*sqrt(sa),12*ra,[],'red','filled');
scatter(12*sqrt(st0),12*rt0,[],'blue','filled');
line(12*sqrt(sVI(1:ind2)),12*rVI(1:ind2),'Color','green');
%line(12*sqrt(sVII),12*rVII,'Color','red');
%hold off;
legend('I,II,III,IV,V','ORP (VIII g_0=1)','GMV (g_0=1)','VI,VII (g_0=1)');
xlim([0 1])
ylim([0 0.5])

%% Generate Data for Efficient Surfaces

g0 = 0.5:0.05:1.5;
% VI and VII
VIr = zeros(length(tVI),length(g0));
VIs = zeros(length(tVI),length(g0));
VIg = zeros(length(tVI),length(g0));
for j = 1:length(g0)
    tVI = (g0(j) - omega) .* t0 + omega .* ta;
    VIr(:,j) = tVI'*a;
    VIs(:,j) = diag(tVI' * S * tVI);
    VIg(:,j) = ones(length(tVI),1)*g0(j);
end
figure;
s1=surf(12*sqrt(VIs),VIg,12*VIr,'FaceColor','g');
%s1.AlphaData = gradient(VIr);
s1.FaceAlpha = 0.5; % 'flat';
zlabel('\alpha_p');
xlabel('\sigma_p');
ylabel('g_0');
hold on;
as1 = gca;
xlim([0 1]) % risk
ylim([0.5 1.5]) % gearing
zlim([0 0.5]) % alpha
% I,II,III,VI,V
Vr = zeros(length(tVI),length(g0));
Vs = zeros(length(tVI),length(g0));
Vg = zeros(length(tVI),length(g0));
for j = 1:length(g0)
    tV = g0(j).* (B/sqrt(C)) * x .* ta;
    Vr(:,j) = tV'*a;
    Vs(:,j) = diag(tV' * S * tV);
    Vg(:,j) = ones(length(tVI),1)*g0(j);
end
s2=surf(as1,12*sqrt(Vs),Vg,12*Vr,'FaceColor','k');
%s2.AlphaData = gradient(Vr);
s2.FaceAlpha = 0.1; %'flat';
hold on;
% plot OR portfolio line: VIII
ta2 = iS{1}*a/B; % optimal risky
t02 = iS{1}*v/A; % global minimum variance
ra2 = ta2' * a;
sa2 = ta2' * S * ta2; % use actual covariance
st02 = t02' * S * t02; % use actual covariance
rt02 = t02' * a;
for j = 1:length(g0)
    r3(:,j) = g0(j)*[ra2,rt02];
    s3(:,j) = g0(j)^2*[sa2;st02];
    g3(:,j) = g0(j)*[1,1];
end
% plot OR portfolio line
plot3(as1,12*sqrt(s3(1,:)),g3(1,:),12*r3(1,:),'-','Color','r','LineWidth',2);
% plot GMV portfolio line
plot3(as1,12*sqrt(s3(2,:)),g3(2,:),12*r3(2,:),'-','Color','b','LineWidth',2);
% legend('VI,VII','V,VIII','ORP','GMV');
% set camera view
view([-60,16]);
