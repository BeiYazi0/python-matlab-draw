close all;clear all;clc;
%% 2
file = 'qsiwell5.csv';
L = importdata(file);


%% 4,5
imgname1 = 'Figure_1_part1.png';
imgname2 = 'Figure_1_part2.png';
plotlog(L, 2100, 2250, 0.3, 0.5, imgname1);
plotlog(L, 2100, 2500, 0.3, 0.5, imgname2);


%% 11
global RHO_o K_o RHO_g K_g K_sh K_qz MU_sh MU_qz RHO_sh RHO_qz K_b RHO_b
RHO_qz = 2.6;   K_qz = 37;  MU_qz = 44;
RHO_sh = 2.8;   K_sh = 15;  MU_sh = 5;
RHO_b = 1.1 ;   K_b = 2.8 ;
RHO_o = 0.8;    K_o = 0.9;
RHO_g = 0.2;    K_g = 0.06; 

Cn = 8;
phic = 0.4;
f = 1;

phi = linspace(0.01, 0.4, 50); % 不填第三个系数时，python默认为50，matlab默认为100


%% 12
K0 = K_qz;
MU0 = MU_qz;
RHO0 =  RHO_qz;
[Kdry, MUdry] = softsand(K0, MU0, phi, phic, Cn, 45);
[vp_ssm, vs_ssm, rho_ssm, ~] = vels(Kdry, MUdry, K0, RHO0, K_b, RHO_b, phi);

[Kdry, MUdry] = stiffsand(K0, MU0, phi, phic, Cn, 45);
[vp_sti, vs_sti, rho_sti, ~]= vels(Kdry, MUdry, K0, RHO0, K_b, RHO_b, phi);


%% 13
z1=2100;
z2=2250;
cutoff_sand = 0.3;
cutoff_shale = 0.5;

data = L.data;
textdata = L.textdata;

index = data(:,(textdata=="DEPTH"));
VSH_data = data(:,(textdata=="VSH"));
PHIE_data = data(:,(textdata=="PHIE"));
VP_data = data(:,(textdata=="VP"));
VS_data = data(:,(textdata=="VS"));
RHO_data = data(:,(textdata=="RHO"));

ss=(index>=z1) & (index<=z2) & (VSH_data<=cutoff_sand);
sh=(index>=z1) & (index<=z2) & (VSH_data>=cutoff_shale);

figure()
hold on
a = plot(PHIE_data(ss), VP_data(ss), '-go', 'markersize', 6);
plot(phi, vp_ssm, '-k');
plot(phi, vp_sti, ':k');
alpha(0.5);
hold off
xlim([0,0.4]);ylim([2e3,4e3]);
xlabel('porosity PHIE');ylabel('velocity VP [m/s]');
title('Soft Sand and Stiff Sand models');

% imgname = 'Figure_2_part1.png'; % '/Users/matt/Dropbox/Agile/SEG/Tutorials/delMonte_Apr2017/Figure_2_part1.png';
% print(gcf, '-dpng', '-r600', imgname);


%% 14
K_Voigt = 0.8*K_qz + 0.2*K_sh;
K_Reuss = 1./ ( 0.8./K_qz + 0.2./K_sh );
K_Hill = (K_Voigt + K_Reuss) / 2;

MU_Voigt = 0.8*MU_qz + 0.2*MU_sh;
MU_Reuss = 1./ ( 0.8./MU_qz + 0.2./MU_sh );
MU_Hill = (MU_Voigt + MU_Reuss) / 2;

fprintf('K_VRH = %.2f, mu_VRH = %.2f\n', K_Hill, MU_Hill);


%% 15
rho = 0.8*RHO_qz + 0.2*RHO_sh;

fprintf('RHO = %.2f\n',rho);


%% 16
NG = linspace(0.6, 1.0, 5);

figure()

labels = [];
for i = NG
    [~, ~, K0] = vrh(i,K_qz,K_sh);
    [~, ~, MU0] = vrh(i,MU_qz,MU_sh);    
    RHO0 = i*RHO_qz+(1-i)*RHO_sh;
    [Kdry, MUdry] = softsand(K0, MU0, phi, 0.5, 12, 45);
    [vp_ssm, vs_ssm, rho_ssm, ~]= vels(Kdry,MUdry,K0,RHO0,K_b,RHO_b,phi);
    [Kdry, MUdry] = stiffsand(K0, MU0, phi, 0.4, 8, 45);
    [vp_sti, vs_sti, rho_sti, ~]= vels(Kdry,MUdry,K0,RHO0,K_b,RHO_b,phi);
    label = ['N:G=', sprintf('%.2f',i)];
    labels=[labels,{label}];
    subplot(1,2,1);
    hold on
    a = plot(phi, vp_ssm, '-k');
    a.Color(4) = i-0.4;
    hold off
    subplot(1,2,2);
    hold on
    b = plot(phi, vp_sti, '-k');
    b.Color(4) = i-0.4;
    hold off
end
labels = [labels,{'best'}];
subplot(1,2,1);
hold on
plot(PHIE_data(ss), VP_data(ss), '-go', 'markersize', 6);
xlim([0, 0.4]);ylim([2e3, 4e3]);
title('Soft Sand model');
xlabel('porosity PHIE');
hold off
legend(labels);
subplot(1,2,2);
hold on
plot(PHIE_data(ss), VP_data(ss), '-go', 'markersize', 6);
xlim([0, 0.4]);ylim([2e3, 4e3]);
title('Stiff Sand model');
xlabel('porosity PHIE');
hold off
legend(labels);

% imgname = 'Figure_2_part2.png'; % '/Users/matt/Dropbox/Agile/SEG/Tutorials/delMonte_Apr2017/Figure_2_part2.png'
% print(gcf, '-dpng', '-r600', imgname);


%% 17
phi_test = [0.2,0.3];

[Kdry, MUdry] = softsand(K0, MU0, phi_test, 0.5, 12, 45);
[p_ssm, ~, ~, ~] = vels(Kdry,MUdry,K0,RHO0,K_b,RHO_b,phi_test);

[Kdry, MUdry] = stiffsand(K0, MU0, phi_test, 0.4, 8, 45);
[vp_sti, ~, ~, ~] = vels(Kdry,MUdry,K0,RHO0,K_b,RHO_b,phi_test);

fprintf('Soft Sand RPM: phi=%.2f --> Vp=%.0f m/s\n', phi_test(1), vp_ssm(1));
fprintf('Soft Sand RPM: phi=%.2f --> Vp=%.0f m/s\n', phi_test(2), vp_ssm(2));
fprintf('SoftSand RPM:increase in Vp after 10%% porosity decrease=%d%%\n',fix((vp_ssm(1)-vp_ssm(2))/vp_ssm(2)*100));

fprintf('Stiff Sand RPM: phi=%.2f --> Vp=%.0f m/s\n', phi_test(1), vp_sti(1));
fprintf('Stiff Sand RPM: phi=%.2f --> Vp=%.0f m/s\n', phi_test(2), vp_sti(2));
fprintf('Stiff Sand RPM: increase in Vp after 10%% porosity decrease=%d%%\n',fix((vp_sti(1)-vp_sti(2))/vp_sti(1)*100));


%% 19,20
[ip_rpt1, vpvs_rpt1] = rpt('stiff', 0.6, 'oil', 0.4, 8, 45, 0.3);

[ip_rpt0, vpvs_rpt0] = rpt('stiff', 0.6, 'oil', 0.4, 8, 45, 0.3);


%% 21
figure()
hold on
plot(ip_rpt0, vpvs_rpt0, '-sk');
plot(VP_data(ss).*RHO_data(ss), VP_data(ss)./VS_data(ss), '-go', 'markersize', 6);
alpha(0.5);
xlim([2e3,11e3]);ylim([1.6,2.8]);
xlabel('IP');ylabel('VP/VS');
set(gca,'linewidth',0.1)

% imgname = 'Figure_3_part2.png'; % '/Users/matt/Dropbox/Agile/SEG/Tutorials/delMonte_Apr2017/Figure_3_part2.png';
% print(gcf, '-dpng', '-r600', imgname);


%% 22
% Original parameters.
RHO_sh=2.8; K_sh=15; MU_sh=5;

[ip_rpt0a,vpvs_rpt0a] = rpt('soft', 0.0, 'oil', 0.5, 12, 45, 0.3, false);
[ip_rpt0b,vpvs_rpt0b] = rpt('soft', 0.6, 'oil', 0.5, 12, 45, 0.3, false);

figure()
hold on
plot(ip_rpt0a, vpvs_rpt0a, '-sk');
plot(ip_rpt0b, vpvs_rpt0b, '-sr');
plot(VP_data(ss).*RHO_data(ss), VP_data(ss)./VS_data(ss), '-go', 'markersize', 6);
alpha(0.5);
xlim([2e3,11e3]);ylim([1.6,2.8]);
set(gca,'linewidth',0.1)

% imgname = 'Figure_3_part2.png'; % '/Users/matt/Dropbox/Agile/SEG/Tutorials/delMonte_Apr2017/Figure_3_part2.png';
% print(gcf, '-dpng', '-r600', imgname);


%% 23
% New parameters.
RHO_sh=2.8; K_sh=21; MU_sh=3;

[ip_rpt0a,vpvs_rpt0a] = rpt('soft', 0.0, 'oil', 0.5, 12, 45, 0.3, false);
[ip_rpt0b,vpvs_rpt0b] = rpt('soft', 0.6, 'oil', 0.5, 12, 45, 0.3, false);

figure()
hold on
% plot(ip_rpt0a, vpvs_rpt0a, '-sk', 'markeredgewidth', 0);
plot(ip_rpt0b, vpvs_rpt0b, '-sb');
plot(VP_data(ss).*RHO_data(ss), VP_data(ss)./VS_data(ss), '-go', 'markersize', 6);
alpha(0.5);
xlim([2e3,11e3]);ylim([1.6,2.8]);
xlabel('IP');ylabel('VP/VS');
set(gca,'linewidth',0.1)

% imgname = 'Figure_3_part3.png'; % '/Users/matt/Dropbox/Agile/SEG/Tutorials/delMonte_Apr2017/Figure_3_part3.png'
% print(gcf, '-dpng', '-r600', imgname);


%% 25,26,27
fm = 400000; %主频
t = 0.25:0.001:25;
% ricker
wavelet = (1-2*(pi*fm*(t-0.02)*1e-3).^2) .* exp(-(pi*fm*(t-0.02)*1e-3).^2);
wavelet = wavelet(1:500);
fprintf('%d\n',length(wavelet));

top=2175;
z0=top-40;
z1=top;
z2=top+40;
z=index;
ss=(z>=z1) & (z<=z2);
sh=(z>=z1) & (z<=z2);

vp0 = mean(VP_data(sh));vs0 = mean(VP_data(sh));rho0 = mean(VP_data(sh));
vp1 = mean(VP_data(ss));vs1 = mean(VP_data(ss));rho1 = mean(VP_data(ss));
samples_sh=length(isfinite(index(sh)));
samples_ss=length(isfinite(index(ss)));

% fprintf('(%.0f-%.0f,%d samples) Vp=%.0f, Vs=%.0f, rho=%.2f, Ip=%.0f, Vp/Vs=%.2f\n',z0,z1,samples_sh,vp0,vs0,rho0,vp0*rho0,vp0/vs0);
% fprintf('(%.0f-%.0f,%d samples) Vp=%.0f, Vs=%.0f, rho=%.2f, Ip=%.0f, Vp/Vs=%.2f\n',z1,z2,samples_ss,vp1,vs1,rho1,vp1*rho1,vp1/vs1);

twolayer(vp0, vs0, rho0, vp1, vs1, rho1);


%% 28
RHO_qz=2.65;   K_qz=37;  MU_qz=15;
RHO_sh=2.8;   K_sh=15;  MU_sh=5;

phi=0.15;
nn=0.9;
[~, ~, K0] = vrh(nn,K_qz,K_sh);
[~, ~, MU0] = vrh(nn,MU_qz,MU_sh);   
RHO0 = nn*RHO_qz+(1-nn)*RHO_sh;

[Kdry, MUdry] = softsand(K0, MU0, phi, 0.5, 12, 45, 0.3);
[vp_rpm, vs_rpm, rho_rpm, ~]= vels(Kdry,MUdry,K0,RHO0,K_g,RHO_g,phi);

twolayer(vp0, vs0, rho0, vp_rpm, vs_rpm, rho_rpm);
title('GAS case, phi=0.15')

% imgname = 'Figure_4.png'; % '/Users/matt/Dropbox/Agile/SEG/Tutorials/delMonte_Apr2017/Figure_4.png'
% print(gcf, '-dpng', '-r600', imgname);


%% 29
phi=0.25;
nn=0.9;
[~, ~, K0] = vrh(nn,K_qz,K_sh);
[~, ~, MU0] = vrh(nn,MU_qz,MU_sh);   
RHO0 = nn*RHO_qz+(1-nn)*RHO_sh;

[Kdry, MUdry] = softsand(K0, MU0, phi, 0.5, 12, 45, 0.3);
[vp_rpm, vs_rpm, rho_rpm, ~]= vels(Kdry,MUdry,K0,RHO0,K_g,RHO_g,phi);

twolayer(vp0, vs0, rho0, vp_rpm, vs_rpm, rho_rpm);
title('GAS case, phi=0.25')


%% 功能函数
% 3
function plotlog(L, z1, z2, cutoff_sand, cutoff_shale, imgname)
data = L.data;
textdata = L.textdata;

index = data(:,(textdata=="DEPTH"));
VSH_data = data(:,(textdata=="VSH"));
IP_data = data(:,(textdata=="IP"));
VPVS_data = data(:,(textdata=="VPVS"));
PHIE_data = data(:,(textdata=="PHIE"));
VP_data = data(:,(textdata=="VP"));
RHO_data = data(:,(textdata=="RHO"));
VS_data = data(:,(textdata=="VS"));

% define filters to select sand (ss) and shale (sh)
ss=(index>=z1) & (index<=z2) & (VSH_data<=cutoff_sand);
sh=(index>=z1) & (index<=z2) & (VSH_data>=cutoff_shale);

index_ss = index(ss);
index_sh = index(sh);

% plot figure
figure()

subplot(1,5,1);
hold on
plot(VSH_data(ss), index_ss, '-go', 'markersize', 6);
plot(VSH_data(sh), index_sh, '-ro', 'markersize', 6);
alpha(.5);
plot(VSH_data, index, '-k', 'linewidth', 1);
plot([0,1], [2153,2153], '-k');
plot([0,1], [2183,2183], '--k');
hold off
xlabel('VSH', 'fontsize', 8);
ylim([z1,z2]);
set(gca,'YTick',linspace(z1,z2,5));

subplot(1,5,2);
hold on
plot(IP_data(ss), index_ss, '-go', 'markersize', 6);
plot(IP_data(sh), index_sh, '-ro', 'markersize', 6);
alpha(.5);
plot(IP_data, index, '-k', 'linewidth', 1);
hold off
xlabel('IP', 'fontsize', 8);
xlim([4e3,8e3]);ylim([z1,z2])
set(gca,'XTick',linspace(4e3,8e3,3),'YTickLabel',{});

subplot(1,5,3);
hold on
plot(VPVS_data(ss), index_ss, '-go', 'markersize', 6);
plot(VPVS_data(sh), index_sh, '-ro', 'markersize', 6);
alpha(.5);
plot(VPVS_data, index, '-k', 'linewidth', 1);
hold off
xlabel('VPVS', 'fontsize', 8);
xlim([1.5,3]);ylim([z1,z2])
set(gca,'XTick',linspace(1.5,3,3),'YTickLabel',{});

subplot(1,5,4);
plot(PHIE_data(ss), VP_data(ss), '-go', 'markersize', 6);
alpha(.5);
xlabel('VP vs phie', 'fontsize', 8);
xlim([0,0.4]);ylim([2e3,4e3]);

subplot(1,5,5);
hold on
% ax4.plot(L.VP*L.RHO[ss], L.VP/L.VS[ss], **sty1) VP,RHO(ss)不等长向量，如何点乘？
plot(VP_data(ss).*RHO_data(ss), VP_data(ss)./VS_data(ss), '-go', 'markersize', 6);
plot(VP_data(sh).*RHO_data(sh), VP_data(sh)./VS_data(sh), '-ro', 'markersize', 6);
plot(VPVS_data(sh), index_sh, '-ro', 'markersize', 6);
alpha(.5);
hold off
xlabel('VP VS vs IP', 'fontsize', 8);
xlim([4e3,8e3]);ylim([1.5,3]);

% plt.subplots_adjust(wspace=.8,left=0.05,right=0.95)

print(gcf, '-dpng', '-r600', imgname);
end

% 6
function [M_Voigt, M_Reuss, M_VRH]=vrh(f, M1, M2)
%     Simple Voigt-Reuss-Hill bounds for 2-components mixture, (C) aadm 2017
% 
%     INPUT
%     f: volumetric fraction of mineral 1
%     M1: elastic modulus mineral 1
%     M2: elastic modulus mineral 2
% 
%     OUTPUT
%     M_Voigt: upper bound or Voigt average
%     M_Reuss: lower bound or Reuss average
%     M_VRH: Voigt-Reuss-Hill average
M_Voigt = f.*M1 + (1-f).*M2;
M_Reuss = 1./ ( f./M1 + (1-f)./M2 );
M_VRH   = (M_Voigt+M_Reuss)/2;
end

% 7
function [vp, vs, rho, K]=vels(K_DRY, G_DRY, K0, D0, Kf, Df, phi)
%     Calculates velocities and densities of saturated rock via Gassmann equation, (C) aadm 2015
% 
%     INPUT
%     K_DRY,G_DRY: dry rock bulk & shear modulus in GPa
%     K0, D0: mineral bulk modulus and density in GPa
%     Kf, Df: fluid bulk modulus and density in GPa
%     phi: porosity
rho = D0.*(1-phi)+Df.*phi;
K = K_DRY + (1-K_DRY./K0).^2 ./ ( (phi./Kf) + ((1-phi)./K0) - (K_DRY./(K0.^2)) );
vp = sqrt((K+4/3*G_DRY)./rho)*1e3;
vs = sqrt(G_DRY./rho)*1e3;
end

% 8
function [K_HM, G_HM]=hertzmindlin(K0, G0, phi, phic, Cn, P, f) % 原代码中没用到phi
%     Hertz-Mindlin model
%     written by aadm (2015) from Rock Physics Handbook, p.246
% 
%     INPUT
%     K0, G0: mineral bulk & shear modulus in GPa
%     phi: porosity
%     phic: critical porosity (default 0.4)
%     Cn: coordination nnumber (default 8.6)
%     P: confining pressure in MPa (default 10)
%     f: shear modulus correction factor
%        1=dry pack with perfect adhesion
%        0=dry frictionless pack
if nargin==3
    phic=0.4; Cn=8.6; P=10; f=1;
elseif nargin==4
    Cn=8.6; P=10; f=1;
elseif nargin==5
    P=10; f=1;
elseif nargin==6
    f=1;
end
P = P/1e3; % converts pressure in same units as solid moduli (GPa)
PR0 =(3*K0-2*G0)./(6*K0+2*G0); % poisson's ratio of mineral mixture
K_HM = (P*(Cn^2*(1-phic)^2*G0.^2) ./ (18*pi^2*(1-PR0).^2)).^(1/3);
G_HM = ((2+3*f-PR0*(1+3*f))./(5*(2-PR0))).*((P*(3*Cn^2*(1-phic)^2*G0.^2)./(2*pi^2*(1-PR0).^2))).^(1/3);
end

% 9
function [K_DRY, G_DRY]=softsand(K0, G0, phi, phic, Cn, P, f)
%     Soft-sand (uncemented) model
%     written by aadm (2015) from Rock Physics Handbook, p.258
% 
%     INPUT
%     K0, G0: mineral bulk & shear modulus in GPa
%     phi: porosity
%     phic: critical porosity (default 0.4)
%     Cn: coordination nnumber (default 8.6)
%     P: confining pressure in MPa (default 10)
%     f: shear modulus correction factor
%        1=dry pack with perfect adhesion
%        0=dry frictionless pack
if nargin==3
    phic=0.4; Cn=8.6; P=10; f=1;
elseif nargin==4
    Cn=8.6; P=10; f=1;
elseif nargin==5
    P=10; f=1;
elseif nargin==6
    f=1;
end
[K_HM, G_HM] = hertzmindlin(K0, G0, phi, phic, Cn, P, f);
K_DRY = -4/3*G_HM+(((phi/phic)./(K_HM+4/3*G_HM))+((1-phi/phic)./(K0+4/3*G_HM))).^(-1);
tmp = G_HM/6.*((9*K_HM+8*G_HM) ./ (K_HM+2*G_HM));
G_DRY = -tmp + ((phi/phic)./(G_HM+tmp) + ((1-phi/phic)./(G0+tmp))).^(-1);
end

% 10
function [K_DRY, G_DRY]=stiffsand(K0, G0, phi, phic, Cn, P, f)
%     Stiff-sand model
%     written by aadm (2015) from Rock Physics Handbook, p.260
% 
%     INPUT
%     K0, G0: mineral bulk & shear modulus in GPa
%     phi: porosity
%     phic: critical porosity (default 0.4)
%     Cn: coordination nnumber (default 8.6)
%     P: confining pressure in MPa (default 10)
%     f: shear modulus correction factor
%        1=dry pack with perfect adhesion
%        0=dry frictionless pack
if nargin==3
    phic=0.4; Cn=8.6; P=10; f=1;
elseif nargin==4
    Cn=8.6; P=10; f=1;
elseif nargin==5
    P=10; f=1;
elseif nargin==6
    f=1;
end
[K_HM, G_HM] = hertzmindlin(K0, G0, phi, phic, Cn, P, f);
K_DRY  = -4/3*G0 + (((phi/phic)./(K_HM+4/3*G0)) + ((1-phi/phic)./(K0+4/3*G0))).^(-1);
tmp = G0/6*((9*K0+8*G0) ./ (K0+2*G0));
G_DRY = -tmp + ((phi/phic)./(G_HM+tmp) + ((1-phi/phic)./(G0+tmp))).^(-1);
end

% 18
function [xx,yy]=rpt(model, vsh, fluid, phic, Cn, P, f, display)
if nargin==0
    model='soft'; vsh=0.0; fluid='gas'; phic=0.4; Cn=8; P=10; f=1; display=true;
elseif nargin==1
    vsh=0.0; fluid='gas'; phic=0.4; Cn=8; P=10; f=1; display=true;
elseif nargin==2
    fluid='gas'; phic=0.4; Cn=8; P=10; f=1; display=true;
elseif nargin==3
    phic=0.4; Cn=8; P=10; f=1; display=true;
elseif nargin==4
    Cn=8; P=10; f=1; display=true;
elseif nargin==5
    P=10; f=1; display=true;
elseif nargin==6
    f=1; display=true;
elseif nargin==7
    display=true;
end
phi = linspace(0.01,phic,10);
sw = linspace(0,1,10);
xx=zeros(size(phi,2),size(sw,2));
yy=zeros(size(phi,2),size(sw,2));

global RHO_o K_o RHO_g K_g K_sh K_qz MU_sh MU_qz RHO_sh RHO_qz K_b RHO_b

if strcmp(fluid,'gas')
    K_hc = K_g; RHO_hc = RHO_g;
else
    K_hc = K_o; RHO_hc = RHO_o;
end

[~, ~, K0] = vrh(vsh,K_sh,K_qz);
[~, ~, MU0] = vrh(vsh,MU_sh,MU_qz);
RHO0 = vsh.*RHO_sh+(1-vsh).*RHO_qz;
if strcmp(model,'soft')
    [Kdry, MUdry] = softsand(K0,MU0,phi,phic,Cn,P,f);
elseif strcmp(model,'stiff')
    [Kdry, MUdry] = stiffsand(K0,MU0,phi,phic,Cn,P,f);
end

for i = 1:length(sw)
    val = sw(i);
    [~, K_f, ~] = vrh(val,K_b,K_hc);
    RHO_f = val*RHO_b + (1-val)*RHO_hc;
    [vp, vs, rho, ~] = vels(Kdry,MUdry,K0,RHO0,K_f,RHO_f,phi);
    xx(:,i) = vp.*rho;
    yy(:,i) = vp./vs;
end
% opt1={'backgroundcolor':'0.9'}
% opt2={'ha':'right','backgroundcolor':'0.9'}
    
if display
    figure();
    hold on
    plot(xx, yy, '-ok');
    plot(xx.', yy.', '-ok');
    hold off
    alpha(0.3)
%     for i,val in enumerate(phi):
%         plt.text(xx[i,-1]+150,yy[i,-1]+.02,'$\phi={:.02f}$'.format(val), **opt1)
%         plt.text(xx[-1,0]-200,yy[-1,0]-0.015,'$S_\mathrm{{w}}={:.02f}$'.format(sw[0]), **opt2)
%         plt.text(xx[-1,-1]-200,yy[-1,-1]-0.015,'$S_\mathrm{{w}}={:.02f}$'.format(sw[-1]),**opt2)
    xlabel('IP');ylabel('VP/VS');
    xlim([2e3,11e3]);ylim([1.6,2.8]);
    title(['RPT ',upper(model),' (N:G=',num2str(1-vsh),', fluid=',fluid,')']);
    % imgname = 'Figure_3_part1.png';% '/Users/matt/Dropbox/Agile/SEG/Tutorials/delMonte_Apr2017/Figure_3_part1.png';
    % print(gcf, '-dpng', '-r600', imgname);
end
end

function avo=shuey2(vp1, vs1, rho1, vp2, vs2, rho2, ang)
 theta1 = real(ang);
 
drho = rho2-rho1;
dvp = vp2-vp1;
dvs = vs2-vs1;
rho = (rho1+rho2)/2;
vp = (vp1+vp2)/2;
vs = (vs1+vs2)/2;

r0 = 0.5 * (dvp./vp + drho./rho);
g = 0.5 * dvp./vp - 2 * (vs.^2./(vp.^2)) * (drho./rho + 2 * dvs./vs);
f = 0.5 * dvp./vp;

term1 = r0;
term2 = g .* (sin(theta1).^2);
term3 = f .* (tan(theta1).^2 - sin(theta1).^2);

avo = squeeze(term1 + term2 + term3);
end

% 24
function twolayer(vp0,vs0,rho0,vp1,vs1,rho1)
n_samples = 500;
interface = fix(n_samples/2);
ang = 0:30;

fm = 400000; %主频
t = 0.25:0.001:25;
% ricker
wavelet = (1-2*(pi*fm*(t-0.02)*1e-3).^2) .* exp(-(pi*fm*(t-0.02)*1e-3).^2);
wavelet = wavelet(1:500);
fprintf('%d\n',length(wavelet));
model_ip = zeros(1,n_samples);
model_vpvs = zeros(1,n_samples);
rc0 = zeros(1,n_samples);
rc1 = zeros(1,n_samples);
model_z = 0:n_samples-1;
model_ip(1:interface) = vp0.*rho0;
model_ip(interface+1:end) = vp1.*rho1;
model_vpvs(1:interface) = vp0./vs0;
model_vpvs(interface+1:end) = vp1./vs1;

% 找不到bruges.reflection.shuey2的文档，根据bruges.reflection.reflection.shuey编写shuey
avo = shuey2(vp0,vs0,rho0,vp1,vs1,rho1,ang);
rc0(interface) = avo(1);
rc1(interface) = avo(end);
synt0 = conv(rc0,wavelet,'same');
synt1 = conv(rc1,wavelet,'same');
clip = max(abs([synt0, synt1]));
% 运行到这clip=0!!!
clip = clip*1.2 + eps;

figure()
subplot(2,3,1);
plot(model_ip, model_z, 'b', 'linewidth', 4);
xlabel('IP', 'fontsize', 8);
ylim([150,350]);
set(gca, 'YTickLabel',{});
% set(gca,'XTick',linspace(min(model_ip),max(model_ip),3), 'YTickLabel',{});

subplot(2,3,2);
plot(model_vpvs, model_z, 'b', 'linewidth', 4);
xlabel('VP/VS', 'fontsize', 8);
ylim([150,350]);
set(gca, 'YTickLabel',{});
% set(gca,'XTick',linspace(min(model_vpvs),max(model_vpvs),3), 'YTickLabel',{});

subplot(2,3,3);
plot(synt0, model_z, 'k', 'linewidth', 2);
X=synt0(synt0>0);Y=model_z(synt0>0);
a = fill(X, Y, 'k', 'linewidth', 0);
alpha(a,0.5);
% xlim([-clip,clip]);
ylim([150,350]);
xlabel(['angle=',sprintf('%.0f',ang(1))], 'fontsize', 8);
set(gca, 'YTickLabel',{});
% set(gca,'XTick',linspace(min(model_vpvs),max(model_vpvs),3), 'YTickLabel',{});

subplot(2,3,4);
plot(synt1, model_z, 'k', 'linewidth', 2);
X=synt1(synt1>0);Y=model_z(synt1>0);
a = fill(X, Y, 'k', 'linewidth', 0);
alpha(a,0.5);
% xlim([-clip,clip]);
ylim([150,350]);
xlabel(['angle=',sprintf('%.0f',ang(end))], 'fontsize', 8);
set(gca, 'YTickLabel',{});
% set(gca,'XTick',linspace(min(model_ip),max(model_ip),3), 'YTickLabel',{});

subplot(2,3,5);
hold on
plot(ang, avo, 'b', 'linewidth', 4)
plot([min(ang),max(ang)], [0,0], 'k', 'linewidth', 2);
hold off
xlabel('angle of incidence', 'fontsize', 8)

% plt.subplots_adjust(wspace=.8,left=0.05,right=0.95)
end

% 30
function [K_DRY, G_DRY]=critpor(K0, G0, phi, phic)
%     Critical porosity, Nur et al. (1991, 1995)
%     written by aadm (2015) from Rock Physics Handbook, p.353
% 
%     INPUT
%     K0, G0: mineral bulk & shear modulus in GPa
%     phi: porosity
%     phic: critical porosity (default 0.4)
if nargin==3
    phic=0.4;
end
K_DRY  = K0 .* (1-phi/phic);
G_DRY  = G0 .* (1-phi/phic);
end

function [K_DRY, G_DRY]=contactcement(K0, G0, phi, phic, Cn, Kc, Gc, scheme)
%     Contact cement (cemented sand) model, Dvorkin-Nur (1996)
%     written by aadm (2015) from Rock Physics Handbook, p.255
% 
%     INPUT
%     K0, G0: mineral bulk & shear modulus in GPa
%     phi: porosity
%     phic: critical porosity (default 0.4)
%     Cn: coordination nnumber (default 8.6)
%     Kc, Gc: cement bulk & shear modulus in GPa
%             (default 37, 45 i.e. quartz)
%     scheme: 1=cement deposited at grain contacts
%             2=uniform layer around grains (default)
if nargin==3
    phic=0.4; Cn=8.6; Kc=37; Gc=45; scheme=2;
elseif nargin==4
    Cn=8.6; Kc=37; Gc=45; scheme=2;
elseif nargin==5
    Kc=37; Gc=45; scheme=2;
elseif nargin==6
    Gc=45; scheme=2;
elseif nargin==7
    scheme=2;
end
PR0 = (3*K0-2*G0)./(6*K0+2*G0);
PRc = (3*Kc-2*Gc)./(6*Kc+2*Gc);
if scheme == 1 % scheme 1: cement deposited at grain contacts
    alpha = ((phic-phi)./(3*Cn*(1-phic))).^(1/4);
else % scheme 2: cement evenly deposited on grain surface
    alpha = ((2*(phic-phi))./(3*(1-phic))).^(1/2);
end
LambdaN = (2*Gc.*(1-PR0).*(1-PRc)) ./ (pi*G0.*(1-2*PRc));
N1 = -0.024153*LambdaN.^(-1.3646);
N2 = 0.20405*LambdaN.^(-0.89008);
N3 = 0.00024649*LambdaN.^(-1.9864);
Sn = N1.*alpha.^2 + N2.*alpha + N3;
LambdaT = Gc./(pi*G0);
T1 = -10^(-2)*(2.26*PR0.^2+2.07*PR0+2.3).*LambdaT.^(0.079*PR0.^2+0.1754*PR0-1.342);
T2= (0.0573*PR0.^2+0.0937*PR0+0.202)*LambdaT.^(0.0274*PR0.^2+0.0529*PR0-0.8765);
 T3=10.^(-4)*(9.654*PR0.^2+4.945*PR0+3.1)*LambdaT.^(0.01867*PR0.^2+0.4011*PR0-1.8186);
St = T1.*alpha.^2 + T2.*alpha + T3;
K_DRY = 1/6*Cn.*(1-phic).*(Kc+(4/3)*Gc).*Sn;
G_DRY = 3/5*K_DRY+3/20*Cn.*(1-phic).*Gc.*St;
end