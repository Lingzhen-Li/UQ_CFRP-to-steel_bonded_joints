%% CFRP-S1-T0.5-1
clear
clc

% improve the quality of EPS images
set(gcf,'renderer','painters')

% Latex format
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex'); 
set(groot,'defaultLegendInterpreter','latex');

load CFRP-S1-T0.5-1.mat
sz = 16;

%% Force-disiplacement
figure(9)

plot(disp_DIC, F_DIC, 'r');
hold on
plot(disp_LUNA, F_LUNA, 'b--');
hold on
for i=1:n_MC
    ld = plot(s_ana, F_ana(i,:),'g.');
    ld.MarkerSize = 2+3*p_tau_max(i)/max(p_tau_max); % density of dots
    hold on
    [~,index_s(i)] = min(abs(F_ana(i,:)-0.7*max(F_ana(i,:))));
end
xlabel('Displacement of the loaded end (mm)','FontSize',sz);
ylabel('Load (kN)','FontSize',sz);
% title('F-\Delta curve');
legend('Measurement-1','Measurement-2','Monte-Carlo simulation','Location','Northwest')
xlim([0 1]); ylim([0 50]);
hold off

% saveas(gcf,'Load_displacement_linear','epsc2')
% saveas(gcf,'Load_displacement_linear','jpeg')


%% Shear stress distribution
figure(10)
for i=1:(length(x_DIC)-1)
    xx_DIC(i) = (x_DIC(i) + x_DIC(i+1))/2;
end

for i=1:(length(x_LUNA)-1)
    xx_LUNA(i) = (x_LUNA(i) + x_LUNA(i+1))/2;
end

i=200;
plot(xx_DIC, tau_DIC(i,:),'r');
hold on
plot(xx_LUNA, tau_LUNA(i+row_difference,:), 'b--');
hold on

for i=1:n_MC
    ss = plot(x_ana+254, tau_x_ana(i,:),'g.');
    ss.MarkerSize = 2+2*p_tau_max(i)/max(p_tau_max); % density of dots
    hold on
end

plot([0 300],[0 0],'k-');
hold off
xlim([0 300]);
ylim([0 35]);
xlabel('Distance to the free end (mm)','FontSize',sz);
ylabel('Shear stress (MPa)','FontSize',sz);
legend('Measurement-1','Measurement-2','Monte-Carlo simulation','Location','Northwest')
% title('Shear stress along the bond line');
% legend('Analytical','DIC','LUNA','location','northwest');
hold off

% saveas(gcf,'Shear_stress_along_bond_line_linear','epsc')
% saveas(gcf,'Shear_stress_along_bond_line_linear','jpeg')

%% Bond-slip curves, both
figure(11)
n=5; %Manually adjust the output interval
slip_DIC2 = []; tau_DIC2 = [];
slip_LUNA2 = []; tau_LUNA2 = [];

for i=1:n_MC
    bs = plot(s_ana, tau_s_ana(i,:),'g.');
    bs.MarkerSize = 2+5*p_tau_max(i)/max(p_tau_max); % density of dots
    hold on
end
h3 = plot(s_ana, tau_s_ana(n_MC,:),'g.');

for j=10:n:column_DIC
    slip_DIC2 = [slip_DIC2; slip_DIC(:,j)]; 
    tau_DIC2 = [tau_DIC2; tau_DIC(:,j)];
end

for j=50:n:column_LUNA-50
    slip_LUNA2 = [slip_LUNA2; slip_LUNA(:,j)]; 
    tau_LUNA2 = [tau_LUNA2; tau_LUNA(:,j)];
end

h1 = plot(slip_DIC2, tau_DIC2,'r.');
hold on

h2 = plot(slip_LUNA2, tau_LUNA2,'b.');
hold on


xlim([0 0.15]);
ylim([0 40]);
xlabel('Slip (mm)','FontSize',sz); 
ylabel('Shear stress (MPa)','FontSize',sz);
% title('Bond-slip from LUNA');
legend([h1,h2,h3],'Measurement-1','Measurement-2','Monte-Carlo simulation','Location','Northeast')
hold off

% saveas(gcf,'Bond_slip_curve_linear','epsc')
% saveas(gcf,'Bond_slip_curve_linear','jpeg')

%% From test
i=202;

figure(12)
plot(x_LUNA, e_CFRP_LUNA_smooth(i,:), 'k');
xlim([0,300]); ylim([0,4000]);
xlabel('Distance from the free end (mm)','FontSize',sz);
ylabel('Tensile strain of CFRP ($\mu\varepsilon$)','FontSize',sz);

% saveas(gcf,'Test_tensile_strain','epsc')

figure(13)
plot(xx_LUNA, tau_LUNA(i,:), 'k');
xlim([0 300]);
ylim([0 30]);
xlabel('Distance to the free end (mm)','FontSize',sz);
ylabel('Shear stress (MPa)','FontSize',sz);
% title('Shear stress along the bond line');
% legend('Analytical','DIC','LUNA','location','northwest');

% saveas(gcf,'Test_shear_stress','epsc')

%% CFRP-A-T0.5-1
clear

% improve the quality of EPS images
set(gcf,'renderer','painters')

load CFRP-A-T0.5-1.mat
sz = 16;

%% Force-disiplacement
figure(14)
for i=1:n_MC
    ld = plot(s_ana, F_ana(i,:),'g.');
    ld.MarkerSize = 2+3*p_tau_max(i)/max(p_tau_max); % density of dots
    hold on
    [~,index_s(i)] = min(abs(F_ana(i,:)-0.9*max(F_ana(i,:))));
end
plot(disp_DIC, F_DIC, 'r');
hold on
plot(disp_LUNA, F_LUNA, 'b--');
hold on

xlabel('Displacement of the loaded end (mm)','FontSize',sz);
ylabel('Load (kN)','FontSize',sz);
xlim([0 3]); ylim([0 140]);
hold off
% title('F-\Delta curve');

% saveas(gcf,'Load_displacement_nonlinear','epsc2')
% saveas(gcf,'Load_displacement_nonlinear','jpeg')


%% Shear stress distribution
figure(15)
for i=1:(length(x_DIC)-1)
    xx_DIC(i) = (x_DIC(i) + x_DIC(i+1))/2;
end

for i=1:(length(x_LUNA)-1)
    xx_LUNA(i) = (x_LUNA(i) + x_LUNA(i+1))/2;
end

for i=1:n_MC
    ss = plot(x_ana+253, tau_x_ana(i,:),'g.');
    ss.MarkerSize = 3+5*p_tau_max(i)/max(p_tau_max); % density of dots
    hold on
end

n=500; %Manually adjust the output step interval    
for i=650:n:row_LUNA %better start from 1
    plot(xx_DIC, tau_DIC(i,:),'r');
    hold on
    plot(xx_LUNA, tau_LUNA(i+row_difference,:), 'b--');
    hold on
end

hold off
xlim([0 300]);
ylim([0 35]);
xlabel('Distance to the free end (mm)','FontSize',sz);
ylabel('Shear stress (MPa)','FontSize',sz);
% title('Shear stress along the bond line');
% legend('DIC','LUNA','location','northwest');

% saveas(gcf,'Shear_stress_along_bond_line_nonlinear','epsc')
% saveas(gcf,'Shear_stress_along_bond_line_nonlinear','jpeg')

%% Bond-slip curve, both
figure(16)
n=10; %Manually adjust the output interval

%output by looping over columns, each curve represents behavior of a point,
%different time

for i=1:n_MC
    bs = plot(s_ana, tau_s_ana(i,:),'g.');
    bs.MarkerSize = 2+5*p_tau_max(i)/max(p_tau_max); % density of dots

    hold on
end

for j=50:n:column_DIC-10
    plot(slip_DIC(:,j), tau_DIC(:,j),'r.');
    hold on
end


for j=100:n:column_LUNA-50
    plot(slip_LUNA(:,j), tau_LUNA(:,j),'b.');
    hold on
end

% title('Bond-slip curves');
xlim([0 2]);
ylim([0 35]);
xlabel('Slip (mm)','FontSize',sz); 
ylabel('Shear stress in the bond (MPa)','FontSize',sz);
hold off

% saveas(gcf,'Bond_slip_curve_nonlinear','epsc')
% saveas(gcf,'Bond_slip_curve_nonlinear','jpeg')
