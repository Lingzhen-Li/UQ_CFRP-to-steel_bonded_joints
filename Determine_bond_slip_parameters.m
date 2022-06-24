%% correlation, shear strength and tensile strength
figure(6)
scatter(tensile_strength, shear_strength,'ro');
hold on
linearfit = fitlm(tensile_strength, shear_strength,'Intercept',false);
kl=table2array(linearfit.Coefficients);
lim = [0 40];
plot(lim,kl(1)*lim,'k--');
xlim(lim); ylim(lim);
xlabel('Tensile strength of adhesive');
ylabel('Maximum shear stress in the bond');
% title('Tensile strength vs. shear strength');
text(25,10,"$\tau_{max}$ = " + sprintf('%.2f',kl(1))+"$\cdot \sigma_{T,a}$");
R2_ratio = 1 - sum((tensile_strength - shear_strength).^2)...
    /sum((tensile_strength - mean(tensile_strength)).^2);
text(25,7,"$R^2$ = " + sprintf('%.4f', linearfit.Rsquared.Adjusted));
hold off
box on

% saveas(gcf,'Shear_tensile_stress','epsc')

%% Histogram/pdf, shear strength/tensile strength
ratio_adhe_strength = shear_strength./tensile_strength;

% nbins = ceil(sqrt(length(shear_strength))); %Square-root choice
% nbins = ceil(1+3.32*log(length(tensile_strength))/log(10)); % Sturges' rule
nbins = ceil(2*(length(shear_strength))^(1/3)); % Rice Rule

[pHat, pCI] = lognfit(ratio_adhe_strength);
x = 0.6:0.01:1.4;
pd1 = makedist('lognormal','mu',pHat(1),'sigma',pHat(2)); %mu and sigma are...
     % the same as those in normal distribution
y1 = pdf(pd1,x); % evaluate pdf values following lognormal distribution
pd2 = makedist('normal','mu',mean(ratio_adhe_strength),'sigma',std(ratio_adhe_strength));
y2 = pdf(pd2,x); % evaluate pdf values folloing normal distribution

figure(7)
histogram(ratio_adhe_strength,nbins,'Normalization','pdf');
hold on
plot(x,y1); 
hold on
% plot(x,y2);
% hold on
text(1.2,2,"$\mu$ = " + sprintf('%.2f',mean(ratio_adhe_strength)));
text(1.2,1.8,"$\sigma$ = " + sprintf('%.2f',std(ratio_adhe_strength)));

xlabel('Ratio of shear strength on adhesive tensile strength');
ylabel('Probability density');
legend('Normalized histogram','Lognormal distribution');
% legend('Normalized histogram','Lognormal distribution','Normal distribution');
hold off

% saveas(gcf,'Stress_ratio_distribution','epsc')

%% Determination of parameter A
n_MC = 10000; % MC runs for each specimen

for i=1:length(F_PCE)/2 % loop over specimens
    tau_max(i,:) = input(i,7) .* lognrnd(pHat(1),pHat(2),[1,n_MC]); % generate random tau around sigma_max
    A(i,:) = F_PCE(i)*1e3./(4*tau_max(i,:)*input(i,2));
end

%% Bond length and bond capacity ratio based on the estimated bond length, MC
Leff = 2*A(1:97,:)*log(39);

[Lm,I] = sort(mean(Leff,2));
nl = (1:97);
Leff_p = Leff(I,:)'; %reorder using Leff from lower to higher
Leff_exp2 = Leff_exp(I);
Leff_th = 250;
for i=1:97
    mu2sigma = mean(Leff_p(:,i))+2*std(Leff_p(:,i));
    if mu2sigma<=Leff_th
        Leff_sg(i) = mu2sigma;
    else
        Leff_sg(i) = Leff_th;
    end
end

figure(8)
H1 = shadedErrorBar(nl,Leff_p,{@mean,@(x) 2*std(x)},'lineprops','-g','patchSaturation',0.33);
hold on
H2 = shadedErrorBar(nl,Leff_p,{@mean,@(x) std(x)},'lineprops','-b','patchSaturation',0.33);
hold on
H3 = shadedErrorBar(nl,Leff_p,{@mean,@(x) 0.5*std(x)},'lineprops','-r','transparent',false,'patchSaturation',0.33);
hold on
H4 = plot(nl,Lm,'k');
hold on
H5 = plot(nl,Leff_sg,'--r','LineWidth',1.5);
hold on
H6 = plot(nl,Leff_exp2,'sm');
hold on



legend('$\mu \pm 2\cdot\sigma$','$\mu \pm \sigma$','$\mu \pm 0.5\cdot\sigma$','$\mu$','Suggested','Experiment',...
    'Location','Northwest');
xlabel("Specimen number");
ylabel("Effective bond length (mm)");
% title("Bond length estimation");
% text(40,450,"Overestimated");
% text(40,420,"(Safe Side)");
% text(300,180,"Underestimated");
% text(300,150,"(Unsafe Side)");
box on
hold off

% saveas(gcf,'Estimation_effective_bond_length','epsc')
