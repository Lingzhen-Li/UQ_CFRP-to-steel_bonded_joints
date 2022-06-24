%% PCE results, evaluation on all data
figure(1)
scatter(F, F_PCE, '.r');
xlabel("Tested bond capacity (kN)"); ylabel("Estimated bond capacity (kN)");
hold on 
xlim([0 160]); ylim([0 160]);
text(100,60,"$R^2$ = " + sprintf('%.4f', R2));
% title('All points included');
xl = xlim; 
plot(xl,xl,'k--');
plot(xl,1.1*xl,'b--',1.1*xl,xl,'b--');
box on
hold off

%% PCE, 10-fold validation
figure(2)
scatter(F_nfold_test, F_nfold_result, '.r');
xlabel("Tested bond capacity (kN)"); ylabel("Estimated capacity from n folds (kN)");
hold on 
xlim([0 160]); ylim([0 160]);
text(100,60,"$R^2$ = " + sprintf('%.4f', R2_nfold));
% title(n_fold+"-fold vadilation");
xl = xlim; 
plot(xl,xl,'k--');
box on
hold off

% saveas(gcf,'10-fold_validation','epsc')

%% draw histogram of F_PCE/F
ratio_F = F_PCE./F; %ratio of estimated capacity on experimental capacity, all data involved
p = 0.03; 
coeff_F = norminv(p,mean(ratio_F),std(ratio_F)); % 3% quantile, according to normal distribution 
quantile_F_table = [0.07:-0.01:0.01, 1e-3, 1e-4]';
nbins = ceil(1+3.32*log(length(ratio_F))/log(10)); % Sturges' rule
% nbins = ceil(2*(length(ratio_F))^(1/3)); % Rice Rule

figure(3)
histogram(ratio_F,nbins,'Normalization','pdf');
hold on
[pHat, pCI] = lognfit(ratio_F); % lognormal index
x = 0.5:0.01:1.5;
pd2 = makedist('normal','mu',mean(ratio_F),'sigma',std(ratio_F));
y2 = pdf(pd2,x); % evaluate pdf values folloing normal distribution

plot(x,y2); %corresponding normal distribution
hold on
plot([coeff_F coeff_F], [0 6],'--');
hold on
text(1.3,3.5,"$\mu$ = " + sprintf('%.2f',mean(ratio_F)));
text(1.3,3.3,"$\sigma$ = " + sprintf('%.2f',std(ratio_F)));

text(0.75,3.5,'$\rightarrow$');
text(0.55,3.6,"$3\%$ quantile");
text(0.55,3.4,"Ratio = " + sprintf('%.2f',coeff_F));

xlabel('Ratio of estimated bond capacity on experimental capacity');
ylabel('Probability density');
legend('Normalized histogram','Normal distribution');

hold off

coeff_F_table = round(100*icdf('Normal',quantile_F_table,mean(ratio_F),std(ratio_F)))/100;
VarNames = {'Coefficient', 'Probability of Survival'};
table(coeff_F_table,1-quantile_F_table, 'VariableNames',VarNames)

% saveas(gcf,'Force_ratio_distribution','epsc')

%% lower 97% bound of bond capacity
figure(4)
plot(polyshape([0 0 160],[0 160 160]),'FaceColor','r');
hold on
plot(polyshape([0 160 160],[0 0 160]),'FaceColor','g');
hold on
scatter(F, F_PCE*coeff_F, 'k.');
hold on
xlabel("Tested bond capacity (kN)"); ylabel("Estimated bond capacity $\times$ safety factor (kN)");
hold on 
% title("Bond capacity estimation with 97% confidence of safety");
xlim([0 160]); ylim([0 160]);
xl = xlim; 
plot(xl,xl,'k--');
text(110,45,"Underestimated");
text(110,35,"(Safe Side)");
text(30,120,"Overestimated");
text(30,110,"(Unsafe Side)");

hold off

% saveas(gcf,'Safe_bond_capacity','epsc')

%% Empirical model
% Model-1: Xia and Teng (2005) 
% Model-2: Fernando thesis
% Model-3: Wang and Wu (2018)

% Index of seen data
index_known_1 = [45:49];
index_known_2 = [1:4 7:12];
index_known_3 = [1:4 7:12 24:33 39:49 68:71];

index_all = [1:91]; 

% estimated interfacial fracture energy, empirical regression equation from Teng
Gf_empi_1 = 0.5*62*(input(:,7)./(input(:,6)/2/(1+0.3))).^0.56.*input(:,5).^0.27;
delta_f = 2*Gf_empi_1./(0.8*input(:,7));
Leff_empi = pi()/2./sqrt(0.8*input(:,7)./(input(:,3).*input(:,1).*delta_f));
corr_L = input(:,4)./Leff_empi;
for ii=1:length(corr_L)
    if corr_L(ii)>1
        corr_L(ii)=1;
    end
end

omega_empi = input(:,7).^2/2./input(:,6); 
% estimating strain energy, 1/2*sigma_max*epsilon_max

for ii=1:length(omega_empi)/2
    if (isnan(strain_energy(ii))==1)
        strain_energy(ii) = omega_empi(ii);
    end
end

Gf_empi_2 = 628 * input(:,5).^0.5.*strain_energy.^2;
Gf_empi_3 = 243 * input(:,5).^0.4.*strain_energy.^1.7;

F_empi_1 = corr_L.*input(:,2).* sqrt(2*input(:,3).*input(:,1).*Gf_empi_1)/1e3; % estimated bond capacity, in kN
F_empi_2 = corr_L.*input(:,2).* sqrt(2*input(:,3).*input(:,1).*Gf_empi_2)/1e3; % estimated bond capacity, in kN
F_empi_3 = corr_L.*input(:,2).* sqrt(2*input(:,3).*input(:,1).*Gf_empi_3)/1e3; % estimated bond capacity, in kN

R2_known_1 = 1 - sum((F(index_known_1) - F_empi_1(index_known_1)).^2)...
    /sum((F(index_known_1) - mean(F(index_known_1))).^2); % only with known data wrt. model
R2_known_2 = 1 - sum((F(index_known_2) - F_empi_2(index_known_2)).^2)...
    /sum((F(index_known_2) - mean(F(index_known_2))).^2); % only with known data wrt. model
R2_known_3 = 1 - sum((F(index_known_3) - F_empi_3(index_known_3)).^2)...
    /sum((F(index_known_3) - mean(F(index_known_3))).^2); % only with known data wrt. model
MAPE_known_1 = 1/length(index_known_1)* sum(abs((F(index_known_1) - F_empi_1(index_known_1))./F(index_known_1)));
MAPE_known_2 = 1/length(index_known_2)* sum(abs((F(index_known_2) - F_empi_2(index_known_2))./F(index_known_2)));
MAPE_known_3 = 1/length(index_known_3)* sum(abs((F(index_known_3) - F_empi_3(index_known_3))./F(index_known_3)));

for no_run = 1:100
    for i=1:10
        %index of samples not involved in training the three models
        n_unknown_1 = round(length(index_known_1)/9); index_unknown_1 = randsample(setdiff(index_all,index_known_1),n_unknown_1);
        n_unknown_2 = round(length(index_known_2)/9); index_unknown_2 = randsample(setdiff(index_all,index_known_2),n_unknown_2);
        n_unknown_3 = round(length(index_known_3)/9); index_unknown_3 = randsample(setdiff(index_all,index_known_3),n_unknown_3);

        % to reduce bias of unknown/known ratio in Xia and Teng model,
        % known data (5) is twice counted while prediction (1) is counted only once
        ind_sample_1 = [index_known_1 index_known_1 index_unknown_1];
        ind_sample_2 = [index_known_2 index_unknown_2];
        ind_sample_3 = [index_known_3 index_unknown_3];
        R2_unknown_1(i) = 1 - sum((F(ind_sample_1) - F_empi_1(ind_sample_1)).^2)...
            /sum((F(ind_sample_1) - mean(F(ind_sample_1))).^2); % contain unknown data wrt. model
        R2_unknown_2(i) = 1 - sum((F(ind_sample_2) - F_empi_2(ind_sample_2)).^2)...
            /sum((F(ind_sample_2) - mean(F(ind_sample_2))).^2); % contain unknown data wrt. model
        R2_unknown_3(i) = 1 - sum((F([index_known_3 index_unknown_3]) - F_empi_3([index_known_3 index_unknown_3])).^2)...
            /sum((F([index_known_3 index_unknown_3]) - mean(F([index_known_3 index_unknown_3]))).^2); % contain unknown data wrt. model

        MAPE_unknown_1(i) = 1/length(ind_sample_1)* sum(abs((F(ind_sample_1) - F_empi_1(ind_sample_1))./F(ind_sample_1)));
        MAPE_unknown_2(i) = 1/length(ind_sample_2)* sum(abs((F(ind_sample_2) - F_empi_1(ind_sample_2))./F(ind_sample_2)));
        MAPE_unknown_3(i) = 1/length(ind_sample_3)* sum(abs((F(ind_sample_3) - F_empi_1(ind_sample_3))./F(ind_sample_3)));
 
    end

    R2_mean_1(no_run) = mean(R2_unknown_1);
    R2_mean_2(no_run) = mean(R2_unknown_2);
    R2_mean_3(no_run) = mean(R2_unknown_3);
    
    MAPE_mean_1(no_run) = mean(MAPE_unknown_1);
    MAPE_mean_2(no_run) = mean(MAPE_unknown_2);
    MAPE_mean_3(no_run) = mean(MAPE_unknown_3);
    
end


figure(17)
plot(R2_mean_1,'b')
hold on
plot(R2_mean_2,'r')
hold on
plot(R2_mean_3,'g')
hold on
plot(R2_nfold_100,'k')
hold on

plot(MAPE_mean_1,'b--')
hold on
plot(MAPE_mean_2,'r--')
hold on
plot(MAPE_mean_3,'g--')
hold on
plot(MAPE_nfold_100,'k--')
hold on

xlim([0 no_run])
ylim([0 1])
xlabel('Number of random model runnings')
ylabel('Solid line: $R^2$ score, dashed line: MAPE')
% legend('Xia and Teng model','Fernando model',' Wang and Wu model','PCE model')

h = legend('Xia and Teng model','Fernando model',' Wang and Wu model','PCE model');
pos = get(h,'Position');
posx = 0.6;
posy = 0.35;
set(h,'Position',[posx posy pos(3) pos(4)]);

hold off

% saveas(gcf,'Model_performance_comparison','epsc')