clear
clc
uqlab % activate uqlab

% improve the quality of EPS images
set(gcf,'renderer','painters')

% set plot format in accordance with latex
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex'); 
set(groot,'defaultLegendInterpreter','latex');

%% Read test data
load test_data.mat
X = input;
F = output;

name = strings(8,2);
name(1,1) = "CFRP thickness";
name(2,1) = "CFRP width";
name(3,1) = "CFRP E-modulus";
name(4,1) = "Bond length";
name(5,1) = "Adhesive thickness";
name(6,1) = "Adhesive E-modulus";
name(7,1) = "Adhesive tensile strength";
name(8,1) = "Adhesive maximum tensile strain";

name(1,2) = "(mm)";
name(2,2) = "(mm)";
name(3,2) = "(MPa)";
name(4,2) = "(mm)";
name(5,2) = "(mm)";
name(6,2) = "(MPa)";
name(7,2) = "(MPa)";
name(8,2) = "(%)";

% assuming uniform distribution in the input variables
input_lower = min(input);
input_upper = max(input);

%% PCE using UQLab
for i=1:8
    InputOpts.Marginals(i).Name = name(i,1);
    InputOpts.Marginals(i).Type = 'Uniform';
    InputOpts.Marginals(i).Parameters = [input_lower(i) input_upper(i)]; % assuming Uniform distribution
end

myInput = uq_createInput(InputOpts);

%% PCE model, all data
[myPCE] = PCE_UQLab(X, F);

uq_print(myPCE)

% Validation
F_PCE = uq_evalModel(myPCE,X);

R2 = 1 - sum((F - F_PCE).^2)/sum((F - mean(F)).^2);
MAPE = 1/length(F)* sum(abs((F - F_PCE)./F));

%% Sobol indices
SobolOpts.Type = 'Sensitivity';
SobolOpts.Method = 'Sobol';
SobolOpts.Sobol.Order = 1;
SobolAnalysis = uq_createAnalysis(SobolOpts);

uq_print(SobolAnalysis)

%% PCE model, n-fold validation
%uncomment only when estimating 100 10-fold validation
for no_run = 1:100
%uncomment only when estimating 100 10-fold validation
n_fold=10;
n_partition = cvpartition(F,'KFold',n_fold,'Stratify',false);
F_nfold_result = [];
F_nfold_test = [];


    for k=1:n_fold % k-th fold
        idx_training = training(n_partition,k); %position index for training
        idx_test = test(n_partition,k); %position index for testing

        %reset these arguments, as they change size with different validation
        F_nfold_train = []; F_nfold_val = []; input_nfold_train = []; input_nfold_val = [];

        F_nfold_train = F(idx_training); %training set
        F_nfold_val = F(idx_test); %validation set

        input_nfold_train(:,:) = X(idx_training,:); %training set
        input_nfold_val(:,:) = X(idx_test,:); %validation set

        % train PCE model without k-th fold
        % 
        [myPCE] = PCE_UQLab(input_nfold_train, F_nfold_train);

        % Validating the k-th model
        F_nfold_temp = uq_evalModel(myPCE,input_nfold_val);

        F_nfold_test = [F_nfold_test;F_nfold_val];
        F_nfold_result = [F_nfold_result; F_nfold_temp];
    end

R2_nfold = 1 - sum((F_nfold_test - F_nfold_result).^2)/sum((F_nfold_test - mean(F_nfold_test)).^2);
MAPE_nfold = 1/length(F_nfold_test)* sum(abs((F_nfold_test - F_nfold_result)./F_nfold_test));

%uncomment only when estimating 100 10-fold validation
R2_nfold_100(no_run) = 1 - sum((F_nfold_test - F_nfold_result).^2)/sum((F_nfold_test - mean(F_nfold_test)).^2);
MAPE_nfold_100(no_run) = 1/length(F_nfold_test)* sum(abs((F_nfold_test - F_nfold_result)./F_nfold_test));
end
%uncomment only when estimating 100 10-fold validation

%% Plot_PCE
Plot_PCE

%% Determine_bond-slip_parameters
Determine_bond_slip_parameters

%% MC simulation
Plot_test
