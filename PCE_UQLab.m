%% create PCE model
function [myPCE] = PCE_UQLab(X, Y) % X is model input, Y is model response

MetaOpts.Type = 'Metamodel';
MetaOpts.MetaType = 'PCE';
MetaOpts.Method = 'OLS';

MetaOpts.ExpDesign.X = X;
MetaOpts.ExpDesign.Y = Y;

MetaOpts.Degree = 2;

MetaOpts.ValidationSet.X = X;
MetaOpts.ValidationSet.Y = Y;

myPCE = uq_createModel(MetaOpts);

end