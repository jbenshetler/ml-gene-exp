function [ hlplot, Beta, FitInfo ] = lassoRegression( y, X, predNames )
%LASSOTABLE Perform lasso analysis and plotting
% The goal is to identify important predictors and discard those that are unnecessary.
% INPUTS
% y is a row vector of length N
% X is an NxP array, where P is the number of predictor variables
% predNames is a 1xN cell of character vectors
% OUTPUTS
% B is the fitted least-squares regression coefficients
% FitInfo is information about the fits.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by Jeff Benshetler, 2016-09-24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Based heavily on the following MATLAB examples/docs
% Machine Learning with MATLAB: Lasso Regularization
% http://www.mathworks.com/solutions/machine-learning/examples.html?file=/products/demos/machine-learning/diabetes/diabetes.html


% Perform Lasso Regularization
[Beta, FitInfo] = lasso(X,y,...
                        'Standardize',true,...
                        'CV',10,...
                        'PredictorNames',predNames,...
                        'Options',statset('UseParallel',true)...
                        );

% Larger values of lambda appear on the left side of the graph, which 
% means that there is increased regularization. As the lambda value 
% increases, the number of nonzero predictors also increases.
lassoPlot(Beta,FitInfo,'PlotType','Lambda','XScale','log');


hlplot = get(gca,'Children');

% Generating colors for each line in the plot
colors = hsv(numel(hlplot));
for ii = 1:numel(hlplot)
    set(hlplot(ii),'color',colors(ii,:));
end

set(hlplot,'LineWidth',2)
%Sset(gcf,'Units','Normalized','Position',[0.2 0.4 0.5 0.35])
legend('Location','Best')

% As a rule of thumb, one standard-error value is often used for 
% choosing a smaller model with a good fit.
lam = FitInfo.Index1SE;
isImportant =Beta(:,lam) ~= 0;
fprintf('Important Predictors\n');
disp(predNames(isImportant))

disp('Fit a Linear Model with the Terms for Comparison')
lmNames = predNames;
lmNames{end+1} = 'y';
tbl = array2table([X y],'VariableNames',lmNames);
mdlFull = fitlm(tbl,'Intercept',false);
%disp(mdlFull)

disp('Compare the MSE for regularized and unregularized models.')
disp(['Lasso MSE: ', num2str(FitInfo.MSE(lam))])
disp(['Full  MSE: ', num2str(mdlFull.MSE)])


end

