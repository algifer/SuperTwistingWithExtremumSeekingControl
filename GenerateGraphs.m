%% Datos de ejemplo
rng(1)
N = 200;
X = rand(N,2);
y = 3*X(:,1).^2 + 2*X(:,2) + 0.3*randn(N,1);

figure('Name','Regression Learner Models','Color','w')

plot_id = 1;
rows = 3; cols = 3;

%% =========================
% 1. MODELOS LINEALES
%% =========================

% Linear simple
mdl = fitlm(X(:,1),y);
subplot(rows,cols,plot_id); plot_id = plot_id + 1;
scatter(y,predict(mdl,X(:,1)),10,'filled')
title('Linear Simple','FontSize',12,'FontWeight','bold'); grid on

% Linear multiple
mdl = fitlm(X,y);
subplot(rows,cols,plot_id); plot_id = plot_id + 1;
scatter(y,predict(mdl,X),10,'filled')
title('Linear Multiple','FontSize',12,'FontWeight','bold'); grid on

% % Linear interactions
% mdl = fitlm(X,y,'interactions');
% subplot(rows,cols,plot_id); plot_id = plot_id + 1;
% scatter(y,predict(mdl,X),10,'filled')
% title('Linear + Interactions','FontSize',12,'FontWeight','bold'); grid on

% Robust linear
mdl = fitlm(X,y,'RobustOpts','on');
subplot(rows,cols,plot_id); plot_id = plot_id + 1;
scatter(y,predict(mdl,X),10,'filled')
title('Robust Linear','FontSize',12,'FontWeight','bold'); grid on

% % Stepwise
% mdl = stepwiselm(X,y);
% subplot(rows,cols,plot_id); plot_id = plot_id + 1;
% scatter(y,predict(mdl,X),10,'filled')
% title('Stepwise Linear','FontSize',12,'FontWeight','bold'); grid on

%% =========================
% 2. GLM
%% =========================
mdl = fitglm(X,y,'Distribution','normal');
subplot(rows,cols,plot_id); plot_id = plot_id + 1;
scatter(y,predict(mdl,X),10,'filled')
title('GLM','FontSize',12,'FontWeight','bold'); grid on

%% =========================
% 3. REGRESSION TREES
%% =========================

mdl = fitrtree(X,y,'MaxNumSplits',10);
subplot(rows,cols,plot_id); plot_id = plot_id + 1;
scatter(y,predict(mdl,X),10,'filled')
title('Tree (Simple)','FontSize',12,'FontWeight','bold'); grid on

mdl = fitrtree(X,y,'MaxNumSplits',100);
subplot(rows,cols,plot_id); plot_id = plot_id + 1;
scatter(y,predict(mdl,X),10,'filled')
title('Tree (Deep)','FontSize',12,'FontWeight','bold'); grid on

%% =========================
% 4. ENSEMBLES
%% =========================

mdl = fitrensemble(X,y,'Method','Bag');
subplot(rows,cols,plot_id); plot_id = plot_id + 1;
scatter(y,predict(mdl,X),10,'filled')
title('Bagged Trees','FontSize',12,'FontWeight','bold'); grid on

mdl = fitrensemble(X,y,'Method','LSBoost');
subplot(rows,cols,plot_id); plot_id = plot_id + 1;
scatter(y,predict(mdl,X),10,'filled')
title('Boosted Trees','FontSize',12,'FontWeight','bold'); grid on

%% =========================
% 5. SVM
%% =========================

mdl = fitrsvm(X,y,'KernelFunction','linear');
subplot(rows,cols,plot_id); plot_id = plot_id + 1;
scatter(y,predict(mdl,X),10,'filled')
title('SVM Linear','FontSize',12,'FontWeight','bold'); grid on

mdl = fitrsvm(X,y,'KernelFunction','polynomial');
subplot(rows,cols,plot_id); plot_id = plot_id + 1;
scatter(y,predict(mdl,X),10,'filled')
title('SVM Polynomial','FontSize',12,'FontWeight','bold'); grid on

mdl = fitrsvm(X,y,'KernelFunction','gaussian');
subplot(rows,cols,plot_id); plot_id = plot_id + 1;
scatter(y,predict(mdl,X),10,'filled')
title('SVM Gaussian','FontSize',12,'FontWeight','bold'); grid on

%% =========================
% 6. GPR
%% =========================

mdl = fitrgp(X,y,'KernelFunction','squaredexponential');
subplot(rows,cols,plot_id); plot_id = plot_id + 1;
scatter(y,predict(mdl,X),10,'filled')
title('Gaussian Process','FontSize',12,'FontWeight','bold'); grid on

%% =========================
% 7. NEURAL NETWORK
%% =========================
% 
% mdl = fitrnet(X,y);
% subplot(rows,cols,plot_id); plot_id = plot_id + 1;
% scatter(y,predict(mdl,X),10,'filled')
% title('Neural Network','FontSize',12,'FontWeight','bold'); grid on

%% =========================
% 8. KERNEL APPROXIMATION
%% =========================

mdl = fitrkernel(X,y);
subplot(rows,cols,plot_id);
scatter(y,predict(mdl,X),10,'filled')
title('Kernel Approximation','FontSize',12,'FontWeight','bold'); grid on
