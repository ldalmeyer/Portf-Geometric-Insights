%% Results - Group of 4 Charts
% Sharpe Ratio, Probabilistic Sharpe Ratio, Alpha-weight angle and Covariance matrix improvement

% Average Sharpe Ratio over 1000 Sims, for each K (shrinkage intensity)
headers = ["Diag ", "CCM", "OPM", "TPM", "OFMM", "G&J", "AOCS"];

figure();
ytickformat('%.2f')
plot(mean(sharpe_adj(:,:,1)), 'Color',[0 0.4470 0.7410],'LineWidth', 2.5) 
hold on
xlim([0 100])
plot(mean(sharpe_adj(:,:,2)),'Color',[0.4660 0.6740 0.1880],'LineWidth', 2.5) 
plot(mean(sharpe_adj(:,:,3)),'Color',[0.4940 0.1840 0.5560],'LineWidth', 2.5) 
plot(mean(sharpe_adj(:,:,4)),'Color',[0.3010 0.7450 0.9330],'LineWidth', 2.5) 
plot(mean(sharpe_adj(:,:,5)),'Color',[0.9290 0.6940 0.1250],'LineWidth', 2.5)
plot(mean(sharpe_adj(:,:,6)),'Color',[0 0 0],'LineWidth', 2.5)
plot(mean(sharpe_adj(:,:,7)),'Color',[0.6350 0.0780 0.1840],'LineWidth', 2.5)
legend(headers(1), headers(2), headers(3), headers(4), headers(5), headers(6), headers(7))
ylabel("Sharpe Ratio", 'Interpreter','latex', 'FontSize',12)
xlabel('k')

% Average Probabilistic Sharpe Ratio over 1000 Sims, for each K (shinkage intensity)
headers = ["Diag ", "CCM", "OPM", "TPM", "OFMM", "G&J", "AOCS"];

figure();
plot(mean(prob_sr(:,:,1)), 'Color',[0 0.4470 0.7410],'LineWidth', 2.5)
hold on
xlim([0 100])
plot(mean(prob_sr(:,:,2)),'Color',[0.4660 0.6740 0.1880],'LineWidth', 2.5) 
plot(mean(prob_sr(:,:,3)),'Color',[0.4940 0.1840 0.5560],'LineWidth', 2.5) 
plot(mean(prob_sr(:,:,4)),'Color',[0.3010 0.7450 0.9330],'LineWidth', 2.5) 
plot(mean(prob_sr(:,:,5)),'Color',[0.9290 0.6940 0.1250],'LineWidth', 2.5) 
plot(mean(prob_sr(:,:,6)),'Color',[0 0 0],'LineWidth', 2.5) 
plot(mean(prob_sr(:,:,7)),'Color',[0.6350 0.0780 0.1840],'LineWidth', 2.5) 
legend(headers(1), headers(2), headers(3), headers(4), headers(5), headers(6), headers(7))
ylabel("Probabilistic Sharpe Ratio", 'Interpreter','latex', 'FontSize',12)
xlabel('k')

% Angle Improvement
improvement_angle = mean((angles_improved_ave)); %average percent improvement of angle over 1000 sims

figure
plot(improvement_angle(:,:,1), 'Color',[0 0.4470 0.7410],'LineWidth', 2.5)
hold on
xlim([0 100])
plot(improvement_angle(:,:,2), 'Color',[0.4660 0.6740 0.1880],'LineWidth', 2.5)
plot(improvement_angle(:,:,3), 'Color',[0.4940 0.1840 0.5560],'LineWidth', 2.5)
plot(improvement_angle(:,:,4), 'Color',[0.3010 0.7450 0.9330],'LineWidth', 2.5)
plot(improvement_angle(:,:,5), 'Color',[0.9290 0.6940 0.1250],'LineWidth', 2.5)
plot(improvement_angle(:,:,6),'Color',[0 0 0],'LineWidth', 2.5)
plot(improvement_angle(:,:,7), 'Color',[0.6350 0.0780 0.1840],'LineWidth', 2.5)
legend(headers(1), headers(2), headers(3), headers(4), headers(5), headers(6), headers(7)) %, headers(8))
ylabel("$\alpha$-weight angle", 'Interpreter','latex', 'FontSize',12)
xlabel('k')

% Condition Number Improvement
cond_improve = mean(cond_improve); %average percent improvement of angle over 1000 sims

figure
plot(cond_improve(:,:,1), 'Color',[0 0.4470 0.7410],'LineWidth', 2.5)
hold on
xlim([0 100])
plot(cond_improve(:,:,2), 'Color',[0.4660 0.6740 0.1880],'LineWidth', 2.5)
plot(cond_improve(:,:,3), 'Color',[0.4940 0.1840 0.5560],'LineWidth', 2.5)
plot(cond_improve(:,:,4), 'Color',[0.3010 0.7450 0.9330],'LineWidth', 2.5)
plot(cond_improve(:,:,5), 'Color',[0.9290 0.6940 0.1250],'LineWidth', 2.5)
plot(cond_improve(:,:,6),'Color',[0 0 0],'LineWidth', 2.5)
plot(cond_improve(:,:,7), 'Color',[0.6350 0.0780 0.1840],'LineWidth', 2.5)
legend(headers(1), headers(2), headers(3), headers(4), headers(5), headers(6), headers(7)) %, headers(8))
ylabel("% Improvement", 'Interpreter','latex', 'FontSize',12)
xlabel('k')


%% Intensity Factor (k)

% Deflated Sharpe Ratio - which shrinkage intensity (k)

colors = [[0 0.4470 0.7410]; [0.4660 0.6740 0.1880]; [0.4940 0.1840 0.5560]; [0.3010 0.7450 0.9330]; [0.9290 0.6940 0.1250]; [0 0 0]; [0.6350 0.0780 0.1840]];

x=1:6;
figure();
headers = ["Diag ", "CCM", "OPM", "TPM", "OFMM", "G&J", "AOCS"];
ax = axes();
hold(ax);
for i = 1:6
    cheaders = categorical(repmat(headers(i),nsim,1));
    boxchart(cheaders, which_pvalue(:,i)-1, 'BoxFaceColor', colors(i,:))
    %boxplot(which_pvalue-1,'Labels',headers,'Colors', ['b' 'g' 'c' 'm' 'y' 'r'], 'BoxStyle', 'filled')
end
ylim([0 102])
xlabel('Shrinkage Methodology')
ylabel('Shrinkage Parameter [0;100]')


% Deflated sharpe ratio - pvalues

x=1:6;
figure();
headers =["Diag ", "CCM", "OPM", "TPM", "OFMM", "G&J", "AOCS"];
ax = axes();
hold(ax);
for i = 1:6
    cheaders = categorical(repmat(headers(i),nsim,1));
    boxchart(cheaders, max_pval(:,i), 'BoxFaceColor', colors(i,:))
end
ylim([0 1])
xlabel('Shrinkage Methodology')
ylabel('P-Value[0;1]')


%% Weight Vector Magnitude

headers = ["Diag ", "CCM", "OPM", "TPM", "OFMM", "G&J", "AOCS"];

figure();
ytickformat('%.2f')
plot(mean(nrms(:,:,1)), 'Color',[0 0.4470 0.7410],'LineWidth', 2.5) 
hold on
xlim([0 100])
plot(mean(nrms(:,:,2)),'Color',[0.4660 0.6740 0.1880],'LineWidth', 2.5)
plot(mean(nrms(:,:,3)),'Color',[0.4940 0.1840 0.5560],'LineWidth', 2.5) 
plot(mean(nrms(:,:,4)),'Color',[0.3010 0.7450 0.9330],'LineWidth', 2.5) 
plot(mean(nrms(:,:,5)),'Color',[0.9290 0.6940 0.1250],'LineWidth', 2.5) 
plot(mean(nrms(:,:,6)),'Color',[0 0 0],'LineWidth', 2.5)
plot(mean(nrms(:,:,7)),'Color',[0.6350 0.0780 0.1840],'LineWidth', 2.5)
legend(headers(1), headers(2), headers(3), headers(4), headers(5), headers(6), headers(7)) %, headers(8))
ylabel("Weight Vector Average Magnitude", 'Interpreter','latex', 'FontSize',12)
xlabel('k')


%% Performance at chosen k

which_k = 4; % chosen intensity factor must be set

colors = [[0 0.4470 0.7410]; [0.4660 0.6740 0.1880]; [0.4940 0.1840 0.5560]; [0.3010 0.7450 0.9330]; [0.9290 0.6940 0.1250]; [0 0 0]; [0.6350 0.0780 0.1840]];

x=1:7;
figure();
headers = ["Diag ", "CCM", "OPM", "TPM", "OFMM", "G&J", "AOCS"];
ax = axes();
hold(ax);
for i = 1:7
    cheaders = categorical(repmat(headers(i),nsim,1));
    boxchart(cheaders, performance(:,which_k, i), 'BoxFaceColor', colors(i,:))
end
ylim([0 10]) % some large outliers distort chart, may vary depending on solution
a=[cellstr(num2str(get(gca,'ytick')'*100))]; 
pct = char(ones(size(a,1),1)*'%'); 
new_yticks = [char(a),pct];
set(gca,'yticklabel',new_yticks)
xlabel('Shrinkage Methodology')
ylabel('5Y Performance')

%% Risk at chosen k

colors = [[0 0.4470 0.7410]; [0.4660 0.6740 0.1880]; [0.4940 0.1840 0.5560]; [0.3010 0.7450 0.9330]; [0.9290 0.6940 0.1250]; [0 0 0]; [0.6350 0.0780 0.1840]];

x=1:7;
figure();
headers = ["Diag ", "CCM", "OPM", "TPM", "OFMM", "G&J", "AOCS"];
ax = axes();
hold(ax);
for i = 1:7
    cheaders = categorical(repmat(headers(i),nsim,1));
    boxchart(cheaders, risk(:,which_k, i)*sqrt(12), 'BoxFaceColor', colors(i,:))
end
ylim([0.005 1])
a=[cellstr(num2str(get(gca,'ytick')'*100))]; 
pct = char(ones(size(a,1),1)*'%'); 
new_yticks = [char(a),pct];
set(gca,'yticklabel',new_yticks)
xlabel('Shrinkage Methodology')
ylabel('Standard Deviation (ann)')

