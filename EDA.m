%% Load Data
clear
clc

load('sample_data.mat'); % monthly sample data

n=size(data,2); % no stocks 

%% Calculate Statitics for Analysis

for i=1:n
A= data(:,i);
A(find(isnan(A)))=[];
    if ~isempty(A)
        outliers = isoutlier(A, "quartiles");
        d_means(i) = mean(A);
        d_medians(i) = median(A);
        d_std(i)= std(A);
        siqr = iqr(A);
        d_25th(i) = prctile(A, 25);
        d_75th(i) = prctile(A, 75); 
        d_min(i) = min(A(~outliers));
        d_max(i) = max(A(~outliers));
        index = A(outliers);
        LOV = index(A(outliers)>0);
        SOV = index(A(outliers)<0);
        iqrs_A(1:length(LOV),i) = (LOV - d_75th(i))/siqr;
        iqrs_B(1:length(SOV),i) = (SOV + d_25th(i))/siqr;
        LO = (A(outliers)>0);
        SO = (A(outliers)<0);
        sum_out(i) = length(index);
        d_ave_outliers_L(i) = mean(index(LO));
        d_ave_outliers_S(i) = mean(index(SO));
        d_counts(i) = nnz(~isnan(A));
        d_skewness(i) = skewness(A);
        d_kurtosis(i) = kurtosis(A);
    end
end

iqrs_A(iqrs_A==0)=NaN;
iqrs_B(iqrs_B==0)=NaN;

%% Histogram of counts 
% Sample data all has same history period

histogram(d_counts,20, 'FaceColor', [0 0.4470 0.7410])
hold on
ylabel('Number of Stocks')
xlabel('Number of Months')


%% Monthly Overview Scatter Chart

figure();
ax = axes;
scatter(1:n, d_25th*100, [], [0.4660 0.6740 0.1880], 'filled')
hold on
scatter(1:n, d_max*100, [], [0 0 0], 'filled')
scatter(1:n, d_75th*100, [], [0.4660 0.6740 0.1880], 'filled')
scatter(1:n, d_min*100, [], [0 0 0], 'filled')
scatter(1:n, d_medians*100, [], [0 0.4470 0.7410], 'filled')
ytickformat(ax, '%g%%');
xlim([1 n])
legend("25th/75th Percentile", "Minimum/Maximum","", "", "Median") 
ylabel("Monthly Return", 'Interpreter','latex', 'FontSize',12)
xlabel('Stocks 1 to n')

%% Monthly Overview Scatter Chart with Outliers

figure();
ax = axes;
scatter(1:n, d_25th*100, [], [0.4660 0.6740 0.1880], 'filled', 'MarkerFaceAlpha', 0.4);
hold on
scatter(1:n, d_75th*100, [], [0.4660 0.6740 0.1880], 'filled', 'MarkerFaceAlpha', 0.4)
scatter(1:n, d_max*100, [], [0 0 0], 'filled', 'MarkerFaceAlpha', 0.4)
scatter(1:n, d_ave_outliers_L*100, [], [1 0 0], 'filled')
scatter(1:n, d_medians*100, [], [0 0.4470 0.7410], 'filled', 'MarkerFaceAlpha', 0.4)
scatter(1:n, d_min*100, [], [0 0 0], 'filled', 'MarkerFaceAlpha', 0.4)
scatter(1:n, d_ave_outliers_S*100, [], [1 0 0], 'filled')
ytickformat(ax, '%g%%');
%ylim([-100 150])
xlim([1 n])
legend("25th/75th Percentile","", "Minimum/Maximum", "Outliers", "Median") %, headers(4), headers(5), headers(6), headers(7)) %, headers(8))
ylabel("Monthly Return", 'Interpreter','latex', 'FontSize',12)
xlabel('Stocks 1 to n')

%% Outlier Investigation

% No of Outliers
figure()
histogram(sum_out,20, 'FaceColor', [0 0.4470 0.7410])
hold on
ylabel('Number of Stocks')
xlabel('Number of Outliers')


% Dispersion of Outliers
x=1:1;
figure();
headers = ["Dispersion on Upper Outliers"];
ax = axes();
hold(ax);
for i = 1:1
        cheaders = categorical(repmat(headers(i),length(iqrs_A(:)),1));
        boxchart(cheaders,iqrs_A(:)) %, 'BoxFaceColor', colors(i,:)
end
ylim([1 26])
ylabel('No. of IQRs above 75th Percentile')

x=1:1;
figure();
headers = ["Dispersion on Lower Outliers"];
ax = axes();
hold(ax);
for i = 1:1
        cheaders = categorical(repmat(headers(i),length(iqrs_B(:)),1));
        boxchart(cheaders,iqrs_B(:)) %, 'BoxFaceColor', colors(i,:)
end
ylabel('No. of IQRs below 25th Percentile')


%% Investigation into Correlations

% Gather correlation statistics
t=1;
for i=1:n
    for j=1:n
        if i~=j
            X = data(:,[i, j]);
            X(any(isnan(X),2),:) = []; 
            if length(X)>35
                cor_matrix(i,j) = corr(X(:,1),X(:,2));
                all_corrs(t) =  corr(X(:,1),X(:,2));
                t=t+1;
            else
                cor_matrix(i,j) = "";
            end
        else
            cor_matrix(i,j) = "";
        end
    end
end

% Correlation Boxplot
x=1:1;
X = [(all_corrs)'];
figure();
headers = ["Variation of Correlations"];
ax = axes();
hold(ax);
for i = 1:1
    cheaders = categorical(repmat(headers(i),length(X),1));
    boxchart(cheaders, X(:,i)) %, 'BoxFaceColor', colors(i,:))
end

% Correlation Histogram
figure()
histogram(all_corrs,20, 'FaceColor', [0 0.4470 0.7410])
hold on
ylabel('Number of Correlations')
xlabel('Correlation Bins')

%% Boxplots of moments

% Mean Boxplot
figure()
ax= axes;
histogram(d_means*100, 50, 'FaceColor', [0 0.4470 0.7410])
hold on
ylabel('Number of Stocks')
xlabel('Mean')
xtickformat(ax, '%g%%');

% Std Dev Boxplot
figure()
ax= axes;
histogram(d_std*100, 50, 'FaceColor', [0 0.4470 0.7410])
hold on
ylabel('Number of Stocks')
xlabel('Standard Deviation')
xtickformat(ax, '%g%%');

% Skewness Boxplot
figure()
ax= axes;
histogram(d_skewness, 50, 'FaceColor', [0 0.4470 0.7410])
hold on
ylabel('Number of Stocks')
xlabel('Skewness')

% Kurtosis Boxplot
figure()
ax= axes;
histogram(d_kurtosis, 50, 'FaceColor', [0 0.4470 0.7410])
hold on
ylabel('Number of Stocks')
xlabel('Kurtosis')

