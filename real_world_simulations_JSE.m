%% Compute Efficient Frontiers
clearvars 
clc;

%% Data Load

%cd 'C:\Users\lara\OneDrive - Abax Investments\Desktop\PHD'
load('sample_data.mat'); % monthly data, simulated from ND Data to demonstrate results. 247 returns (January 1996 to August 2016), for 194 names

warning('off', 'all'); %removes warning from poorly conditioned matrices that arent used (checks rcond)

%% Compute the key variables

% Parameters varied in the thesis
for_ExpRet = 1; % Expected Returns are foreeable (=1), or Expected Returns are unknown (historical mean is used) (=0)
solution = "con_7"; % choose "uncon" (unconstrained or constrained solutions V and VIII), "con_6" or "con_7" (constrained solutions VI or VII) 

% Parameters unvaried for the thesis
roll = 5; % number of years to calculate stock mean and covariance matrices off of
tp = 10; % size of total performance period used (5 years to calcualte mean and covariance; 5 years to compare performance from)
nsim = 5; % number of simulations
ns = 10; % number of stocks in simulation

n = size(data,1);

% Parameters for constrained solutions VI or VII
if solution == "con_6"
    alpha_0 = 0.01; % minimum ER of portfolio
    g_0 = 1; % sum of weights
elseif solution == "con_7"
    lambda = 1; % risk aversion parameter
    g_0 = 1; % sum of weights
end

for q = 1:7 %number of tested shrinkage methodologies
t=1; % tracking the number of simulations
    while t <= nsim
        
        %simulate random 10 stock data
        starts = randi([1,n-(12*tp-1)]);
        ends = starts + (tp*12-1);

        data_sub = data(starts:ends,:); % construct 10Y data history
        data_sub = data_sub(:,all(~isnan(data_sub))); % remove stocks where there is not data available in chosen random history (wont be the case with simulated ND data)
          
        r10 = randi([1  size(data_sub,2)],1,ns);  % select random 10 stocks
        Data = data_sub(:,r10); 
       
        k = 1; % keep track of shrinkage intensity k
        for z= 0:1:100 % z is the shrinkage intensity, z=0 is closest to origional cov matrix. When z=0 we get t_a soln, but when z=100, we get ExpRet unitized sln

            for i=1:(tp*12-(12*roll)-1) % use 5Y period to calculate covariance matrix and expected returns. Keep stepping 1M forward.

                Data_R = Data(i:i+(12*roll)-1,:); % 5Y historical data
                one = ones(size(Data_R,2),1); % ones vector
                covar_orig = cov(Data_R); % origional historical covariance matrix

                if q==1
                    I_adj = covDiag(Data_R,0); % Diag Shrinkage
                elseif q==2
                    I_adj = covCor(Data_R,0); % CCM Shrinkage
                elseif q==3
                    I_adj = cov1Para(Data_R,0); % OPM Shrinkage
                elseif q==4
                    I_adj = cov2Para(Data_R,0); % TPM Shrinkage
                elseif q==5
                    I_adj = covMarket(Data_R,0); % OFMM Shrinkage
                end
                % q = 6 G&J Shrinkage
                % q = 7 AOCS Shrinkage

                % Which ExpRet to use
                if for_ExpRet == 1 % Expected returns are forseeable
                    ExpRet = Data(i+(12*roll)+1,:)';
                elseif for_ExpRet == 0 
                    ExpRet = mean(Data_R)'; % Using historical means, assume no forecast edge
                end
                exp_ret_fwd = Data(i+(12*roll)+1,:)'; % used to calculate alpha-weight angle

                I = eye(size(Data_R,2));
                
                %__________________________________________________________
                % Origional Solution
               
                A = one'*inv(covar_orig)*one;
                B = transpose(one)*inv(covar_orig)*ExpRet;
                C = ExpRet'*inv(covar_orig)*ExpRet;
                D = (A*C) - B^2;

                %Constrained solutions VI or VII
                if (solution == "con_6" || solution == "con_7")
                    ta = inv(covar_orig)*ExpRet*(1/(B));
                    t0 = inv(covar_orig)*one.*(1/(one'*inv(covar_orig)*one));
                    if solution == "con_6"
                        iw = (B/D)*((alpha_0*A)-(g_0*B));
                        t_a = (g_0 - iw)*t0 + iw*(ta);
                    elseif solution == "con_7"
                        mu8 = 1/lambda *B;    
                        t_a = (g_0 - mu8)*t0 + mu8*ta;
                    end
                end

                %Unconstrained solutions and constained solutions V and
                %VIII
                if solution == "uncon"
                    t_a = inv(covar_orig)*ExpRet*(1/(B));
                    t_0 = inv(covar_orig)*one.*(1/(one'*inv(covar_orig)*one));
                end

                cos_phi = (ExpRet'*t_a)/(norm(ExpRet)*norm(t_a)); %use ExpRet here becasue you dont know it neces
                angle = acosd(cos_phi); % alpha-weight angle
                angles_orig(i,k,q) = angle;
 
                if q <6 % For Diag, CCM, OPM, TPM, OFMM shrinkages
                   covar_adj = ((z/100))*I_adj + (1-(z/100))*covar_orig;
                elseif q==6 % For G&J shrinkage
                   if angle> 90
                        new_angle = 89.99999999999;
                        new_cos_phi = cosd(new_angle);
                        covar_adj = ((z/100)/new_cos_phi)*I + (1-((z/100)/new_cos_phi))*covar_orig;
                 
                   else
                        covar_adj = ((z/100)/cos_phi)*I + (1-((z/100)/cos_phi))*covar_orig;
                   end
                elseif q==7 % For AOCS shrinkage
                    covar_adj = covarianceShrinkage(Data_R);
                end

                condition_o(i) = cond(covar_orig); %condition number sample covariance matrix
                condition_t(i) = cond(covar_adj); % condition number of shrunk covariance matrix

                % _________________________________________________________
                % Adjusted portfolio solution
                 
                A = one'*inv(covar_adj)*one;
                B = transpose(one)*inv(covar_adj)*ExpRet;
                C = ExpRet'*inv(covar_adj)*ExpRet;
                D = (A*C) - B^2;

                %Constrained solutions VI or VII
                if (solution == "con_6" || solution == "con_7")
                    w_a = inv(covar_adj)*ExpRet*(1/(transpose(one)*inv(covar_adj)*ExpRet));
                    w_0 = inv(covar_adj)*one.*(1/(one'*inv(covar_adj)*one));
                    if solution == "con_6"
                        iw = (B/D)*((alpha_0*A)-(g_0*B));
                        w1 = (g_0 - iw)*w_0 + iw*(w_a);
                    elseif solution == "con_7"
                        mu8 = 1/lambda *B;    
                        w1 = (g_0 - mu8)*w_0 + mu8*w_a;
                    end
                end

                %Unconstrained solutions and constained solutions V and
                %VIII
                w1 = inv(covar_adj)*ExpRet*(1/(transpose(one)*inv(covar_adj)*ExpRet)); % Optimal Risky Portfolio (Sigma_alpha)
                w0 = inv(covar_adj)*one.*(1/(one'*inv(covar_adj)*one)); % Minimum Variance Portfolio

                if rcond(covar_orig)<0.0000001 % Error Handling
                    return_adj_portfolio(i,k) = NaN;
                    break
                elseif round(sum(t_a),1)~=1
                    return_adj_portfolio(i,k) = NaN;
                    break
                else
         
                    cos_phi_adj = (exp_ret_fwd'*w1)/(norm(exp_ret_fwd)*norm(w1)); %can use exp_ret_fwd to show how close it is to actual rets
                    agl = acosd(cos_phi_adj); %actual alpha-weight angle
                    angles_improved(i,k,q) = acosd(cos_phi_adj);
       
                    return_adj_portfolio(i,k) = Data(i+(12*roll)+1,:)*w1*1; % return calculated for ORP
          
                    norm_w1(i,k)= norm(w1); % magnitude of weight vector
                  
                end
            end

            if rcond(covar_orig)<0.0000001 || round(sum(t_a),1)~=1 % Error Handling
                t=t-1; 
                break
            else     
               
               sharpe_adj(t,k,q) = sharpe_L(return_adj_portfolio(:,k)); % Sharpe Ratio        
               performance(t,k,q) = compound_returns(return_adj_portfolio(:,k)); % subsequent 5Y performance, compounded return   
               nrms(t,k,q) = mean(norm_w1(:,k)); % average weight vector magnitude
               risk(t,k,q) = std(return_adj_portfolio(:,k));  % std deviation of 5Y subsequent returns         
               prob_sr(t,k,q) = computePSR(return_adj_portfolio(:,k),0); % Probabilistic Sharpe Ratio
               cond_improve(t,k,q) = mean(condition_o-condition_t)/mean(condition_o); %average improvement in matrix condition number
               angles_orig_ave(t,k,q) = median(angles_orig(:,k,q)); 
               angles_improved_ave(t,k,q) = median(angles_improved(:,k,q)); % average improvement in alpha-weight angle
            end

            k=k+1;
        end
        
        % Calculation of deflated sharpe ratio, to assess optimal shrinkage
        % intensity
        if rcond(covar_orig)<0.0000001 || round(sum(t_a),1)~=1
        else
            compare_returns = [return_adj_portfolio]; 
            compare_sharpe = [sharpe_adj(:,:,q)];
            N = size(compare_returns,2);
            T = size(compare_returns,1);
            var_Sharpe = var(compare_sharpe(t,:));
            skews = skewness(compare_returns);
            kurts = kurtosis(compare_returns);
            for i=1:size(compare_returns,2)
                [deflatedSharpe(i,1), pvalue(i,1)] = dsr(N, T, compare_sharpe(t,i) , var_Sharpe, skews(1,i), kurts(1,i));
                [max_dsr(t,q), index] = max(deflatedSharpe);
                which_pvalue(t,q) = index;
                max_pval(t,q) = pvalue(index,1);
            end
        end
        t = t+1;
    end

end