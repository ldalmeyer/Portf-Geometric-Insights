 clearvars -except ave_improvement_in_angle_B_Un ave_improvement_in_angle_U_Un ave_improvement_in_angle_Un perc_sims_AWA_improved_Un ...
    ave_improvement_in_angle_B_VI ave_improvement_in_angle_U_VI ave_improvement_in_angle_VI perc_sims_AWA_improved_VI...
    ave_improvement_in_angle_B_VII ave_improvement_in_angle_U_VII ave_improvement_in_angle_VII perc_sims_AWA_improved_VII 
    

for a = 1:3 % a=1 is unconstrained solution, a=2 is constrained solution 6, a=3 is constrained solution 7

    clc
    
    if a==2
        alpha_0 = 0.05; % minimum ER of portfolio
    else    
        g_0 = 1; % sum of weights
    end
    
    k = logspace(-2,2,101);
    
    for t=1:101
    
    i=1;
    
    while i < 1001
    %% Generate Stock Returns
    
    for j=1:10
    randmean = randi(50)/(100*12); %random mean between 0% and 50% per year
    randstd = randi(20)/(100*sqrt(12)); %random std between 0% and 20% per year
    Returns(:,j) = normrnd(randmean,randstd,120,1);
    end
    
    ERet = mean(Returns)';
    covar = cov(Returns);
    one = ones(size(Returns,2),1);
    
    %% Calc Alpha W Angle
    
    wd = inv(covar)*ERet.*(1/(transpose(ERet)*inv(covar)*one)); %change to wd unconstrained
    w0 = inv(covar)*one.*(1/(one'*inv(covar)*one));
    
    if a>1
    
        w_a = wd;
        w_0 = w0;
        
        A = one'*inv(covar)*one;
        B = transpose(one)*inv(covar)*ERet;
        C = ERet'*inv(covar)*ERet;
        D = (A*C) - B^2;
                            
        if a==2 % Sln VI
            iw = (B/D)*((alpha_0*A)-(g_0*B));
            wd = (g_0 - iw)*w_0 + iw*(w_a);
        
        elseif a==3 % Sln VII
            lambda = 1; % risk aversion parameter
            mu8 = 1/lambda *B;                    
            wd = (g_0 - mu8)*w_0 + mu8*w_a;
        end
    end
   
    ew = repmat(1/10,10,1);
    
    A_W(i,t) = acosd((ERet'*wd)/(sqrt(ERet'*ERet)*sqrt(wd'*wd)));
    cosA_W(i,t) = (ERet'*wd)/(sqrt(ERet'*ERet)*sqrt(wd'*wd));
    
    A_EW(i,t) = acosd((ERet'*ew)/(sqrt(ERet'*ERet)*sqrt(ew'*ew)));
    cosA_EW(i,t) = (ERet'*ew)/(sqrt(ERet'*ERet)*sqrt(ew'*ew));
    
    %z=t-1;
    I = eye(size(Returns,2));
    B = transpose(one)*inv(covar)*ERet;
    cos_phi = (ERet'*wd)/(norm(ERet)*norm(wd));                   
    covar_adj2 = ((k(t)/100)/cos_phi)*I + (1-((k(t)/100)/cos_phi))*covar;
    wd2 = inv(covar_adj2)*ERet.*(1/(transpose(one)*inv(covar_adj2)*ERet));
    
        if isnan(sum(wd2))
        xxx=1;
        else
        A_AW(i,t) = acosd((ERet'*wd2)/(sqrt(ERet'*ERet)*sqrt(wd2'*wd2)));
        cosA_AW(i,t) = (ERet'*wd2)/(sqrt(ERet'*ERet)*sqrt(wd2'*wd2));
        
        TrueParam(i,t) = A_AW(i,t)<A_EW(i,t);
        degree_improvement(i,t) = A_EW(i,t)-A_AW(i,t);
        i=i+1; 
        end 
    end
    
    
    end
    
    perc_sims_AWA_improved = mean(TrueParam)*100;
    ave_improvement_in_angle = nanmean(degree_improvement,1);
    stdTP = std(ave_improvement_in_angle);
    ave_improvement_in_angle_U = ave_improvement_in_angle+stdTP;
    ave_improvement_in_angle_B = ave_improvement_in_angle-stdTP;
    
    if a==1
        perc_sims_AWA_improved_Un = perc_sims_AWA_improved;
        ave_improvement_in_angle_Un = ave_improvement_in_angle;
        ave_improvement_in_angle_U_Un = ave_improvement_in_angle_U;
        ave_improvement_in_angle_B_Un = ave_improvement_in_angle_B;
    elseif a==2
        perc_sims_AWA_improved_VI = perc_sims_AWA_improved;
        ave_improvement_in_angle_VI = ave_improvement_in_angle;
        ave_improvement_in_angle_U_VI = ave_improvement_in_angle_U;
        ave_improvement_in_angle_B_VI = ave_improvement_in_angle_B;
    elseif a==3
        perc_sims_AWA_improved_VII = perc_sims_AWA_improved;
        ave_improvement_in_angle_VII = ave_improvement_in_angle;
        ave_improvement_in_angle_U_VII = ave_improvement_in_angle_U;
        ave_improvement_in_angle_B_VII = ave_improvement_in_angle_B;
    end    

end




%% Charts


figure()
%subplot(1,2,1)
axis square
yyaxis left
semilogx(k,(perc_sims_AWA_improved_Un),'o','MarkerFaceColor', 'b')
%ylim([0 100])
ytickformat('%g%%')
ax = gca; % Get handle to current axes.
ax.XColor = 'k'; % Red
ax.YColor = 'k'; % Blue
ylabel("$n_{\phi_{\omega*} < \phi_{1/N}} [\%]$", 'Interpreter','latex', 'FontSize',20)
xlabel("k")
yyaxis right
semilogx(k,(ave_improvement_in_angle_Un),'^','MarkerFaceColor', 'm','MarkerEdgeColor', 'm', 'Color','m')
hold on
semilogx(k,(ave_improvement_in_angle_U_Un),'--','MarkerFaceColor', 'k','MarkerEdgeColor', 'm', 'Color','m')
semilogx(k,(ave_improvement_in_angle_B_Un),'--','MarkerFaceColor', 'k','MarkerEdgeColor', 'm', 'Color','m')
ylim([-20 50]);
degreetick 'y'
ax.XColor = 'k'; % Red
ax.YColor = 'k'; % Blue
ylabel("$\mu_{\Delta\phi_k} \pm \sigma_{\Delta\phi_k} [{\circ}]$", 'Interpreter','latex','FontSize',20)
xlabel("log k", 'FontSize',20)
title("Unconstrained Portfolio Simulations",'Interpreter','latex','FontSize',18)

%subplot(1,2,2)
figure()
axis square
yyaxis left
semilogx(k,(perc_sims_AWA_improved_VI),'o','MarkerFaceColor', 'b')
%ylim([0 100])
ytickformat('%g%%')
ax = gca; % Get handle to current axes.
ax.XColor = 'k'; % Red
ax.YColor = 'k'; % Blue
ylabel("$n_{\varphi_{\omega*} < \varphi_{1/N}} [\%]$", 'Interpreter','latex', 'FontSize',20)
xlabel("k")
yyaxis right
semilogx(k,(ave_improvement_in_angle_VI),'^','MarkerFaceColor', 'm','MarkerEdgeColor', 'm', 'Color','m')
hold on
semilogx(k,(ave_improvement_in_angle_U_VI),'--','MarkerFaceColor', 'k','MarkerEdgeColor', 'm', 'Color','m')
semilogx(k,(ave_improvement_in_angle_B_VI),'--','MarkerFaceColor', 'k','MarkerEdgeColor', 'm', 'Color','m')
ylim([-20 50]);
degreetick 'y'
ax.XColor = 'k'; % Red
ax.YColor = 'k'; % Blue
ylabel("$\mu_{\Delta\varphi_k} \pm \sigma_{\Delta\varphi_k} [{\circ}]$", 'Interpreter','latex','FontSize',20)
xlabel("log k", 'FontSize',20)
title("Constrained Portfolio Simulations",'Interpreter','latex','FontSize',18)

hold on
yyaxis left
semilogx(k,(perc_sims_AWA_improved_VII),'o')
yyaxis right
semilogx(k,(ave_improvement_in_angle_VII),'^','MarkerEdgeColor', 'm', 'Color','m')
semilogx(k,(ave_improvement_in_angle_U_VII),'--','MarkerEdgeColor', 'm', 'Color','m')
semilogx(k,(ave_improvement_in_angle_B_VII),'--','MarkerEdgeColor', 'm', 'Color','m')
ylim([-20 50]);

qw{1} = plot(nan, 'bo','MarkerFaceColor', 'b');
qw{2} = plot(nan, 'bo');
qw{3} = plot(nan, 'm^','MarkerFaceColor', 'm');
qw{4} = plot(nan, 'm^'); % You can add an extra element too
qw{5} = plot(nan, 'm--'); % You can add an extra element too
%legend([qw{:}], {'Unconstrained Sln (count; left axis) [%]','Unconstrained Sln (diff angle; right axis) [^{\circ}]','\pm 1 std deviation'}, 'location', 'best')
legend([qw{:}], {'Soln VI (count; left axis) [%]','Soln VII (count; left axis) [%]','Soln VI (diff angle; right axis) [^{\circ}]', 'Soln VII (diff angle; right axis) [^{\circ}]','\pm 1 std deviation'}, 'location', 'best')

