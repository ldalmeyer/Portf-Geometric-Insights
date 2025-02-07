# Portf-Geometric-Insights
Code for Project: Linearly constrained robust portfolio construction using Bauer-Householder bounds on the condition number

Matlab: Simulation and Result Generation

The simulation and computational work of this thesis was performed in Matlab R2023a, using an academic license.

The only calculations performed outside of Matlab, was to generate figures 7.14 and 7.18; this was performed in Microsoft excel. The excel spreadsheet and figures are provided.

Below, we describe each of the scripts provided. We make reference to where this script was utilised in the project, and any instructions relevent to run the script and reproduce results.

**Code: Simulated Pareto Surfaces**
The script "Eff_Frontier_Surf.m" generates figure 7.1, in section 7.2. 

Code Function : The code generates figures which visualizes and compares the efficient frontiers for optimizations problems I, II, III and VI in Table 1.1, with those that have explicit gearing constraints, namely optimisations V, VI, VII and VIII from table 4.1.

Instructions : The code, when run, will generate figures automatically. The dataset Index20052016 is utilised for its expected return (α) and covariance inputs (Σ) illustratively (note, any asset return data would suffice). It is monthly total return data from 2005 to 2016, of South African cash, bonds, and eight different equity sectors.

**Code: Simulated alpha-angle dependency**
The script "ND_simulations_charts.m" generates the simulations performed in section 7.2, and creates figure 7.2.

Code Function : This script simulates 10 years worth of monthly stock return data, for 10 stock portfolios, from normal distributions with random mean (between 0% and 50% p.a) and random standard deviation (between 0% and 20% p.a.) with 1000 replications per simulation. The α-weight angle for the equally weighted portfolio is calculated and compared to the α-weight angle of the optimal risk portfolio on every simulation. The script then generates two figures that show the percentage number of performance simulations with weights better aligned with those of the returns, and the average difference between the conditioned angles and the unconditioned angles with confidence intervals. One figure is for the unconstrained optimisations (and constrained solutions V and VIII), and the other for the constrained cases.

Instructions : The code will run and generate figures automatically, except for one file that needs to be downloaded by the user:

1. The function degreetick, developed by Green (2025), is used to adds degree symbols to tick labels on an axis.

Parameters are set as was used in the thesis

1. lambda: This is the risk aversion parameter required for constrained solution VII (see Table 4.1).
2. alpha_0 : This is the minimum return parameter required for constrained solution VI (see Table 4.1).

**Code: Exploratory Data Analysis**
The script "EDA.m" performs the exploratory data analysis performed in subsection 7.3.0.1, and generates Figures 7.3 to 7.10.

The data supplied by Peresec for the simulation work is not published, but can be easily recreated to anyone with access to a Bloomberg terminal. The 194 stocks that have a minimum 10 year performance history between January 1996 to August 2016 is listed in Section C of the project. The monthly total return data of the 194 stocks are used, and retrieved as described in Section 7.3.0.1 in Bloomberg. We also describe there how the data is adjusted for dividends and corporate actions within Bloomberg. For the purpose of testing the script, sample data was provided in order to allow the user to run the code. The sample data is called sample data.mat. From the original Peresec data, a random start date was selected, and the following 210 months was selected, for 50 random stocks. Note that all 50 stocks have full data availability over the selected time period. While the results will differ from those presented in this thesis, the gist of the results can be noted from this smaller dataset.

Code Function : The script generates the following charts, in performing exploratory data analysis:
1. A histogram that demonstrates the available history of the considered stocks, over the available dates (figure 7.3).
2. A scatter plot showing the median, 25th and 75th percentiles, and minimum and maximum values of the data (excluding outliers) (figure 7.4).
3. A scatter plot showing the median, 25th and 75th percentiles, the minimum and maximum values of the data, and the average of the outliers per counter (figure 7.5).
4. Histogram displaying the number of outliers (figure 7.6).
5. Two boxplots, one demonstrating the number of inter-quartile ranges the outliers are above the 75th percentile of the data (dispersion of upper quartiles), and the other demonstrating the number of inter-quartile ranges the outliers are below the 25th percentile of the data (dispersion of lower quartiles) (figure 7.7).
6. Four plots, each a histogram showing the distribution of the first four moments of the data (the mean, standard deviation, skewness and kurtosis) (figures 7.8 and 7.9).
7. A boxplot showing the variation of the correlations present within the data (figure 7.10a).
8. A histogram showing the frequency of the correlation within the data (figure 7.10b).
   
Instructions : The code will run and generate figures automatically; no parameter specification is required.

**Code: Simulation Shrinkage Results**
The script "real_world_simulations_JSE.m" performs the simulation work that is performed on historical total return stock data, described in seen in section 7.3. The figures illustrating these results are performed in script JSE sim charts.m. These simulations are performed for unconstrained solutions (see Table 1.1) and constrained solutions V and VII (see Table 4.1). More particularly, this code generates the results seen in:
1. Subsection 7.3.2; the results are illustrated in Charts 7.11, 7.12 and 7.13.
2. Subsection 7.3.2.1; the results are illustrated in Charts 7.15, 7.16 and 7.17.
3. Subsection 7.3.2.2; the results are illustrated in Charts 7.19, 7.20 and 7.21.
4. Subsection 7.3.4; the results are illustrated in Charts 7.22, 7.23, 7.24 and 7.25

Code Function : The script compares the Golts and Jones (2009) methodology against other shrinkage methodologies, by performing a set of 1000 simulations, and assessing the out of sample results of each method. Each simulation follows the steps outlined in Table 7.1. The shrinkage methodologies investigated in this script are those developed by Ledoit (1995) and Ledoit and Wolf (2003b, 2004b). The seven different methodologies we consider are provided in section 7.3.1 in detail.

Instructions : When running the code to produce results seen in this these, the following parameters are described, and need to be set as required to reproduce simulation results.
1. for ExpRet: This variable defines how your expected return vector is calculated. This should be set to 1 when expected returns are perfectly foreseeable, and forward returns are used as the expected return vector (this generates results seen in subsection 7.3.4), and the variable should be set to 0 when one has no insight into expected returns and historical average returns are used (this generates results seen in subsection 7.3.2).
2. solution: This variable specifies which portfolio optimisation solution you are solving for, and should be set to uncon when one wants to solve for unconstrained solutions (see Table 1.1) or constrained solution V and VII , to con 6 when one wants to solve for constrained solution VI, and to con 7 when one wants to solve for constrained solution VII (see Table 4.1).

Other parameters the user is able to set at the start of the script are as follows, but these are set to parameters as used in the thesis and are not required to be changed:
1. tp: The total performance period over which the methodology is assessed; this paper has elected 10 year periods. This paper calculates 5 years worth of out-of-sample performance data, using 5 years historical rolling data.
2. roll : The number of years to use as the historical data, off which to calculated the sample covariance matrix and expected returns (assuming no foreseeability of expected returns); this paper has elected 5 year periods.
3. nsim : The number of simulations to perform; this paper elects 1000 simulations.
4. ns: The number of stocks to include in the portfolio; this paper chose 10 random stocks.

Users need to download other scripts, that are made reference to in real world simulations JSE.m.
These scripts are as follows:
1. Code adjusting the sample covariance matrix for the different shrinkage methodologies was sourced from the CovShrinkage package, developed by Ledoit. The Github link is provided in the reference. This package provides the code for the diagonal matrix shrinkage methodology (Diag), constant correlation shrinkage methodology (CCM), one parameter matrix shrinkage methodology (OPM), two parameter ma- trix shrinkage methodology (TPM) and the one factor model shrinkage methodology (OFMM). Matlab function covarianceShrinkage, part of the financial toolbox, is used to calculate the asymptotically optimal convex shrunk (AOCS) covariance matrix.
2. The function "Sharpe_L.m" is used to calculate standard Sharpe ratios, and was developed for use in this thesis, provided in the code and data depository.
3. The function "compound_returns.m" is used to calculate compounded returns, and was developed for use in this thesis, provided in the code and data depository.
4. The function "computePSR.m" is used to calculate probabilistic sharpe ratios, and was developed for use in this thesis, provided in the code and data depository.
5. The function "dsr.m" is used to calculate deflated sharpe ratios, and was developed for use in this thesis, provided in the code and data depository. The methodology for this code was sourced from Matare, who developed this in R.

Once script "real_world_simulations_JSE.m" has been run, with the appropriate parameters set, the charts can be generated from "JSE_sim_charts.m". This script is broken up into 4 sections:
1. Section 1 "Results- Group of 4 Charts": This generates the Sharpe ratio, probabilistic Sharpe ratio, average α-weight angle and covariance matrix condition number improvement charts. More specifically, figures 7.11, 7.12, 7.13, 7.23 and 7.25.
2. Section 2 "Intensity Factor (k)": This generates the boxplots of deflated sharpe ratios and associated p-values, showing which intensity factor is optimal. More specifically this generates figures 7.15, 7.16 and 7.17.
3. Section 3 "Weight Vector Magnitude": This generates the line charts showing the average magnitude of the weight vectors. More specifically, this generates figure 7.22.
4. Section 4 "Performance and Risk at chosen k" : This generates the boxplots showing the 5 year performance and risk metric of each shrinkage methodology over the simulations, at a specified shrinkage intensity. More specifically, this generates figures 7.19, 7.20 and 7.21.

**References**

Golts, M., Jones, G.C., 2009. A sharper angle on optimization. SSRN (October 5, 2009). URL: https://ssrn.com/abstract=1483412, doi:http://dx.doi.org/10.2139/ssrn.1483412.

Ledoit, O., covShrinkage. URL: https://github.com/oledoit/covShrinkage/releases/tag/1.1.0.

Ledoit, O., 1995. Essays on risk and return in the stock market. PhD thesis, Massachusetts Institute of Technology, Sloan School of Management. 88

Ledoit, O., Wolf, M., 2003b. Improved estimation of the covariance matrix of stock returns with an application to portfolio selection. Journal of Empirical Finance 10, 603–621.

Ledoit, O., Wolf, M., 2004b. A well-conditioned estimator for large-dimensional covariance matrices. Journal of Multivariate Analysis 88, 365–411.

Matare, N., dsr: Deflated sharpe ratio. URL: https://rdrr.io/github/nmatare/quanttools/man/dsr.html.

Chad Greene (2025). degreetick (https://www.mathworks.com/matlabcentral/fileexchange/50792-degreetick), MATLAB Central File Exchange. Retrieved February 2, 2025.
