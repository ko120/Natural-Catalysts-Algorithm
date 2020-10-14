function [sse_price, sst_price, r_price]=M4_Regression_001_17(data_price,enz_a,enz_b,enz_c,enz_d,enz_e)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENGR 132 
% Program Description 
% This program find the regression model for NaturalCatalysts' price data
% with considering enzyme performance. Display the linear regression
% equation, plot the Michaelis constant for all enzymes, and plot the
% regression model. Lastly, find the SSE, SST, and r^2 for the linearized
% data.
%
% Function Call
% [sse_price,sst_price,r_price]=M4_Regression_001_17(data_price,enz_a,enz_b,enz_c,enz_d,enz_e_)
%
% Input Arguments
% data_price: The price data excel that includes km and price
% enz_a: km value for enzyme A(uM)
% enz_b: km value for enzyme B(uM)
% enz_c: km value for enzyme C(uM)
% enz_d: km value for enzyme D(uM)
% enz_e: km value for enzyme E(uM)
%
% Output Arguments
% sse_price: Error Sum of Squares of the linear regression model
% sst_price: Total sum of squares of the linear regression model
% r_price: Coefficient of determination of the linear regression model
%
% Assignment Information
%   Assignment:     M4, Problem part 3
%   Team member:    John Groves,  jsgroves@purdue.edu
%                   Caleb Hancock hancoc25@purdue.edu
%                   Danny Kim     kim3255@purdue.edu  
%                   Kyung Min Ko, ko120@purdue.edu
%   Team ID:        001-17
%   Academic Integrity:
%     [X] We worked with one or more peers but our collaboration
%        maintained academic integrity.
%     Peers we worked with: John Groves,  jsgroves@purdue.edu
%                           Caleb Hancock hancoc25@purdue.edu
%                           Danny Kim     kim3255@purdue.edu  
%                           Kyung Min Ko, ko120@purdue.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ____________________
%% INITIALIZATION

km_data = data_price(:,1); % Michaelis Constant (uM)
price_data = data_price(:,2); % price (USD $ dollar)


%% ____________________
%% CALCULATIONS

%% Linearizing data
price_log = log10(price_data); % log scailling y axis
coeff_lin = polyfit(km_data,price_log,1); % finding coefficient for linear model
predicted_price = polyval(coeff_lin,km_data); % finding predicted price using linear model


%% Generalized model
slope_exp = coeff_lin(1); % slope of the exponential function 
y_intercept = 10^coeff_lin(2); % taking anti-log for the y intercept to find the general equation
predicted_price_gen = y_intercept.*10.^(slope_exp.*km_data); % finding predicted price using general exponentail equation


%% finding SSE SST r^2

sse_price = sum((price_log - predicted_price).^2); % sse using linearized model
sst_price = sum((price_log - mean(price_log)).^2); % sst ussing linearized model
r_price = 1- sse_price/sst_price; % r^2 using linearized model


%% finding price for each enzyme
enz_a_price = y_intercept.*10.^(slope_exp.*enz_a); % enzyme A
enz_b_price = y_intercept.*10.^(slope_exp.*enz_b); % enzyme B
enz_c_price = y_intercept.*10.^(slope_exp.*enz_c); % enzyme C
enz_d_price = y_intercept.*10.^(slope_exp.*enz_d); % enzyme D
enz_e_price = y_intercept.*10.^(slope_exp.*enz_e); % enzyme E


%% ____________________
%% FORMATTED TEXT/FIGURE DISPLAYS

%% plotting in different scales

figure(1)
% Linear scaled
subplot(2,2,1)
plot(km_data,price_data,'b*');
grid on
title('Linear scaled plot')
xlabel('Michaelis constant km (uM)');
ylabel('price (dollar per lb)');
sgtitle('Relation between price and Km in various scaled axis')

% Log scaled
subplot(2,2,2)
semilogx(km_data,price_data,'b*');
grid on
title('Log scaled plot');
xlabel('Michaelis constant km (uM)');
ylabel('price (dollar per lb)');

% Exponential scaled
subplot(2,2,3)
semilogy(km_data,price_data,'b*');
grid on
title('Exponential scaled plot');
xlabel('Michaelis constant km (uM)');
ylabel('price (dollar per lb)');

% Power scaled
subplot(2,2,4)
loglog(km_data,price_data,'b*');
grid on
title('power scaled plot');
xlabel('Michaelis constant km (uM)');
ylabel('price (dollar per lb)');

%% Ploting Regression models

% Linear equation model
figure(2)
plot(km_data, price_log,'*r');
hold on
plot(km_data, predicted_price,'k');
grid on
xlabel('Michaelis constant km (uM)');
ylabel('Log scaled price (dollar per lb)');
title('Linear equation model for relation between price and km');
legend('data','model','Location','best');
caption1 = 'log(price)= -0.0045 * km + 3.2902';
text(170, 1.8 ,caption1, 'FontSize', 11, 'Color', 'k');

% General equation model
figure(3)
plot(km_data,price_data,'*r');
hold on
plot(km_data,predicted_price_gen,'k');
grid on
xlabel('Michaelis constant km (uM)');
ylabel('price (dollar per lb)');
title('General exponential equation model for relation between price and km');
legend('data','model','Location','best');
caption2 = 'price= 1950.9224 * 10^{-0.0045 * km}';
text(220, 350 ,caption2, 'FontSize', 11, 'Color', 'k');

% Equation display
lin_txt = 'The linear model euqation is log(price)= %.4f * km + %.4f\n';
fprintf(lin_txt, coeff_lin(1),coeff_lin(2));
exp_txt = 'The general exponential equation is price= %.4f * 10^(%.4f * km)\n';
fprintf(exp_txt, y_intercept,slope_exp);
stat_txt ='The sse is %.4f sst is %.4f and r^2 is %.4f\n';
fprintf(stat_txt, sse_price, sst_price, r_price); 
pr_txt = 'The price of enzyme A is %.2f $ enzyme B is %.2f $ enzyme C is %.2f $ enzyme D is %.2f $ enzyme E is $%.2f\n';
fprintf(pr_txt,enz_a_price,enz_b_price,enz_c_price,enz_d_price,enz_e_price);
%% ____________________
%% RESULTS
% The linear model euqation is log(price)= -0.0045 * km + 3.2902
% The general exponential equation is price= 1950.9224 * 10^(-0.0045 * km)
% The sse is 0.1367 sst is 4.4451 and r^2 is 0.9693
% The price of enzyme A is 432.08 $ enzyme B is 42.97 $ enzyme C is 248.52 $ enzyme D is 110.05 $ enzyme E is $382.93

%% ____________________
%% ACADEMIC INTEGRITY STATEMENT
% We have not used source code obtained from any other unauthorized
% source, either modified or unmodified. Neither have we provided
% access to my code to another. The program we are submitting
% is our own original work.



