function M4_main_001_17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENGR 132 
% Program Description 
% This program load the data from PGO-X50 to find the parameters using
% M4_Algorithm_001_17.m file. Plot the Michaelis-Menten model and initial
% velocity tangent line with original data. Find the SSE to use it for the
% error calculation process
%
% Function Call
% M4_main_001_17
%
% Input Arguments
% None
%
% Output Arguments
% None
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;
close all;
%% ____________________
%% INITIALIZATION

%PGOX50
data_pgo = readmatrix('Data_PGOX50_enzyme');% importing the data

s_initial = data_pgo(2,2:11); % initial substrate concentration (uM)
time_pgo = data_pgo(5:end,1); % time (s)
subs_pgo = data_pgo(5:end,2:end); % product substrate (uM) 

% loop variable
% test1_3_size = 3; % column size for test 1 through 3
% test4_6_size = 3; % column size for test 4 through 6
% test7_10_size = 4; % column size for test 7 through 10

% Next Generation Enzyme
data_enz_raw = readmatrix('Data_nextGen_KEtesting_allresults'); % importing data
s_initial_enz = data_enz_raw(1,2:11); % iniitial substrate concentration (uM)
time_enz = data_enz_raw(3:end,1); % time (s)
data_enz = data_enz_raw(3:end,2:end); % matrix that only contains enzymes (uM)

% pricing data
data_price = readmatrix('Data_NaturalCatalysts_priceCatalog');
%% ____________________
%% CALCULATIONS
%calling algorithm function PGOX50
[V_init,vmax_pgo,km_pgo] = M4_Algorithm_001_17(subs_pgo,time_pgo,s_initial);

v0_tangent = V_init.*time_pgo; % initial velocity tangent line vector (uM/s)
v_mich = (vmax_pgo.*s_initial)./(km_pgo + s_initial); % velocity vector using Michaelis-Menten model (uM/s)
y_lim_v0= (s_initial_enz+s_initial_enz.*0.2); % proper window size to display initial velocity tangent line for y axis

%finding SSE for PGOX50 Michaelis Menten Model
sse_pgo = sum((V_init - v_mich).^2);




% finding parameters for next generation enzymes

ctr = 1; % setting the loop controll variable
V_init_total = zeros(10,10); % pre-allocating vector with zeros for total initial velocity vector
vmax_enz = zeros(1,10); % pre-allocating vector with zeros for vamx
km_enz = zeros(1,10); % pre-allocating vector with zeros for km
sse_enz = zeros(1,10); % pre-allocating vector with zerose for SSE
v_mich_enz = zeros(10,10); % pre- allocating vector with zeros for velocity founded by using model 

% loop to find the parameter for Enzyme with duplicates 
for ind = 10:10:100 % counting from enzyme A up to enzyme E duplicate

    [V_init_enz,vmax_enz(ctr),km_enz(ctr)] = M4_Algorithm_001_17(data_enz(:,ind-9:ind),time_enz,s_initial_enz); % calling functions
    V_init_total(ctr,:) = V_init_enz; % plugging in vector into matrix
    v_mich_enz(ctr,:) = (vmax_enz(ctr).*s_initial_enz)./(km_enz(ctr) + s_initial_enz); % velocity vector using Michaelis-Menten model (uM/s)
    sse_enz(ctr)= sum((V_init_total(ctr,:) - v_mich_enz(ctr,:)).^2); % finding sse
    ctr = ctr + 1; % increasing ctr by 1 for each loop
end

%% improvement 2
% Selecting less SSE enzyme between duplicate and original 
less_sse = zeros(1,5); % pre-allocating vectors for less sse enzymes
ctr2= 1; % loop controll variable
for ind2= 2:2:10 % going through all the sse
    less_sse(ctr2) = min(sse_enz(ind2-1),sse_enz(ind2)); % finding less sse test
    ctr2 = ctr2 + 1; % increase the controll variable
end



%% ____________________
%% FORMATTED TEXT/FIGURE DISPLAYS
% Plots with initial velocity tangent line lies on original data using for
% loop to display to display in two different figure (We commented plot
% since we don't need it for M4)

% Plots with initial velocity tangent line lies on original data using for
% loop to display to display in two different figure

% % test 1 through 3
% for ind1 = 1:1:test1_3_size
%     
%     figure(1)
%     subplot(1,test1_3_size,ind1);
%     plot(time_pgo,subs_pgo(:,ind1),'k');
%     hold on
%     plot(time_pgo,v0_tangent(:,ind1),'r--');
%     xlim([0 inf]);
%     ylim([0 y_lim_v0(ind1)]); % adjusting index size to get proper y-axis limit
%     grid on 
%     title(['Initial substrate concentration ' num2str(s_initial(ind1)) 'uM']);% using built in num2str function to display initial substrate concentration
%     xlabel('Time (s)');
%     ylabel('Product concentration (uM)');
%     legend('data','initial velocity','Location','best');
%     sgtitle('The original data with initial velocity tangent line of PGOX50 test 1 - 3');
%   
% end
% 
% 
% % test 4 through 6
% for ind2 = 1:1:test4_6_size
%     
%     figure(2)
%     subplot(1,test4_6_size,ind2);
%     plot(time_pgo,subs_pgo(:,ind2+3),'k'); % adjusting index size to get proper graph
%     hold on
%     plot(time_pgo,v0_tangent(:,ind2+3),'r--'); % adjusting index size to get proper graph
%     xlim([0 inf]);
%     ylim([0 y_lim_v0(ind2+3)]); % adjusting index size to get proper y-axis limit
%     grid on 
%     title(['Initial substrate concentration ' num2str(s_initial(ind2+3)) 'uM']);% using built in num2str function to display initial substrate concentration
%     xlabel('Time (s)');
%     ylabel('Product concentration (uM)');
%     legend('data','initial velocity','Location','best');
%     sgtitle('The original data with initial velocity tangent line of PGOX50 test 4 - 6');
%   
% end
% 
% % test 7 through 10
% for ind3 = 1:1:test7_10_size
%     
%     figure(3)
%     subplot(1,test7_10_size,ind3);
%     plot(time_pgo,subs_pgo(:,ind3+6),'k'); % adjusting index size to get proper graph
%     hold on
%     plot(time_pgo,v0_tangent(:,ind3+6),'r--'); % adjusting index size to get proper graph
%     xlim([0 inf]);
%     ylim([0 y_lim_v0(ind3+6)]); % adjusting index size to get proper y-axis limit
%     grid on 
%     title(['Initial substrate concentration ' num2str(s_initial(ind3+6)) 'uM']);% using built in num2str function to display initial substrate concentration
%     xlabel('Time (s)');
%     ylabel('Product concentration (uM)');
%     legend('data','initial velocity','Location','best');
%     sgtitle('The original data with initial velocity tangent line of PGOX50 test 7 - 10');
%   
% end
%     





% Text display PGOX50
fprintf('\n------------------------------------------\n')
pg_txt = 'PGOX50 has vmax:%.3f km:%.3f sse: %.4f\n';
fprintf(pg_txt,vmax_pgo,km_pgo,sse_pgo);


% Printing average value of km ,vmax , and sse for each enzymes
% Enzyme A
enz_a_txt ='Enzyme A has vmax:%.3f km:%.3f sse: %.4f\n';
loc_A = find(less_sse(1)== sse_enz); % finding less SSE enzyme between original and duplicate
fprintf(enz_a_txt,vmax_enz(loc_A), km_enz(loc_A),sse_enz(loc_A));

% Enzyme B
loc_B = find(less_sse(2)== sse_enz);% finding less SSE enzyme between original and duplicate
enz_b_txt ='Enzyme B has vmax:%.3f km:%.3f sse: %.4f\n';
fprintf(enz_b_txt,vmax_enz(loc_B), km_enz(loc_B),sse_enz(loc_B));

% Enzyme C
loc_C = find(less_sse(3)== sse_enz);% finding less SSE enzyme between original and duplicate
enz_c_txt = 'Enzyme C has vmax:%.3f km:%.3f sse: %.4f\n';
fprintf(enz_c_txt,vmax_enz(loc_C), km_enz(loc_C),sse_enz(loc_C));

% Enzyme D
loc_D = find(less_sse(4)== sse_enz);% finding less SSE enzyme between original and duplicate
enz_d_txt = 'Enzyme D has vmax:%.3f km:%.3f sse: %.4f\n';
fprintf(enz_d_txt,vmax_enz(loc_D), km_enz(loc_D),sse_enz(loc_D));


% Enzyme E
loc_E = find(less_sse(5)== sse_enz);% finding less SSE enzyme between original and duplicate
enz_e_txt = 'Enzyme E has vmax:%.3f km:%.3f sse: %.4f\n';
fprintf(enz_e_txt,vmax_enz(loc_E), km_enz(loc_E),sse_enz(loc_E));

% pricing non-linear regression
[sse_price,sst_price,r_price]=M4_Regression_001_17(data_price,km_enz(loc_A),km_enz(loc_B),km_enz(loc_C),km_enz(loc_D),km_enz(loc_E));

% % plot the Michaelis Menten Model
figure(5)
plot(s_initial_enz, v_mich_enz(loc_A,:)); % Enzyme A
hold on 
grid on
plot(s_initial_enz, v_mich_enz(loc_B,:)); % Enzyme B
plot(s_initial_enz, v_mich_enz(loc_C,:)); % Enzyme C
plot(s_initial_enz, v_mich_enz(loc_D,:)); % Enzyme D
plot(s_initial_enz, v_mich_enz(loc_E,:)); % Enzyme E
legend('enzyme A','enzyme B','enzyme C','enzyme D','enzyme E','Location','best');
xlabel('substrate concentration (uM)');
ylabel('initial velocity (uM/s)');
title('Michaelis Menten Kinetics for Next Generation Enzymes');


    

%% ____________________
%% RESULTS
% M3
% PGOX50
% The initial velocity (uM/s) vector is
%     0.0264    0.0509    0.1081    0.1842    0.3288    0.6038    0.9400    1.1719    1.2925    1.5465
% 
% km: 240.431 (uM) vmax: 1.715 (uM/s)
% The SSE for this model is 0.0145

% M4
% ------------------------------------------
% PGOX50 has vmax:1.718 km:247.266 sse: 0.0110
% Enzyme A has vmax:0.892 km:146.373 sse: 0.0040
% Enzyme B has vmax:0.918 km:370.487 sse: 0.0040
% Enzyme C has vmax:1.209 km:200.077 sse: 0.0002
% Enzyme D has vmax:1.481 km:279.174 sse: 0.0159
% Enzyme E has vmax:1.523 km:158.097 sse: 0.0190

%% ____________________
%% ACADEMIC INTEGRITY STATEMENT
% We have not used source code obtained from any other unauthorized
% source, either modified or unmodified. Neither have we provided
% access to my code to another. The program we are submitting
% is our own original work.




