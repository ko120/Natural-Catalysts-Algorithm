function [V_init,vmax,km] = M4_Algorithm_001_17(enz,time_enz,s_initial)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENGR 132 
% Program Description 
% This program smooth the original curve by using taking average of adjacent elements of vector, then find 
% the initial velocity, vmax, and km using linearized model for Michaelis-Menten model
% with using given data that includes each initial substrate concentration
% as time changing. 
%
% Function Call
% [V_init,vmax,km] = M4_Algorithm_001_17(enz,time_enz,s_initial)
%
% Input Arguments
% enz: NextGen enzyme vector that contains 10 tests.
% time_enz: time vector that contains time (s)
% s_initial: initial substrate concentration for each test (uM)
%
% Output Arguments
% V_init: initial velocity vector (uM/s)
% vmax: maximum velocity that is founded by Lineweaver Burk plot (uM/s)
% km: Michaelis constant which demonstrates the affinity of enzyme for the
% substrate. (uM)
% 
%
% Assignment Information
%   Assignment:     M04, Part 2
%   Team member:    Kyung Min Ko, ko120@purdue.edu
%                   Danny Kim     kim3255@purdue.edu 
%                   Caleb Hancock hancoc25@purdue.edu
%                   John Groves,  jsgroves@purdue.edu               
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
%% Improvement 1
% Changing the window size
 %window = 7; % window size of smooting curve loop(we found this size by trying various size and this gave us the ideal fit of the graph)
 window_smooth = 24; % new window size that fits to the data (24)
 %% Improvement 2
 % Changing the initial percentage of vector
 %lin_percent = 0.0209; % initial percentage of total length of the vector that represent linear region of graph which is about 3 percent (we have found that this percentage gave us best fit tangent line by experiment)
 lin_percent = 0.0129; % new percentage that gives us desired values (0.0129)
%% ____________________
%% CALCULATIONS

% smoothing the curve by taking average values of vector within given window size
half_win = floor(window_smooth/2); % half size of the window 
enz_smooth = enz; % pre-allocating the data to increase the compilation speed 
 
 for indx_c = 1:1:size(enz_smooth,2) % column loop
    for indx_r = half_win + 1 : size(enz_smooth,1) - half_win % row loop
        windowVal = enz_smooth(indx_r - half_win: indx_r + half_win,indx_c); % finding value in the range of window
        newVal = mean(windowVal); % find the average value within given window
        enz_smooth(indx_r - half_win,indx_c) = newVal; % assigning average value to new vector
     end
 end



% finding V0

V_init = zeros(1,size(enz_smooth,2)); % pre-allocating vectors with zeros to increase compilation speed

for indx2 = 1:1:size(enz_smooth,2) % column size of total enzyme vector
    
        data_column= enz_smooth(:,indx2); % extracting one column from data
        vec_no_nan =data_column(~isnan(data_column)); % removing NaN from vector 
        linear_length = round(length(vec_no_nan)*lin_percent); % finding linear region of vector 0.0209 is referencing 2 percent of total length (we have found that this percentage gave us best fit tangent line by experiment)
        linear_sub_y = vec_no_nan(1:linear_length); % initial substrate concentration within linear region (uM)
        linear_time_x = time_enz(1:linear_length); % time within linear portion (s)
        %% Improvement 1 
        % Finding the initial velocity using slope between (1,1) and 1
        % percent points instead of using polyfit to minimize SSE since we
        % have noticed that using polyfit ruins the output for PGOX50
        % enzyme
        %coeff_lin = polyfit(linear_time_x,linear_sub_y,1); %finding slope and y intercept
        %V_init(indx2) = coeff_lin(1); % initial velocity vector(uM/min)
        V_init(indx2) =(linear_sub_y(linear_length) - linear_sub_y(1))/(linear_time_x(linear_length) - linear_time_x(1));
end



% finding Lineweaver Burk plot coefficient by linearizing Michaelis-Menten model
burk_x = 1./ s_initial; % x axis of Lineweaver plot (1/uM)
burk_y = 1./V_init; % y axis of Lineweaver plot (s/uM)

coeff_burk = polyfit(burk_x,burk_y,1); % finding coefficient
vmax = 1/coeff_burk(2); % v max (uM/s)
km = coeff_burk(1) * vmax; % km (uM)


%% ____________________
%% FORMATTED TEXT/FIGURE DISPLAYS

v0_txt = 'The initial velocity (uM/s) vector is\n';
fprintf(v0_txt);
disp(V_init);
km_vmax_txt = 'km: %.3f (uM) vmax: %.3f (uM/s)\n';
fprintf(km_vmax_txt, km,vmax);
    
    
%% ____________________
%% RESULTS
% Reuslt from M2
% Enz A original
% The initial velocity (uM/s) vector is
%   Columns 1 through 8
% 
%     0.0225    0.0440    0.0813    0.1589    0.2652    0.4049    0.5577    0.6835
% 
%   Columns 9 through 10
% 
%     0.7951    0.8756
% 
% km: 141.16 (uM) vmax: 0.87 (uM/s)

% Enz A duplicate 

% The initial velocity (uM/s) vector is
%   Columns 1 through 8
% 
%     0.0219    0.0443    0.0819    0.1576    0.2687    0.4068    0.5596    0.6815
% 
%   Columns 9 through 10
% 
%     0.7972    0.8758
% 
% km: 159.83 (uM) vmax: 0.96 (uM/s)

% Enz B original 
% The initial velocity (uM/s) vector is
%   Columns 1 through 8
% 
%     0.0091    0.0186    0.0357    0.0694    0.1350    0.2219    0.3513    0.4832
% 
%   Columns 9 through 10
% 
%     0.6597    0.7163
% 
% km: 378.45 (uM) vmax: 0.93 (uM/s)

% Enz B duplicate
% The initial velocity (uM/s) vector is
%   Columns 1 through 8
% 
%     0.0093    0.0184    0.0357    0.0674    0.1303    0.2203    0.3508    0.4860
% 
%   Columns 9 through 10
% 
%     0.6611    0.7195
% 
% km: 292.10 (uM) vmax: 0.74 (uM/s)

% Enz C original 
% The initial velocity (uM/s) vector is
%   Columns 1 through 8
% 
%     0.0224    0.0439    0.0839    0.1570    0.2900    0.4554    0.6671    0.8536
% 
%   Columns 9 through 10
% 
%     0.9942    1.0955
% 
% km: 186.19 (uM) vmax: 1.14 (uM/s)

% Enz C duplicate
% The initial velocity (uM/s) vector is
%   Columns 1 through 8
% 
%     0.0246    0.0441    0.0838    0.1565    0.2935    0.4591    0.6724    0.8675
% 
%   Columns 9 through 10
% 
%     0.9824    1.0974
% 
% km: 132.83 (uM) vmax: 0.88 (uM/s)

% Enz D original
% The initial velocity (uM/s) vector is
%   Columns 1 through 8
% 
%     0.0204    0.0385    0.0762    0.1455    0.2745    0.4648    0.7177    0.9752
% 
%   Columns 9 through 10
% 
%     1.1970    1.3744
% 
% km: 236.08 (uM) vmax: 1.29 (uM/s)

% Enz D duplicate
% The initial velocity (uM/s) vector is
%   Columns 1 through 8
% 
%     0.0204    0.0389    0.0764    0.1487    0.2743    0.4633    0.7011    0.9792
% 
%   Columns 9 through 10
% 
%     1.1844    1.3751
% 
% km: 238.13 (uM) vmax: 1.31 (uM/s)

% Enz E original
% The initial velocity (uM/s) vector is
%   Columns 1 through 8
% 
%     0.0372    0.0706    0.1402    0.2284    0.4401    0.6602    0.9671    1.1883
% 
%   Columns 9 through 10
% 
%     1.3442    1.4937
% 
% km: 135.27 (uM) vmax: 1.37 (uM/s)

% Enz E duplicate

% The initial velocity (uM/s) vector is
%   Columns 1 through 8
% 
%     0.0383    0.0726    0.1448    0.2326    0.4559    0.6688    0.9632    1.2038
% 
%   Columns 9 through 10
% 
%     1.3622    1.4848
% 
% km: 131.95 (uM) vmax: 1.38 (uM/s)
%% ____________________
%% ACADEMIC INTEGRITY STATEMENT
% We have not used source code obtained from any other unauthorized
% source, either modified or unmodified. Neither have we provided
% access to my code to another. The program we are submitting
% is our own original work.

