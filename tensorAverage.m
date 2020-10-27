function stressAverage
% average of several stress tensors in reduced stress tensor form
% @ 2020 by Andrea Bistacchi
% code distributed under the GNU AGPL v3.0 license
% created 14/07/2016, last modified 13/05/2020

% Initialize figures etc.
set(0,'DefaultFigureWindowStyle','docked','DefaultFigureColor','w');
set(0,'DefaultAxesFontName','Times','DefaultAxesFontSize',16);

clear all; close all; clc;

figure('Name','Stereoplot'); circle; hold on;  axis equal; axis off; % figure(1)
figure('Name','Mohr plot'); % figure(2)
drawnow; commandwindow % focus on command window

% Load or input data
disp(' ');
disp('Average of several stress tensors in reduced stress tensor form');
disp(' ');
disp('Load input data from [1] previously saved .mat file, [2] table in .txt file, or [3] terminal.')
newData = '0';
while not(newData == '1' | newData == '2' | newData == '3')
    newData = input('1, 2, 3 > ','s');
end

switch newData
    case '1'
        % load data
        [file_name,file_path] = uigetfile('*.mat','Select .mat file with tensor data');
        file_full = [file_path file_name];
        [file_path, file_name, ~] = fileparts(file_full);
        load(file_full,'sigmaTarray','sigma1Plunge','sigma1Trend','sigma2Plunge','sigma2Trend','sigma3Plunge','sigma3Trend');
                
        % number of input data rows
        n = length(sigmaTarray);
                
    case '2'
        % load data
        [file_name,file_path] = uigetfile('*.txt','Select .txt file with tensor data.\n Comumn names must be sigma1Plunge, sigma1Trend sigma3Plunge sigma3Trend R in any order');
        file_full = [file_path file_name];
        [file_path, file_name, ~] = fileparts(file_full);
        delimeters = {',', ' ', ';', '\t'};
        inTable = readtable(file_full,'Delimiter',delimeters);
        
        % number of input data rows
        n = height(inTable);
        
        % initialize stress tensor and axes array
        sigmaTarray = zeros(n,6);
        sigma1Plunge = zeros(n,1);
        sigma1Trend = zeros(n,1);
        sigma2Plunge = zeros(n,1);
        sigma2Trend = zeros(n,1);
        sigma3Plunge = zeros(n,1);
        sigma3Trend = zeros(n,1);
        
        for i = 1:n
            % shape ratio and stress values
            shapeRatio = inTable.R(n);
            sigma1 = 1;
            sigma2 = shapeRatio; % this is so simple since sigma1 = 1 and sigma3 = 0
            sigma3 = 0;
            
            % input and calculate orientation
            sigma1PlungeIn = inTable.sigma1Plunge(i);
            sigma1TrendIn  = inTable.sigma1Trend(i);
            sigma3PlungeIn = inTable.sigma3Plunge(i);
            sigma3TrendIn  = inTable.sigma3Trend(i);
            
            tol = 2; % tol = 2 degrees corrensponds to an angle of at least 88 degrees between input axes
            [sigma1Vers,sigma2Vers,sigma3Vers] = threeOrthoAxes(sigma1PlungeIn,sigma1TrendIn,sigma3PlungeIn,sigma3TrendIn,tol);
            
            if sigma1Vers(1,1) == -999,  % exits while loop if versors are OK according to tolerance tol
                disp(' ');
                disp([' ERROR - principal axes are not orthogonal on tensor: ' num2str(i)]);
                disp(' ');
            end
                        
            [sigma1Plunge(i),sigma1Trend(i)] = plungetrend(sigma1Vers');
            [sigma2Plunge(i),sigma2Trend(i)] = plungetrend(sigma2Vers');
            [sigma3Plunge(i),sigma3Trend(i)] = plungetrend(sigma3Vers');
            
            % build tensor and rotate it
            rotMatrix = [sigma1Vers' ; sigma2Vers' ; sigma3Vers'];
            sigmaT = rotMatrix' * [sigma1 0 0 ; 0 sigma2 0 ; 0 0 sigma3] * rotMatrix;
            sigmaTarray(i,:) = [sigmaT(1,1) sigmaT(2,2) sigmaT(3,3) sigmaT(1,2) sigmaT(1,3) sigmaT(2,3)];
            
        end
        
        % clear table
        clear inTable
        
        % save to disk
        [file_path '\' file_name '.mat']
        [file_name,file_path] = uiputfile('*.mat','Save file',[file_path '\' file_name '.mat']);
        file_full = [file_path file_name];
        save(file_full,'sigmaTarray','sigma1Plunge','sigma1Trend','sigma2Plunge','sigma2Trend','sigma3Plunge','sigma3Trend');             

    case '3'
        disp(' ');
        disp('Number of input tensors');
        n = input(' > ');
        
        % initialize stress tensor and axes array
        sigmaTarray = zeros(n,6);
        sigma1Plunge = zeros(n,1);
        sigma1Trend = zeros(n,1);
        sigma2Plunge = zeros(n,1);
        sigma2Trend = zeros(n,1);
        sigma3Plunge = zeros(n,1);
        sigma3Trend = zeros(n,1);
        
        for i = 1:n
            disp(' ');
            disp(['Input tensor n. ' num2str(i)]);
            
            % shape ratio and stress values
            disp(' ');
            shapeRatio = input('Shape Ratio R (sigma2 - sigma3)/(sigma1 - sigma3) [] [0.5] > '); if isempty(shapeRatio), shapeRatio = 0.5; end;
            sigma1 = 1;
            sigma2 = shapeRatio; % this is so simple since sigma1 = 1 and sigma3 = 0
            sigma3 = 0;
            
            % input and calculate orientation
            sigma1Vers = [1;0;0]; % initialize
            while 1
                disp(' ');
                sigma1PlungeIn = input(' sigma1 Plunge > ');
                sigma1TrendIn  = input(' sigma1 Trend  > ');
                sigma3PlungeIn = input(' sigma3 Plunge > ');
                sigma3TrendIn  = input(' sigma3 Trend  > ');
                
                tol = 2; % tol = 2 degrees corrensponds to an angle of at least 88 degrees between input axes
                [sigma1Vers,sigma2Vers,sigma3Vers] = threeOrthoAxes(sigma1PlungeIn,sigma1TrendIn,sigma3PlungeIn,sigma3TrendIn,tol);
                
                if sigma1Vers(1,1) ~= -999,  % exits while loop if versors are OK according to tolerance tol
                    break
                else
                    disp(' ');
                    disp(' ERROR - principal axes are not orthogonal');
                    disp(' ');
                end
            end
            
            [sigma1Plunge(i),sigma1Trend(i)] = plungetrend(sigma1Vers');
            [sigma2Plunge(i),sigma2Trend(i)] = plungetrend(sigma2Vers');
            [sigma3Plunge(i),sigma3Trend(i)] = plungetrend(sigma3Vers');
            
            % build tensor and rotate it
            rotMatrix = [sigma1Vers' ; sigma2Vers' ; sigma3Vers'];
            sigmaT = rotMatrix' * [sigma1 0 0 ; 0 sigma2 0 ; 0 0 sigma3] * rotMatrix;
            sigmaTarray(i,:) = [sigmaT(1,1) sigmaT(2,2) sigmaT(3,3) sigmaT(1,2) sigmaT(1,3) sigmaT(2,3)];
            
            % stereoplot
            figure(3); axis equal; axis off;
            [sigma1x,sigma1y] = plotLine(sigma1Plunge(i),sigma1Trend(i));
            [sigma2x,sigma2y] = plotLine(sigma2Plunge(i),sigma2Trend(i));
            [sigma3x,sigma3y] = plotLine(sigma3Plunge(i),sigma3Trend(i));
            plot(sigma1x,sigma1y,'o','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',8);
            plot(sigma2x,sigma2y,'o','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',8);
            plot(sigma3x,sigma3y,'o','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',8);
            
            drawnow; commandwindow % focus on command window
            
        end
        
        % save to disk
        [file_name,file_path] = uiputfile('*.mat');
        file_full = [file_path file_name];
        save(file_full,'sigmaTarray','sigma1Plunge','sigma1Trend','sigma2Plunge','sigma2Trend','sigma3Plunge','sigma3Trend');             

end

% info
disp(' ');
disp(['File name: ' file_name]);
disp(['Number of input tensors: ' num2str(n)]);

% stereoplot
figure(1);
[sigma1x,sigma1y] = plotLine(sigma1Plunge,sigma1Trend);
[sigma2x,sigma2y] = plotLine(sigma2Plunge,sigma2Trend);
[sigma3x,sigma3y] = plotLine(sigma3Plunge,sigma3Trend);
plot(sigma1x,sigma1y,'o','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',8);
plot(sigma2x,sigma2y,'o','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',8);
plot(sigma3x,sigma3y,'o','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',8);

% calculate mean tensor
outVector = meanReducedTensor(sigmaTarray);

% reformat
sigma1avr = outVector(1);
sigma2avr = outVector(2);
sigma3avr = outVector(3);
sigma1avrPlunge = outVector(4);
sigma1avrTrend = outVector(5);
sigma2avrPlunge = outVector(6);
sigma2avrTrend = outVector(7);
sigma3avrPlunge = outVector(8);
sigma3avrTrend = outVector(9);
shapeRatioAvr = outVector(10);

% _________________________________________________________________________
% calculate dispersion of input data

% [sigma1dispersAxMin,sigma1dispersAxMax,sigma1dispersEllipse] = confEllipse(sigma1Plunge,sigma1Trend,sigma1avrPlunge,sigma1avrTrend,0);
% [sigma2dispersAxMin,sigma2dispersAxMax,sigma2dispersEllipse] = confEllipse(sigma2Plunge,sigma2Trend,sigma2avrPlunge,sigma2avrTrend,0);
% [sigma3dispersAxMin,sigma3dispersAxMax,sigma3dispersEllipse] = confEllipse(sigma3Plunge,sigma3Trend,sigma3avrPlunge,sigma3avrTrend,0);

% _________________________________________________________________________


% simulate nboot mean tensors with bootstrapping
% bootstat gives values simulated from function @meanReducedTensor
% ci gives 95% confidence intervals - these can be difficult to interpret
% due to nonlinear relationship of Plunge/Trend in polar coords
nboot = 100;
[ci,bootstat] = bootci(nboot,@meanReducedTensor,sigmaTarray);

% reformat bootstrapping output

% bootstat is the bootstrapped statistic computed for each of the nboot bootstrap
% replicate samples. Each one of the nboot rows of bootstat contains the results
% of applying the meanReducedTensor function to one bootstrap sample. Columns are
% the standard output from meanReducedTensor.
sigma1avrBS = bootstat(:,1);
sigma2avrBS = bootstat(:,2);
sigma3avrBS = bootstat(:,3);
sigma1avrPlungeBS = bootstat(:,4);
sigma1avrTrendBS = bootstat(:,5);
sigma2avrPlungeBS = bootstat(:,6);
sigma2avrTrendBS = bootstat(:,7);
sigma3avrPlungeBS = bootstat(:,8);
sigma3avrTrendBS = bootstat(:,9);
shapeRatioAvrBS = bootstat(:,10);

% _________________________________________________________________________
% calculate dispersion of bootstrapped averages

[sigma1confidAxMin,sigma1confidAxMax,sigma1confidEllipsePlunge,sigma1confidEllipseTrend] = confEllipse(sigma1avrPlungeBS,sigma1avrTrendBS,sigma1avrPlunge,sigma1avrTrend,0);
[sigma2confidAxMin,sigma2confidAxMax,sigma2confidEllipsePlunge,sigma2confidEllipseTrend] = confEllipse(sigma2avrPlungeBS,sigma2avrTrendBS,sigma2avrPlunge,sigma2avrTrend,0);
[sigma3confidAxMin,sigma3confidAxMax,sigma3confidEllipsePlunge,sigma3confidEllipseTrend] = confEllipse(sigma3avrPlungeBS,sigma3avrTrendBS,sigma3avrPlunge,sigma3avrTrend,0);


% _________________________________________________________________________

%ELIMINA ESCLUSO SHAPE RATIO_______________________________________________
% ci includes upper and lower confidence intervals (scalars) obtained from bootci, but they
% are difficult to interpret due to nonlinear relationship of Plunge/Trend in polar coords
sigma1avrBSlowCI = ci(1,1);
sigma2avrBSlowCI = ci(1,2);
sigma3avrBSlowCI = ci(1,3);
sigma1avrPlungeBSlowCI = ci(1,4);
sigma1avrTrendBSlowCI = ci(1,5);
sigma2avrPlungeBSlowCI = ci(1,6);
sigma2avrTrendBSlowCI = ci(1,7);
sigma3avrPlungeBSlowCI = ci(1,8);
sigma3avrTrendBSlowCI = ci(1,9);
shapeRatioAvrBSlowCI = ci(1,10);

sigma1avrBShighCI = ci(2,1);
sigma2avrBShighCI = ci(2,2);
sigma3avrBShighCI = ci(2,3);
sigma1avrPlungeBShighCI = ci(2,4);
sigma1avrTrendBShighCI = ci(2,5);
sigma2avrPlungeBShighCI = ci(2,6);
sigma2avrTrendBShighCI = ci(2,7);
sigma3avrPlungeBShighCI = ci(2,8);
sigma3avrTrendBShighCI = ci(2,9);
shapeRatioAvrBShighCI = ci(2,10);

%ELIMINA___________________________________________________________________
% get differences between sample mean and bootstrapped mean
deltaSigma1Plunge = sigma1avrPlungeBS - sigma1avrPlunge;
deltaSigma1Trend  = sigma1avrTrendBS  - sigma1avrTrend;
deltaSigma2Plunge = sigma2avrPlungeBS - sigma2avrPlunge;
deltaSigma2Trend  = sigma2avrTrendBS  - sigma2avrTrend;
deltaSigma3Plunge = sigma3avrPlungeBS - sigma3avrPlunge;
deltaSigma3Trend  = sigma3avrTrendBS  - sigma3avrTrend;

deltaSigma1PlungeLow = quantile(deltaSigma1Plunge,0.05);
deltaSigma1TrendLow  = quantile(deltaSigma1Trend,0.05);
deltaSigma2PlungeLow = quantile(deltaSigma2Plunge,0.05);
deltaSigma2TrendLow  = quantile(deltaSigma2Trend,0.05);
deltaSigma3PlungeLow = quantile(deltaSigma3Plunge,0.05);
deltaSigma3TrendLow  = quantile(deltaSigma3Trend,0.05);

deltaSigma1PlungeHigh = quantile(deltaSigma1Plunge,0.95);
deltaSigma1TrendHigh  = quantile(deltaSigma1Trend,0.95);
deltaSigma2PlungeHigh = quantile(deltaSigma2Plunge,0.95);
deltaSigma2TrendHigh  = quantile(deltaSigma2Trend,0.95);
deltaSigma3PlungeHigh = quantile(deltaSigma3Plunge,0.95);
deltaSigma3TrendHigh  = quantile(deltaSigma3Trend,0.95);

deltaRatio = shapeRatioAvrBS - shapeRatioAvr;

deltaRatioLow = quantile(deltaRatio,0.05);
deltaRatioHigh = quantile(deltaRatio,0.95);


% Display values
disp(' ');
disp('Principal stresses (Plunge / Trend / Value):');
disp(['Sigma1: ' num2str(sigma1avrPlunge) ' / ' num2str(sigma1avrTrend) ' / '  num2str(sigma1avr,'%+10.5g')]);
disp(['Sigma2: ' num2str(sigma2avrPlunge) ' / ' num2str(sigma2avrTrend) ' / '  num2str(sigma2avr,'%+10.5g')]);
disp(['Sigma3: ' num2str(sigma3avrPlunge) ' / ' num2str(sigma3avrTrend) ' / '  num2str(sigma3avr,'%+10.5g')]);
disp(' ');
disp(['Shape Ratio (sigma2 - sigma3)/(sigma1 - sigma3): ' num2str(shapeRatioAvr)]);
disp(' ');
disp('95% confidence intervals from bootstrapping');
disp(' ');
disp('Principal stresses confidence ellipse (Short Axis / Long Axis):');
disp(['Sigma1: ' num2str(sigma1confidAxMin) ' ÷ ' num2str(sigma1confidAxMax)]);
disp(['Sigma2: ' num2str(sigma2confidAxMin) ' ÷ ' num2str(sigma2confidAxMax)]);
disp(['Sigma3: ' num2str(sigma3confidAxMin) ' ÷ ' num2str(sigma3confidAxMax)]);
disp(' ');
disp(['Shape Ratio (sigma2 - sigma3)/(sigma1 - sigma3): ' num2str(shapeRatioAvrBSlowCI) ' ÷ ' num2str(shapeRatioAvrBShighCI)]);
disp(' ');

% stereoplot
figure(1);
sigma1avrXBS = zeros(size(sigma1avrPlungeBS));
sigma1avrYBS = zeros(size(sigma1avrPlungeBS));
sigma2avrXBS = zeros(size(sigma1avrPlungeBS));
sigma2avrYBS = zeros(size(sigma1avrPlungeBS));
sigma3avrXBS = zeros(size(sigma1avrPlungeBS));
sigma3avrYBS = zeros(size(sigma1avrPlungeBS));

for i=1:length(sigma1avrPlungeBS)
    [sigma1avrXBS(i),sigma1avrYBS(i)] = plotLine(sigma1avrPlungeBS(i),sigma1avrTrendBS(i));
    [sigma2avrXBS(i),sigma2avrYBS(i)] = plotLine(sigma2avrPlungeBS(i),sigma2avrTrendBS(i));
    [sigma3avrXBS(i),sigma3avrYBS(i)] = plotLine(sigma3avrPlungeBS(i),sigma3avrTrendBS(i));
end

plot(sigma1avrXBS,sigma1avrYBS,'h','LineWidth',1,'MarkerEdgeColor','b','MarkerFaceColor','w','MarkerSize',4);
plot(sigma2avrXBS,sigma2avrYBS,'p','LineWidth',1,'MarkerEdgeColor','g','MarkerFaceColor','w','MarkerSize',4);
plot(sigma3avrXBS,sigma3avrYBS,'d','LineWidth',1,'MarkerEdgeColor','r','MarkerFaceColor','w','MarkerSize',4);

[sigma1avrX,sigma1avrY] = plotLine(sigma1avrPlunge,sigma1avrTrend);
[sigma2avrX,sigma2avrY] = plotLine(sigma2avrPlunge,sigma2avrTrend);
[sigma3avrX,sigma3avrY] = plotLine(sigma3avrPlunge,sigma3avrTrend);

plot(sigma1avrX,sigma1avrY,'h','LineWidth',2,'MarkerEdgeColor','b','MarkerFaceColor','w','MarkerSize',15);
plot(sigma2avrX,sigma2avrY,'p','LineWidth',2,'MarkerEdgeColor','g','MarkerFaceColor','w','MarkerSize',15);
plot(sigma3avrX,sigma3avrY,'d','LineWidth',2,'MarkerEdgeColor','r','MarkerFaceColor','w','MarkerSize',12);

[sigma1confidEllipseX,sigma1confidEllipseY] = plotLine(sigma1confidEllipsePlunge,sigma1confidEllipseTrend);
[sigma2confidEllipseX,sigma2confidEllipseY] = plotLine(sigma2confidEllipsePlunge,sigma2confidEllipseTrend);
[sigma3confidEllipseX,sigma3confidEllipseY] = plotLine(sigma3confidEllipsePlunge,sigma3confidEllipseTrend);

plot(sigma1confidEllipseX,sigma1confidEllipseY,'-b','LineWidth',2);
plot(sigma2confidEllipseX,sigma2confidEllipseY,'-g','LineWidth',2);
plot(sigma3confidEllipseX,sigma3confidEllipseY,'-r','LineWidth',2);

drawnow; commandwindow % focus on command window

% MOHR PLOT

% Mohr circles
[unitCircleX,unitCircleY] = pol2cart(linspace(0,2*pi,90),ones(1,90));
circle13X = unitCircleX*(sigma1avr-sigma3avr)/2+(sigma1avr+sigma3avr)/2;
circle12X = unitCircleX*(sigma1avr-sigma2avr)/2+(sigma1avr+sigma2avr)/2;
circle23X = unitCircleX*(sigma2avr-sigma3avr)/2+(sigma2avr+sigma3avr)/2;
circle13Y = unitCircleY*(sigma1avr-sigma3avr)/2;
circle12Y = unitCircleY*(sigma1avr-sigma2avr)/2;
circle23Y = unitCircleY*(sigma2avr-sigma3avr)/2;

% plot
figure(2); hold on;
title('Mohr plot'); axis equal; axis off;
plot(circle13X,circle13Y,'-b','LineWidth',2);
plot(circle12X,circle12Y,'-b','LineWidth',2);
plot(circle23X,circle23Y,'-b','LineWidth',2);
plot([-0.2 1.2],[0 0],'-k','LineWidth',2);  % horizontal axis
plot([0 0],[-0.6 0.6],'-k','LineWidth',2);  % vertical axis
drawnow; commandwindow % focus on command window


end


%%
function outVector = meanReducedTensor(sigmaTarray)

% calculate average
sigmaTavrArray = mean(sigmaTarray,1);

% reformat
sigmaTavr = [sigmaTavrArray(1) sigmaTavrArray(4) sigmaTavrArray(5);...
             sigmaTavrArray(4) sigmaTavrArray(2) sigmaTavrArray(6);...
             sigmaTavrArray(5) sigmaTavrArray(6) sigmaTavrArray(3)];

% eigenvectors / eigenvalues
[sigma1avr,sigma2avr,sigma3avr,sigma1avrPlunge,sigma1avrTrend,sigma2avrPlunge,sigma2avrTrend,sigma3avrPlunge,sigma3avrTrend] = formattedEig(sigmaTavr);

if sigma1avrPlunge < 0,
    sigma1avrPlunge = -sigma1avrPlunge;
    sigma1avrTrend = sigma1avrTrend + 180;
    if sigma1avrTrend > 360,
        sigma1avrTrend = sigma1avrTrend - 360;
    end
end

if sigma2avrPlunge < 0,
    sigma2avrPlunge = -sigma2avrPlunge;
    sigma2avrTrend = sigma2avrTrend + 180;
    if sigma2avrTrend > 360,
        sigma2avrTrend = sigma2avrTrend - 360;
    end
end

if sigma3avrPlunge < 0,
    sigma3avrPlunge = -sigma3avrPlunge;
    sigma3avrTrend = sigma3avrTrend + 180;
    if sigma3avrTrend > 360,
        sigma3avrTrend = sigma3avrTrend - 360;
    end
end

% shape ratio
shapeRatioAvr = (sigma2avr - sigma3avr)/(sigma1avr - sigma3avr);

% arrange output in one single row vector as required by bootci
outVector = [sigma1avr,sigma2avr,sigma3avr,sigma1avrPlunge,sigma1avrTrend,sigma2avrPlunge,sigma2avrTrend,sigma3avrPlunge,sigma3avrTrend,shapeRatioAvr];

end


%%
function circle;  % empty circular plots

[XC,YC] = pol2cart(linspace(0,2*pi,360/4),ones(1,360/4));

plot(XC,YC,'-k','LineWidth',2);

end


%%
function [vers1,vers2,vers3] = threeOrthoAxes(vers1Plunge,vers1Trend,vers3Plunge,vers3Trend,tol)

% returns a DEXTRAL set of orthogonal axes
% vers1 and vers3 are the first and last axes resp. expressed as plunge/trend
% tol could be 2 degrees corrensponding to an angle of at least 88 degrees
% between input axes
%
% Andrea Bistacchi 27/10/2016

% column versors from plunge/trend
vers1 = lineation(vers1Plunge,vers1Trend)';
vers3 = lineation(vers3Plunge,vers3Trend)';

% vers2 from normalized cross product
vers2 = -cross(vers1,vers3);
modVers2 = sqrt(vers2(1,1).^2 + vers2(2,1).^2 + vers2(3,1).^2);
vers2 = vers2/modVers2;

% recalculate vers 3 to ensure orthogonality
vers3 = cross(vers1,vers2);

% if modulus of vers2 is shorter than tol, gives an error code
if modVers2 < sin((90-tol)*pi/180),
    vers1 = [-999; -999;-999]; % error code
    vers2 = [-999; -999;-999]; % error code
    vers3 = [-999; -999;-999]; % error code
end

end

%%
function directionCosines = lineation(plunges,trends)

% directionCosines are row unit vector format, one row for each vector
%
% cos1 = East
% cos2 = North
% cos3 = Upward
%
% Andrea Bistacchi 27/10/2016, modified 6/2/2019 to process several inputs
% at once and obtain row vectors

% check plunge trend format
if size(plunges,2)>1 , plunges = plunges'; end
if size(trends,2)>1 , trends = trends'; end

directionCosines = [cos(plunges*pi/180).*sin(trends*pi/180)   cos(plunges*pi/180).*cos(trends*pi/180)   -sin(plunges*pi/180)];

end

%%
function [plunge,trend] = plungetrend(directionCosines)

% directionCosines must be in row unit vector format, one row for each vector, with
%
% cos1 = East
% cos2 = North
% cos3 = Upward
%
% Andrea Bistacchi 27/10/2016, modified 6/2/2019 to process several inputs

trend = atan2(directionCosines(:,1),directionCosines(:,2))*180/pi;  % IMPORTANT. DO NOT USE ATAN! IT IS LIMITED TO -PI/2 PI/2.
plunge = -asin(directionCosines(:,3))*180/pi;

% check NaNs and negative trends
trend(isnan(trend)) = 0;
trend(trend<0) = trend(trend<0) + 360;

end


%%
function [Xl,Yl] = plotLine(Plunge,Trend);

% check for upwards vectors (with negative Plunge)
Trend = Trend + 180*(sign(Plunge)==-1);
Plunge = Plunge .* sign(Plunge);

rho = sqrt(2).*sin(pi/4-Plunge*pi/180./2);   %projected distance from origin in equiareal Lambert poj. (Schmidt net)
Xl = rho .* sin(Trend*pi/180);
Yl = rho .* cos(Trend*pi/180);

% find reflections in stereoplot and insert NaN to avoid plotting them
reflectId = [find(abs(Trend(1:end-1) - Trend(2:end))>90) find(not(sign(Trend(1:end-1)) == sign(Trend(2:end))))];

if isreal(reflectId)
    for i = 1:length(reflectId)
        Xl = [Xl(1:reflectId(i)) ; NaN ; Xl(reflectId(i)+1:end)];
        Yl = [Yl(1:reflectId(i)) ; NaN ; Yl(reflectId(i)+1:end)];
        reflectId = reflectId + 1;
    end
end

end

%%
function [eigVal1,eigVal2,eigVal3,eigVect1Plunge,eigVect1Trend,eigVect2Plunge,eigVect2Trend,eigVect3Plunge,eigVect3Trend] = formattedEig(Tensor);

[eigVect,eigVal] = eig(Tensor);

[eigVal1,column1] = max([eigVal(1,1) eigVal(2,2) eigVal(3,3)]);
[eigVal3,column3] = min([eigVal(1,1) eigVal(2,2) eigVal(3,3)]);
column2 = 6 - column1 - column3;
eigVal2 = eigVal(column2,column2);

eigVect1 = eigVect(:,column1);
eigVect2 = eigVect(:,column2);
eigVect3 = eigVect(:,column3);

[eigVect1Plunge,eigVect1Trend] = plungetrend(eigVect1');
[eigVect2Plunge,eigVect2Trend] = plungetrend(eigVect2');
[eigVect3Plunge,eigVect3Trend] = plungetrend(eigVect3');

end

%%
function R = rotMatrix(A,B)
% function to calculate a rotation matrix that can be used to rotate in 3D
% a row vector A in order to align it to row vector B
%
% to rotate any vector according to this transformation use vRotated = R*v
%
% to rotate a tensor use tRotated = R'*t*R
%

% check row vectors
if size(A,2)>size(A,1), A = A'; end
if size(B,2)>size(B,1), B = B'; end

% rotation vector
rv = cross(A,B);

% skew-symmetric cross-product matrix
sscm = [ 0    -rv(3)  rv(2);...
       rv(3)     0   -rv(1);...
      -rv(2)  rv(1)     0  ];

% rotation matrix
R = eye(3) + sscm + sscm^2*(1-dot(A,B))/(norm(rv))^2;
     
end

%%
function [cntLineDistMin,cntLineDistMax,confidEllipsePlunge,confidEllipseTrend] = confEllipse(plunges,trends,avrPlunge,avrTrend,plotToggle)

if max(size(plunges))<20, return; end

% check that plunges & trends are column vectors
if size(plunges,2) > size(plunges,1), plunges = plunges'; end
if size(trends,2) > size(trends,1), trends = trends'; end

% convert plunges trends to unit vectors in radians
dataVectors = lineation(plunges,trends);
avrVector = lineation(avrPlunge,avrTrend);

% check that single vectors are oriented as mean vector
sameSense = (dataVectors * avrVector')>=0;
dataVectors = dataVectors.*sameSense - dataVectors.*not(sameSense);

% rotation matrix A used to bring data in a reference centered on
% the average vector avrVector with U V W = 0 0 1
Rot = rotMatrix(avrVector,[0 0 1]);

% rotate data to average vector reference
avrVectorRot = (Rot*avrVector')';
dataVectorsRot = (Rot*dataVectors')';

% compute arc distances in radians along new rotated axes U and V
arcDistances = asin([dataVectorsRot(:,1) dataVectorsRot(:,2)]);

%______________________________________________Mahalanobis contour

% calculate Mahalanobis distance
mahalDistance = mahal(arcDistances,arcDistances);

% calculate 95th percentile of d2_mahal
pc95val = prctile(mahalDistance, 95);

% interpolate on regular grid
minU = min(arcDistances(:,1));
maxU = max(arcDistances(:,1));
minV = min(arcDistances(:,2));
maxV = max(arcDistances(:,2));

interpFunc = scatteredInterpolant(arcDistances, mahalDistance);
interpFunc.Method = 'natural';
[U,V] = meshgrid( linspace(minU,maxU,1000), linspace(minV,maxV,1000) );
interpMahalDistance = interpFunc(U,V);

% calculate 95% contour
cntLine = contourc(linspace(minU,maxU,1000), linspace(minV,maxV,1000), interpMahalDistance, [pc95val pc95val]);
cleanId = find(not(cntLine(1,:)==pc95val));
cntLine = cntLine(:,cleanId)';

% find max and min distance contour to mean
cntLineDist = sqrt(cntLine(:,1).^2 + cntLine(:,2).^2);
cntLineDistMin = min(abs(cntLineDist));
cntLineDistMax = max(abs(cntLineDist));


%______________________________________________Covariance ellipse

% Calculate the eigenvectors and eigenvalues
covariance = cov(arcDistances);
[eigenvec, eigenval ] = eig(covariance);

% Get the index of the largest eigenvector
[largest_eigenvec_ind_c, r] = find(eigenval == max(max(eigenval)));
largest_eigenvec = eigenvec(:, largest_eigenvec_ind_c);

% Get the largest eigenvalue
largest_eigenval = max(max(eigenval));

% Get the smallest eigenvector and eigenvalue
if(largest_eigenvec_ind_c == 1)
    smallest_eigenval = max(eigenval(:,2));
    smallest_eigenvec = eigenvec(:,2);
else
    smallest_eigenval = max(eigenval(:,1));
    smallest_eigenvec = eigenvec(1,:);
end

% Calculate the angle between the x-axis and the largest eigenvector
angle = atan2(largest_eigenvec(2), largest_eigenvec(1));

% This angle is between -pi and pi.
% Let's shift it such that the angle is between 0 and 2pi
if(angle < 0)
    angle = angle + 2*pi;
end

% Get the coordinates of the data mean
avg = mean(arcDistances);

% Get the 95% confidence interval error ellipse
chisquare_val = 2.4477;
theta_grid = linspace(0,2*pi,72);
phi = angle;
X0=avg(1);
Y0=avg(2);
a=chisquare_val*sqrt(largest_eigenval);
b=chisquare_val*sqrt(smallest_eigenval);

% the ellipse in x and y coordinates 
ellipse_x_r  = a*cos( theta_grid );
ellipse_y_r  = b*sin( theta_grid );

%Define a rotation matrix
R = [ cos(phi) sin(phi);...
     -sin(phi) cos(phi) ];

%let's rotate the ellipse to some angle phi
r_ellipse = [ellipse_x_r;ellipse_y_r]' * R;

%______________________________________________Calculate corrected ellipse
%
% In nice situations the Mahalanobis contour and Covariance ellipse are
% almost the same. For very flattened distributions however the ellipse
% might overestimate errors. The corrected ellipse calculated here has axes
% oriented as the covariance ellipse and axes obtained from min and max
% axes of the Mahalanobis contour.

% the corrected ellipse in x and y coordinates 
correctedEllipseX  = cntLineDistMax*cos( theta_grid );
correctedEllipseY  = cntLineDistMin*sin( theta_grid );

%Define a rotation matrix
R = [ cos(phi) sin(phi);...
     -sin(phi) cos(phi) ];

%let's rotate the ellipse to some angle phi
correctedEllipse = [correctedEllipseX;correctedEllipseY]' * R;

%______________________________________________Plot

if plotToggle == 1
    
    % Plot the original data
    figure; hold on
    scatter(arcDistances(:,1), arcDistances(:,2), 20, mahalDistance,'o','filled')
    
    % plot gridded data
    %scatter(U(:), V(:), 20, interpMahalDistance(:),'o')
    
    %plot contour
    plot(cntLine(:,1),cntLine(:,2),'LineWidth',2);
    
    % plot max and min circle
    a = linspace(0,2*pi);
    plot(cos(a)*cntLineDistMin,sin(a)*cntLineDistMin)
    plot(cos(a)*cntLineDistMax,sin(a)*cntLineDistMax)
    
    % Draw the covariance ellipse
    plot(r_ellipse(:,1) + X0,r_ellipse(:,2) + Y0,'-')
    
    % Draw the corrected ellipse
    plot(correctedEllipse(:,1),correctedEllipse(:,2),'-','LineWidth',2)
    
    % cosmetics
    hb = colorbar;
    ylabel(hb,'Mahalanobis Distance')
    axis equal
    xlabel('U');
    ylabel('V');
    
end

%______________________________________________Return corrected ellipse

% add W coordinate = 1 (unit sphere) and rotate corrected ellipse back to
% input reference frame
correctedEllipseX = sin(correctedEllipse(:,1));
correctedEllipseY = sin(correctedEllipse(:,2));
correctedEllipseZ = cos(asin(sqrt(correctedEllipseX.^2 + correctedEllipseY.^2)));

correctedEllipseGeol = (Rot'*[correctedEllipseX correctedEllipseY correctedEllipseZ]')';

if plotToggle == 1
    figure; hold on
    plot3(correctedEllipseGeol(:,1),correctedEllipseGeol(:,2),correctedEllipseGeol(:,3),'k-','LineWidth',2)
    plot3(dataVectors(:,1),dataVectors(:,2),dataVectors(:,3),'o')
end

% convert to plunge trend
[confidEllipsePlunge,confidEllipseTrend] = plungetrend(correctedEllipseGeol);

% close line
confidEllipsePlunge = [confidEllipsePlunge ; confidEllipsePlunge(1,:)];
confidEllipseTrend = [confidEllipseTrend ; confidEllipseTrend(1,:)];

% convert to degrees
cntLineDistMin = cntLineDistMin*180/pi;
cntLineDistMax = cntLineDistMax*180/pi;

end