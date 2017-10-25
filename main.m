clear; close all;

load('Resistance+ThermalVoltage.mat')
load('bode.mat')

%% Finding Gain Bandwidth Product
Vin                 = 0.2; %V Set in experiment
Attenutation_Gain   = -85;  %db Set in experiment
RMS_Scale           = 1 / (2 * sqrt(2));
Voltage_Gain_Scale  = 10^( Attenutation_Gain / 20 );
VinScaled           = Vin * RMS_Scale * Voltage_Gain_Scale; %V

T                   = 296; %K Measured
deltaT              = 2;   %K
T_MAX               = 296 + deltaT;
T_MIN               = 296 - deltaT;

z = cumtrapz(FrequencyHz, (MeanVV/VinScaled) .^2 );
gamma = z(length(z));
disp(['The gain-bandwidth product is: ', num2str(gamma)]);

zmax        = cumtrapz(FrequencyHz, ((MeanVV+STDV)/VinScaled) .^2 );
zmin        = cumtrapz(FrequencyHz, ((MeanVV-STDV)/VinScaled) .^2 );
gammaMax    = zmax(length(zmax));
gammaMin    = zmin(length(zmin));
deltaGamma  =  (gammaMax  - gammaMin) / 2;

%V2overR = 0.0003069;

%gamma = 5008693.576 ;

%% Comparison of Gamma to the more rigourous method
load('vivout.mat')
rigor_Vin = [0.05 0.1 0.15 0.2 0.25];
rigor_Vin_Scaled = rigor_Vin * RMS_Scale * Voltage_Gain_Scale;
rigor_f = [100 1e3 5e3 10e3 100e3];
G100    = Vout100Hz     .^ 2 ./ rigor_Vin_Scaled .^2;
G1k     = Vout1kHz      .^ 2 ./ rigor_Vin_Scaled .^2;
G5k     = Vout5kHz      .^ 2 ./ rigor_Vin_Scaled .^2;
G10k    = Vout10kHz     .^ 2 ./ rigor_Vin_Scaled .^2;
G100k   = Vout100kHz    .^ 2 ./ rigor_Vin_Scaled .^2;
G = [mean(G100); mean(G1k); mean(G5k); mean(G10k); mean(G100k)];
rigor_gamma = G(1)*(rigor_f(1)/2 - 0) + G(2) * (rigor_f(2)/2 - rigor_f(1)/2) + G(3) * (rigor_f(3)/2 - rigor_f(2)/2) + G(4) * (rigor_f(4)/2 - rigor_f(3)/2) + ...
    G(5) * (rigor_f(5)/2 - rigor_f(4)/2);
disp(['The alternate gain-bandwidth product is: ', num2str(rigor_gamma)]);


%% Resistors A-J

%% Finding Voltage to Resistance Relation
resistanceAJ        = Resistances(1:10);
squareVoltageAJ     = (ThermalVoltage(1:10)) .^ 2;
errorsResistanceAJ  = STDR(1:10);
errorsVoltageAJ     = STDTV(1:10);
weightsAJ           = 1 ./ (errorsVoltageAJ).^2;
PAJ                 = polyfit(resistanceAJ,squareVoltageAJ, 1 );
V2R_AJ              = PAJ(1);
V2R_AJ_MAX          = 0.0003085; %from Weighted Least Square Fit
V2R_AJ_MIN          = 0.0003054; %from Weighted Least Square Fit

kb_1                =  V2R_AJ / (4*T*gamma);
kb_1Max             =  V2R_AJ_MAX / (4 * T_MIN * gammaMin);
kb_1Min             =  V2R_AJ_MIN / (4 * T_MAX * gammaMax);
delta_kb_1          = ( kb_1Max - kb_1Min ) / 2;
disp(['The Boltzmann constant from the first set of resistors is: ', num2str(kb_1), ' plus/minus: ', num2str(delta_kb_1)]);

%% Resistors 1-9

%% Finding Voltage to Resistance Relation
resistanceNum        = Resistances(11:19);
squareVoltageNum     = (ThermalVoltage(11:19)) .^ 2;
errorsResistanceNum  = STDR(11:19);
errorsVoltageNum     = STDTV(11:19);
weightsNum           = 1 ./ (errorsVoltageNum).^2;
PNUM                 = polyfit(resistanceNum,squareVoltageNum, 1 );
V2R_NUM              = PNUM(1);
V2R_NUM_MAX          = 0.0003081; %from Weighted Least Square Fit
V2R_NUM_MIN          = 0.0003035; %from Weighted Least Square Fit

kb_2                =  V2R_NUM / (4*T*gamma);
kb_2Max             =  V2R_NUM_MAX / (4 * T_MIN * gammaMin);
kb_2Min             =  V2R_NUM_MIN / (4 * T_MAX * gammaMax);
delta_kb_2          = ( kb_2Max - kb_2Min ) / 2;
disp(['The Boltzmann constant from the second set of resistors is: ', num2str(kb_2), ' plus/minus: ', num2str(delta_kb_2)]);


%% Plots 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Voltage vs Resistance A-J %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f1 = figure(1);
h1 = plot(resistanceAJ, squareVoltageAJ, 'r.');
hold on
e  = errorbar(resistanceAJ, squareVoltageAJ, errorsVoltageAJ/2, errorsVoltageAJ/2, errorsResistanceAJ/2, errorsResistanceAJ/2, 'r.');
e.LineWidth = 1.5;
h2 = plot(resistanceAJ, PAJ(1)*resistanceAJ + PAJ(2), 'Linewidth', 1.5);

ax = gca;
ax.FontName = 'LaTeX';
ax.TickLabelInterpreter = 'LaTeX';
ax.FontSize = 18;
ax.XColor = 'k';
ax.YColor = 'k';

ax.YLabel.String = ('Voltage (V)');
ax.YLabel.FontSize = 16;
ax.YLabel.FontWeight = 'bold';
ax.YLabel.Color = 'k';

ax.XLabel.String = ('Resistance (\Omega)');
ax.XLabel.FontSize = 16;
ax.XLabel.FontWeight = 'bold';
ax.XLabel.Color = 'k';
ax.XLim = [0, 10050];

ax.Box = 'off';
ax.LineWidth = 1.5;
ax.YGrid = 'on';
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';

t = title('Thermal Voltage Vs. Resistance A-J');
t.Color = 'k';
t.Interpreter = 'LaTeX';
t.FontSize = 24;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Voltage vs Resistance 1-9 %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f2 = figure(2);
h3 = plot(resistanceNum, squareVoltageNum, 'or');
hold on
e  = errorbar(resistanceNum, squareVoltageNum, errorsVoltageNum/2, errorsVoltageNum/2, errorsResistanceNum/2, errorsResistanceNum/2, 'r.');
e.LineWidth = 1.5;
h4 = plot(resistanceNum, PNUM(1)*resistanceNum + PNUM(2), 'Linewidth', 1.5);
ax = gca;
ax.FontName = 'LaTeX';
ax.TickLabelInterpreter = 'LaTeX';
ax.FontSize = 18;
ax.XColor = 'k';
ax.YColor = 'k';

ax.YLabel.String = ('Voltage (V)');
ax.YLabel.FontSize = 16;
ax.YLabel.FontWeight = 'bold';
ax.YLabel.Color = 'k';

ax.XLabel.String = ('Resistance (\Omega)');
ax.XLabel.FontSize = 16;
ax.XLabel.FontWeight = 'bold';
ax.XLabel.Color = 'k';

ax.Box = 'off';
ax.LineWidth = 1.5;
ax.YGrid = 'on';
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';

t = title('Thermal Voltage Vs. Resistance 1-9');
t.Color = 'k';
t.Interpreter = 'LaTeX';
t.FontSize = 24;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Gain vs Frequency%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f3 = figure(3);
b1 = plot(FrequencyHz, (MeanVV/VinScaled).^2);
b1.LineWidth = 2;
ax = gca;
ax.FontName = 'LaTeX';
ax.TickLabelInterpreter = 'LaTeX';
ax.FontSize = 18;
ax.XColor = 'k';
ax.YColor = 'k';

ax.YLabel.String = ('Gain (V/V)');
ax.YLabel.FontSize = 16;
ax.YLabel.FontWeight = 'bold';
ax.YLabel.Color = 'k';

ax.XLabel.String = ('Frequency (Hz)');
ax.XLabel.FontSize = 16;
ax.XLabel.FontWeight = 'bold';
ax.XLabel.Color = 'k';


ax.Box = 'off';
ax.LineWidth = 1.5;
ax.YGrid = 'on';
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
t = title('Gain vs Frequency');
t.Color = 'k';
t.Interpreter = 'LaTeX';
t.FontSize = 24;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Bode Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f4 = figure(4);
b2 = plot(log10(FrequencyHz), (MeanVV/VinScaled).^2);
b2.LineWidth = 2;
ax = gca;
ax.FontName = 'LaTeX';
ax.TickLabelInterpreter = 'LaTeX';
ax.FontSize = 18;
ax.XColor = 'k';
ax.YColor = 'k';

ax.YLabel.String = ('Gain (V/V)');
ax.YLabel.FontSize = 16;
ax.YLabel.FontWeight = 'bold';
ax.YLabel.Color = 'k';

ax.XLabel.String = ('Log_{10}(Frequency) (Hz)');
ax.XLabel.FontSize = 16;
ax.XLabel.FontWeight = 'bold';
ax.XLabel.Color = 'k';

ax.Box = 'off';
ax.LineWidth = 1.5;
ax.YGrid = 'on';
ax.XGrid = 'on';
ax.XMinorGrid ='on';
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
t = title('Bode Plot');
t.Color = 'k';
t.Interpreter = 'LaTeX';
t.FontSize = 24;

%[fits, gof] = createFits(resistanceNum, squareVoltageNum, weightsNum, resistanceAJ, squareVoltageAJ, weightsAJ); 

