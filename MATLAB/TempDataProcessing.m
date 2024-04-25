% Team 11 Data Processing File
% Contributers: Naomi Horiguchi '26, Tali Smith '26

% TempDataProcessing.m
% This script reads Teensy-collected data from an SD card and processes it.

clear;
clf;

%% Pull correct .bin and .inf files
filenum = '021'; % file number for the data you want to read
infofile = strcat('INF', filenum, '.TXT');
datafile = strcat('LOG', filenum, '.BIN');

%% Map from datatype to length in bytes
dataSizes.('float') = 4;
dataSizes.('ulong') = 4;
dataSizes.('int') = 4;
dataSizes.('int32') = 4;
dataSizes.('uint8') = 1;
dataSizes.('uint16') = 2;
dataSizes.('char') = 1;
dataSizes.('bool') = 1;

%% Read from info file to get log file structure
fileID = fopen(infofile);
items = textscan(fileID,'%s','Delimiter',',','EndOfLine','\r\n');
fclose(fileID);
[ncols,~] = size(items{1});
ncols = ncols/2;
varNames = items{1}(1:ncols)';
varTypes = items{1}(ncols+1:end)';
varLengths = zeros(size(varTypes));
colLength = 256;
for i = 1:numel(varTypes)
    varLengths(i) = dataSizes.(varTypes{i});
end
R = cell(1,numel(varNames));

%% Read column-by-column from datafile
fid = fopen(datafile,'rb');
for i=1:numel(varTypes)
    %# seek to the first field of the first record
    fseek(fid, sum(varLengths(1:i-1)), 'bof');
    
    %# % read column with specified format, skipping required number of bytes
    R{i} = fread(fid, Inf, ['*' varTypes{i}], colLength-varLengths(i));
    eval(strcat(varNames{i},'=','R{',num2str(i),'};'));
end
fclose(fid);

%% Thermistor data processing

% Convert thermistor voltage data to temperatures using calibration curve
vTherm = vpa(A01)*3.3/1024;
tempTherm = zeros([length(vTherm),1]);

for i = 1:length(vTherm)
    tempTherm(i) = -3.98*(vTherm(i)) + 21.4;
end

%% Thermocouple data processing

% NIST table (temperature in Celsius, thermoelectric voltage in mV)
tempNIST = 0:1:30;
mvNIST = [0.000, 0.059, 0.118, 0.176, 0.235, 0.294, 0.354, 0.413, 0.472, 0.532, 0.591,...
          0.651, 0.711, 0.770, 0.830, 0.890, 0.950, 1.010, 1.071, 1.131, 1.192,...
          1.252, 1.313, 1.373, 1.434, 1.495, 1.556, 1.617, 1.678, 1.740, 1.801];

NIST = [tempNIST; mvNIST]';

% Get reference temp from IC sensor
vIC = vpa(A03)*3.3/1024;
tempIC = zeros([1,length(vIC)]);

for i = 1:length(vIC)
    tempIC(i) = vIC(i)*50.9-22.1;
end

%% UNCOMMENT IF YOU WANT TO ACCOUNT FOR THE IC SENSOR

% Find the thermoelectric voltage of that temperature
% icNIST = zeros([1,length(vIC)]);
% counter = 1;
% for i = 1:length(vIC)
%     icNIST(i) = NIST(round(tempIC(i))+1,2);
% end

%% MORE THERMOCOUPLE DATA PROCESSING

% Read in data (vector of voltages)

vTcouple = double(A00)*3.3/1024;

% Get mV value read from the thermocouple to put into the NIST table
mvTcouple = zeros([1,length(vIC)]);
for i = 1:length(vTcouple)
    mvTcouple(i) = ((vTcouple(i) - 2.5)/2001)*1000;
end

icNIST = 1.373; % COMMENT OUT IF YOU ARE USING THE IC SENSOR PROCESSING LOOP ABOVE

mvTot = mvTcouple + icNIST;

% Find which thermoelectric voltage our measured value is closest to in the NIST table

tempThermocouple = zeros([1,length(vIC)]);
for i = 1:length(vIC)
    tempThermocouple(i) = 16.7*mvTot(i) + 0.111;
end

%% LEGACY CODE, DO NOT UNCOMMENT

% index = zeros([1,length(vIC)]);
% for i = 1:length(mvTot)
%     [m,k] = min(abs(mvNIST-mvTot(i)));
%     index(i) = k;
% end
% 
% tempThermocouple = zeros([1,length(vIC)]);
% for i = 1:length(index)
%     tempThermocouple(i) = NIST(index(i),1);
% end

%% Pressure sensor data processing

% Get depth from pressure
teensyPressure = vpa(A02);
depth = zeros([length(teensyPressure),1]);

for i = 1:length(teensyPressure)
    depth(i) = 0.0103*(teensyPressure(i)) - 0.889 - 0.61;
end

%% Extra processing

pwmIndex = find(~motorC)';

for i = 1:(((length(pwmIndex)-mod(length(pwmIndex), 30))/30))
    avThermocouple(i) = mean(tempThermocouple(pwmIndex(((i-1)*30)+1):pwmIndex(i*30)));
    avThermistor(i) = mean(tempTherm(pwmIndex(((i-1)*30)+1):pwmIndex(i*30)));
    sampMedian(i) = round(median(pwmIndex(((i-1)*30)+1):pwmIndex(i*30)));
end

%% Plotting 
sampNumVec = 1:length(A00);

% %Thermistor plots
% figure(1)
% 
% subplot(2,1,1);
% plot(sampNumVec,tempTherm);
% xlabel("Sample number");
% ylabel("Temperature (C)");
% title("Temperature (C) vs. Sample number for Thermistor")
% 
% subplot(2,1,2);
% plot(vTherm,tempTherm);
% xlabel("Voltage (V)");
% ylabel("Temperature (C)")
% title("Temperature (C) vs. Voltage (V) for Thermistor")
% 
% % Pressure plots
% figure(2)
% 
% subplot(2,1,1)
% plot(teensyPressure,depth)
% xlabel("Teensy Units")
% ylabel("Depth (m)")
% title("Depth (m) vs. Teensy Units")
% 
% subplot(2,1,2)
% plot(sampNumVec,depth)
% xlabel("Sample number")
% ylabel("Depth (m)")
% title("Depth (m) vs. Sample number")
% 
% % IC vs. Thermocouple (Box temp. vs outside temp.)
% figure(3)
% 
% plot(sampNumVec, tempIC)
% xlabel("Sample number")
% ylabel("Temperature (C)")
% title("Temperature (C) vs. Sample number for IC")
% 
% % Thermistor vs. Thermocouple
% figure(4)
% 
% plot(sampNumVec,tempTherm)
% hold on
% plot(sampNumVec,tempThermocouple)
% hold on
% xlabel("Sample Number")
% ylabel("Temperature (C)")
% title("Thermistor vs. Thermocouple")

%% Plots for final report

% Motor interference (File 014)
figure(5)

yyaxis left
plot(sampNumVec,tempThermocouple,'LineWidth',2)
xlabel("Sample Number")
ylabel("Temperature (C)")
yyaxis right
plot(sampNumVec, motorC,'LineWidth',3)
ylabel("PWM")
title("Interference from Vertical Motor")

% Depth/Thermistor Temperature Overlay (File 021)
figure(6)

yyaxis left
plot(sampNumVec,tempTherm,'LineWidth',2)
xlabel("Sample Number")
ylabel("Temperature (C)")
yyaxis right
plot(sampNumVec, -depth,'LineWidth',2)
ylabel("Depth (m)")
title("Thermistor Change in Temperature over Time vs. Depth over Time")
xlim([0, length(vIC)])

% Thermistor vs. Depth (File 021)
figure(7)

p = polyfit(depth, tempTherm, 1);
px = [min(depth) max(depth)];
py = polyval(p, px);
scatter(depth,tempTherm,14,"filled")
hold on
plot(px, py, 'LineWidth', 2)
xlabel("Depth (m)")
ylabel("Temperature(C)")
title("Thermistor temperature (C) vs. Depth (m)")
xlim([0, max(depth)])

% Sinking behavior for File 021, Bobbing behavior for file 020 (CHANGE FILE # AND PLOT TITLE TO GO FROM ONE TO THE OTHER)
figure(8)
plot(sampNumVec,depth)
xlabel("Sample Number")
ylabel("Depth (m)")
title("Depth (m) over Time: Sinking Behavior")

% Thermocouple vs. Thermistor Discrete Points w/ Thermocouple with no Interference
figure (9)
scatter(depth(sampMedian), avThermocouple,"filled")
hold on
scatter(depth(sampMedian), avThermistor,"filled")
xlabel("Depth (m)")
ylabel("Temperature (C)")
title("Temperature (C) vs. Depth (m) for Thermocouple")
