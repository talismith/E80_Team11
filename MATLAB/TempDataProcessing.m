% TempDataProcessing.m
% This script reads Teensy-collected data from an SD card and processes it.

clear;
clf;

%% Pull correct .bin and .inf files
filenum = '004'; % file number for the data you want to read
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
    tempIC(i) = (vIC(i) - 0.4)/0.0195;
end

% Find the thermoelectric voltage of that temperature
icNIST = zeros([1,length(vIC)]);
counter = 1;
for i = 1:length(vIC)
    icNIST(i) = NIST(round(tempIC(i))+1,2);
end

% CAUTION BELOW HERE

% Read in data (vector of voltages)

vTcouple = vpa(A00)*3.3/1024;

% Get mV value read from the thermocouple to put into the NIST table
mvTcouple = zeros([1,length(vIC)]);
for i = 1:length(vTcouple)
    mvTcouple(i) = ((vTcouple(i) - 2.5)/2001)*1000;
end

mvTot = mvTcouple + icNIST; %+1 is to adjust for incorrect linear sensor, take out next time

% Find which thermoelectric voltage our measured value is closest to in the NIST table
index = zeros([1,length(vIC)]);
for i = 1:length(mvTot)
    [m,k] = min(abs(mvNIST-mvTot(i)));
    index(i) = k;
end

tempThermocouple = zeros([1,length(vIC)]);
for i = 1:length(index)
    tempThermocouple(i) = NIST(index(i),1);
end

% CAUTION ABOVE HERE

%% Pressure sensor data processing

% Get depth from pressure
teensyPressure = vpa(A02);
depth = zeros([length(teensyPressure),1]);

for i = 1:length(teensyPressure)
    depth(i) = 0.0103*(teensyPressure(i)) - 0.889;
end

%% Plotting
sampNumVec = 1:length(A00);

%Thermistor plots
figure(1)

subplot(2,1,1);
plot(sampNumVec,tempTherm);
xlabel("Sample number");
ylabel("Temperature (C)");
title("Temperature (C) vs. Sample number for Thermistor")

subplot(2,1,2);
plot(vTherm,tempTherm);
xlabel("Voltage (V)");
ylabel("Temperature (C)")
title("Temperature (C) vs. Voltage (V) for Thermistor")

% Pressure plots
figure(2)

subplot(2,1,1)
plot(teensyPressure,depth)
xlabel("Teensy Units")
ylabel("Depth (m)")
title("Depth (m) vs. Teensy Units")

subplot(2,1,2)
plot(sampNumVec,depth)
xlabel("Sample number")
ylabel("Depth (m)")
title("Depth (m) vs. Sample number")

% IC vs. Thermocouple (Box temp. vs outside temp.)
figure(3)

plot(sampNumVec, tempIC)
xlabel("Sample number")
ylabel("Temperature (C)")
title("Temperature (C) vs. Sample number for IC")

% Thermistor vs. Thermocouple
figure(4)

plot(sampNumVec,tempTherm)
hold on
plot(sampNumVec,tempThermocouple)
hold on
xlabel("Sample Number")
ylabel("Temperature (C)")
title("Thermistor vs. Thermocouple")


