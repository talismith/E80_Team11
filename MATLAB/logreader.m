% TempDataProcessing.m
% This script reads Teensy-collected data from an SD card and processes it.

clear;
clf;

%% Pull correct .bin and .inf files
filenum = '038'; % file number for the data you want to read
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
thermTemp = zeros([length(VARHERE),1]);
for counter = 1:length(VARHERE)
    thermTemp(counter) = -3.98*(VARHERE(counter)) + 21.4;
end

%% Thermocouple data processing

% NIST table (temperature in Celsius, thermoelectric voltage in mV)
tempNIST = 0:1:30;
mvNIST = [0.000  0.059  0.118  0.176  0.235  0.294  0.354  0.413  0.472 0.532  0.591...
  0.651  0.711  0.770  0.830  0.890  0.950  1.010  1.071  1.131  1.192...
  1.252  1.313  1.373  1.434  1.495  1.556  1.617  1.678  1.740  1.801];

NIST = [tempNIST; mvNIST]';

% Get reference temp from IC sensor
vIC = [0.860, 0.860]; %NEED TO INSERT CODE TO READ IN
tIC = zeros([1,length(vIC)]);
for i = 1:length(vIC)
    tIC(i) = (vIC(i) - 0.4)/0.0195;
end

% Find the thermoelectric voltage of that temperature
icNIST = zeros([1,length(vIC)]);
for i = 1:length(tIC)
    icNIST(i) = NIST(round(tIC(i))+1,2);
end

% Read in data (vector of voltages)
vTcouple = [1.473,0.12]; %NEED TO INSERT CODE TO READ IN

% Get mV value read from the thermocouple to put into the NIST table
mvTcouple = zeros([1,length(vIC)]);
for i = 1:length(vTcouple)
    mvTcouple(i) = (vTcouple(i) - 2.5)/2001*1000;
end

mvTot = mvTcouple + icNIST;

% Find which thermoelectric voltage our measured value is closest to in the NIST table
index = zeros([1,length(vIC)]);
for i = 1:length(mvTot)
    [m,k] = min(abs(mvNIST-mvTot(i)));
    index(i) = k;
end

thermocoupleTemp = zeros([1,length(vIC)]);
for i = 1:length(index)
    thermocoupleTemp(i) = NIST(index(i),1);
end

%% Plotting
sampNumVec = 1:length(accelX);

% Thermistor plots
figure 1;

subplot(3,1,1);
plot(sampNumVec,thermTemp);
xlabel("Sample number");
ylabel("Temperature (C)");

subplot(3,1,2);
plot(thermVolt,thermTemp);
xlabel("Voltage (UNIT)");
ylabel("Temperature (C)");

subplot(3,1,3);
plot(timeVec,thermTemp);

% Pressure plots
figure 2;

subplot(3,1,1);
plot(pressVolt,pressAmp);
xlabel("Voltage (UNIT)");
ylabel("Pressure (UNIT)");

subplot(3,1,2);
plot(pressVolt, depth);
xlabel("Voltage (UNIT)");
ylabel("Depth (UNIT)");

subplot(3,1,3);
plot(pressAmp,depth);
xlabel("Pressure (UNITS)");
ylabel("Depth (UNITS)")

% Thermistor vs. Thermocouple
figure 3;

plot(sampNumVec,thermTemp);
hold on;
plot(sampNumVec,thermocoupleTemp);
xlabel("Sample Number");
ylabel("Temperature (C)")

plot(depth,thermTemp);
hold on;
plot(depth,thermocoupleTemp);
xlabel("Depth (UNITS)");
ylabel("Temperature (C)")

