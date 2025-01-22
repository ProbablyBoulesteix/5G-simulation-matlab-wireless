%NOTE: based on https://nl.mathworks.com/help/phased/ug/sinr-map-for-a-5g-urban-micro-cell-test-environment.html
close all;
clear;
maxNumCompThreads("automatic"); %set as desired depending on cpu available, or set as "automatic"

%% globals params

name_of_saved_filestate = "run_10dbm500m160kmh-80dbmInt.mat";


%inter-site-distance
isd = 500; %m by default 200m or 250m is fine
cell_range_factor = 1.0; %range of cell wrt isd. factor of 1.0 assumes no overlap (best case)
sinr_map_resolution = 30; %fraction of isd, higher = better res -> CAUTION: computationnaly demanding, default was 20

%transmiter info
tx_power = 10; %dbm
fq = 4e9; %Hz, carrier frequency, 4ghz is default as defined in the matlab example

%reciever info
guardband_size = 900e3;  %see http://howltestuffworks.blogspot.com/2019/11/5g-nr-resource-blocks.html or TR 38.101 -> each set ressource blocks has preset gardbands on each side depending on bw and subcarrier spacing
rx_bw = 40e6; %reciever bandwidth, Hz
rx_bw = rx_bw - 2*guardband_size;  %remove useless/guardband from each side of bandwidth

rx_noiseF = 7; %dbm, noise figure of reciever-> RX antenna gain is defined lowerin RX section at 0db (unity gain)

ambient_interference = -80; %dbm -> adjustable, can be used to simulated baseline network load through interference SHOULD BE MUCH LOWER THAN TX POWER
%pathloss_model_type = "longley-rice"; %refer to https://nl.mathworks.com/help/comm/ref/rfprop.freespace.pathloss.html -> not implemented for now

%moving info
moving_speed = 160 / 3.6; %speed of vehicle on route, m/s

%interpolation position subdivision
%CAUTION: computation time increases rapidly with decreasing length here,
%1-5m is generally good
position_subdivide_length = 5; %m -> trajectory is divided into simulated steps spaced this far apart


%symbol info
slot_period = 1e-3; %time per ofdm slot (symbols + prefixes)
symbols_per_slot = 14; % default is 14. ofdm symbols per slot period

symbol_period = slot_period / symbols_per_slot;
symbol_rate =  1/ symbol_period;

%modulator params
modulation_order = 2; % width of subcarrier = 2^MO * 15khz -> MO 0-2 for <6ghz, 2-5 for mmwave, see https://www.keysight.com/us/en/assets/9921-03326/training-materials/Understanding-the-5G-NR-Physical-Layer.pdf 

%additionnal params
spectral_eff_type = 0; %if 0, uses 90% BLER, 99.9% if 1 ->currently not implemented, 0.1 is used


%% CELL INFORMATION 
% Define center location site (cells 1-3)
centerSite = txsite('Name','MathWorks Glasgow', ...
    'Latitude',55.862787,...
    'Longitude',-4.258523);

% Initialize arrays for distance and angle from center location to each cell site, where
% each site has 3 cells
numCellSites = 19;
siteDistances = zeros(1,numCellSites);
siteAngles = zeros(1,numCellSites);

c = physconst('LightSpeed');

%find number of cells per layer here
layer_counts = zeros(1);


% Define distance and angle for inner ring of 6 sites (cells 4-21)

siteDistances(2:7) = isd;
siteAngles(2:7) = 30:60:360;

% Define distance and angle for middle ring of 6 sites (cells 22-39)
siteDistances(8:13) = 2*isd*cosd(30);
siteAngles(8:13) = 0:60:300;

% Define distance and angle for outer ring of 6 sites (cells 40-57)
siteDistances(14:19) = 2*isd;
siteAngles(14:19) = 30:60:360;

%siteDistances(20:25) = 3*isd;
%siteAngles(20:25) = 0:60:300;

%siteDistances(26:31) = 3*isd*cosd(30);
%siteAngles(26:31) = 30:60:360;



%% CELL PARAMS
% Initialize arrays for cell transmitter parameters
numCells = numCellSites*3;
cellLats = zeros(1,numCells);
cellLons = zeros(1,numCells);
cellNames = strings(1,numCells);
cellAngles = zeros(1,numCells);


% Define cell sector angles
cellSectorAngles = [30 150 270];

% For each cell site location, populate data for each cell transmitter
cellInd = 1;
for siteInd = 1:numCellSites
    % Compute site location using distance and angle from center site
    [cellLat,cellLon] = location(centerSite, siteDistances(siteInd), siteAngles(siteInd));
    
    % Assign values for each cell
    for cellSectorAngle = cellSectorAngles
        cellNames(cellInd) = "Cell " + cellInd;
        cellLats(cellInd) = cellLat;
        cellLons(cellInd) = cellLon;
        cellAngles(cellInd) = cellSectorAngle;
        cellInd = cellInd + 1;
    end
end

%% TX SITES
% Define transmitter parameters using Table 8-2 (b) of Report ITU-R M.[IMT-2020.EVAL]
fq = fq; % Carrier frequency (4 GHz) for Dense Urban-eMBB
antHeight = 40; % m
txPowerDBm = tx_power; % Total transmit power in dBm
txPower = 10.^((txPowerDBm-30)/10); % Convert dBm to W

% Create cell transmitter sites
txs = txsite('Name',cellNames, ...
    'Latitude',cellLats, ...
    'Longitude',cellLons, ...
    'AntennaAngle',cellAngles, ...
    'AntennaHeight',antHeight, ...
    'TransmitterFrequency',fq, ...
    'TransmitterPower',txPower);

% Launch Site Viewer
%viewer = siteviewer;

% Show sites on a map
%show(txs);
%viewer.Basemap = 'topographic';

%% ANTENNA INFO 
% Define pattern parameters
azvec = -180:180;
elvec = -90:90;
Am = 30; % Maximum attenuation (dB)
tilt = 0; % Tilt angle
az3dB = 65; % 3 dB bandwidth in azimuth
el3dB = 65; % 3 dB bandwidth in elevation

% Define antenna pattern
[az,el] = meshgrid(azvec,elvec);
azMagPattern = -12*(az/az3dB).^2;
elMagPattern = -12*((el-tilt)/el3dB).^2;
combinedMagPattern = azMagPattern + elMagPattern;
combinedMagPattern(combinedMagPattern<-Am) = -Am; % Saturate at max attenuation
phasepattern = zeros(size(combinedMagPattern));

% Create antenna element
antennaElement = phased.CustomAntennaElement(...
    'AzimuthAngles',azvec, ...
    'ElevationAngles',elvec, ...
    'MagnitudePattern',combinedMagPattern, ...
    'PhasePattern',phasepattern);
   
% Display radiation pattern
f = figure;
pattern(antennaElement,fq);

%% RX 
% Assign the antenna element for each cell transmitter
for tx = txs
    tx.Antenna = antennaElement;
end

% Define receiver parameters using Table 8-2 (b) of Report ITU-R M.[IMT-2020.EVAL] 
bw = rx_bw; % 20 MHz bandwidth
rxNoiseFigure = rx_noiseF; % 7 dB in spec
rxNoisePower = -174 + 10*log10(bw) + rxNoiseFigure
rxNoisePowerInterference = sum_dbm_power(rxNoisePower, ambient_interference)
rxGain = 0; % dBi
rxAntennaHeight = 1.5; % m

% Display SINR map
if isvalid(f)
    close(f)
end
sinrValues = sinr(txs,'freespace', ...
    'ReceiverGain',rxGain, ...
    'ReceiverAntennaHeight',rxAntennaHeight, ...
    'ReceiverNoisePower',rxNoisePower, ...    
    'MaxRange',cell_range_factor*isd, ...
    'Resolution',isd/sinr_map_resolution);

datalength = height(sinrValues.Data)

%uncomment if want to display SINR map. Takes time!

%sinr(txs,'freespace', 'ReceiverGain',rxGain, 'ReceiverAntennaHeight',rxAntennaHeight, 'ReceiverNoisePower',rxNoisePower, ...    
%    'MaxRange',cell_range_factor*isd, ...
%    'Resolution',isd/sinr_map_resolution)
%% now define MIMO
% Define array size
nrow = 8;
ncol = 8;

% Define element spacing
lambda = physconst('lightspeed')/fq;
drow = lambda/2;
dcol = lambda/2;

% Define taper to reduce sidelobes 
dBdown = 30;
taperz = chebwin(nrow,dBdown);
tapery = chebwin(ncol,dBdown);
tap = taperz*tapery.'; % Multiply vector tapers to get 8-by-8 taper values

% Create 8-by-8 antenna array
cellAntenna = phased.URA('Size',[nrow ncol], ...
    'Element',antennaElement, ...
    'ElementSpacing',[drow dcol], ...
    'Taper',tap, ...
    'ArrayNormal','x');
    
% Display radiation pattern
f = figure;
pattern(cellAntenna,fq);
% Assign the antenna array for each cell transmitter, and apply downtilt.
% Without downtilt, pattern is too narrow for transmitter vicinity.
downtilt = 15;
for tx = txs
    tx.Antenna = cellAntenna;
    tx.AntennaAngle = [tx.AntennaAngle; -downtilt];
end

% Display SINR map
if isvalid(f)
    close(f)
end
sinrValues = sinr(txs,'freespace', ...
    'ReceiverGain',rxGain, ...
    'ReceiverAntennaHeight',rxAntennaHeight, ...
    'ReceiverNoisePower',rxNoisePower, ...    
    'MaxRange',cell_range_factor*isd, ...
    'Resolution',isd/sinr_map_resolution)

%uncomment if want display
% sinr(txs,'freespace', ...
%     'ReceiverGain',rxGain, ...
%     'ReceiverAntennaHeight',rxAntennaHeight, ...
%     'ReceiverNoisePower',rxNoisePower, ...    
%     'MaxRange',cell_range_factor*isd, ...
%     'Resolution',isd/sinr_map_resolution)

datalength = height(sinrValues.Data)


%%

% define trajectory
%critical waypoints go here (turns) in lattitude/longitude degree/decimal
%format
waypoints = [55.859444, -4.233611; 55.860833, -4.245278; 55.861389, -4.25; 55.864444, -4.248889; 55.865278, -4.256944; 55.860833, -4.258333;  55.861389,  -4.263056; 55.864444, -4.261944;  55.865278, -4.270278; 55.867222, -4.271111; 55.869444,-4.267778];%; 55.871389, -4.27 ;55.874722, -4.279444]

lats = zeros(1, length(waypoints));
longs = zeros(1, length(waypoints));

for k = 1:1:length(waypoints)
    lats(k) = waypoints(k, 1);
    longs(k) = waypoints(k, 2);
end

%h = geoplot(lats, longs);

%distance between each point specified in waypoint using great circle
%method (used later for speed calc)
wgs84 = wgs84Ellipsoid("m");
distances = zeros(1, length(waypoints)-1);
running_distance = 0; %total distance travelled during trajectory
for k = 1:1:length(waypoints)-1
    arclen = distance(lats(k),longs(k),lats(k+1),longs(k+1),wgs84);
    running_distance = running_distance + arclen;
    distances(k) = arclen;
end

%let's now create new waypoints closer together between the specified
%waypoints
subdivide_length = position_subdivide_length; 
interpolated_lats = zeros(0);
interpolated_longs = zeros(0);



for k = 1:1:length(waypoints) -1
    lat1 = lats(k);
    lat2 = lats(k+1);
    long1 = longs(k);
    long2 = longs(k+1);

    lats_i = [lat1; lat2];
    longs_i = [long1; long2];
    
    dist = distances(k);
    
    chunk_count = floor(dist/subdivide_length);

    interpolated_lats = [interpolated_lats; lat1];
    interpolated_longs = [interpolated_longs; long1];

    for i = 1:1:chunk_count
        newlong_fact = (long2 - long1)/chunk_count;
        newlong = long1 + newlong_fact * i;
        newlat = intrplat(longs_i, lats_i, newlong, 'pchip', 'degrees');
        interpolated_longs = [interpolated_longs; newlong];
        interpolated_lats = [interpolated_lats; newlat];
        
    end
    interpolated_longs = [interpolated_longs; long2];
    interpolated_lats = [interpolated_lats; lat2];
end

geoplot(interpolated_lats, interpolated_longs, lats, longs, 'p');


%% Driving trajectory and SINR MAPPING
%for each waypoint, let's find the base SINR at the location

base_velocity = moving_speed; %magnitude of the UE velocity vector , m/s

N = numel(interpolated_longs);
SINR_at_waypoint = zeros(1,N);
cellID_at_waypoint = zeros(1,numCells);
angle_wrt_cell_at_waypoint= zeros(1,numCells);
radial_velocity_wrt_cell = zeros(1,numCells);

%get base data from map
SINR_lats = sinrValues.Data.Latitude;
SINR_longs = sinrValues.Data.Longitude;
SINR_coords = [SINR_lats, SINR_longs];
SINR_DB = sinrValues.Data.SINR; % in DB

cell_coords = [cellLats; cellLons].';

%for each waypoint, find the closest sinr datapoint and note measured sinr ->
%simple and stupid algorithm but works if dataset is large enough (and
%often is !)
% also take note of closet cell ID as well as the angle wrt that cell for
% later radial speed calculation
for k = 1:1:N
    %first find sinr
    point_query = [interpolated_lats(k), interpolated_longs(k)];
    [n, dist] = dsearchn(SINR_coords, point_query);
    
    SINR_at_waypoint(k) = SINR_DB(n);

    %compute cell index and distance to tell
    [n_cell, d_cell] = dsearchn(cell_coords, point_query);
    cellID_at_waypoint(k) = n_cell;
    
    %tower = txsite("Latitude",cell_coords(n_cell, 1), "Longitude", cell_coords(n_cell, 2),"Antenna", antennaElement ,TransmitterFrequency=fq);
    %site = rxsite("Latitude",point_query(1), "Longitude",point_query(2));
    %pm = propagationModel("freespace");
    %pl = pathloss(pm,site, tower);


    %now find angle of motion wrt cell by computing vector between current
    %and next point
    link_vector = [cellLats(n_cell) - interpolated_lats(k), cellLons(n_cell)- interpolated_longs(k)];
    if k < N
        mov_vec = [interpolated_lats(k+1)-interpolated_lats(k); interpolated_longs(k+1) - interpolated_longs]; 
    else
        mov_vec = [0, 0];
    end
    ang = atan2d(link_vector(1)*mov_vec(2) - mov_vec(1)*link_vector(2), mov_vec(1)*mov_vec(2)+link_vector(1)*link_vector(2));
    angle_wrt_cell_at_waypoint(k) = ang;
    %now finally compute the radial velocity of user wrt to the nearest
    %cell
    radial_speed = base_velocity * cosd(ang);
    radial_velocity_wrt_cell(k) = radial_speed;


end


d_vec = zeros(1, numel(interpolated_longs));
d_run = zeros(1, numel(interpolated_longs));
for k = 1:1:numel(d_vec)
    if (k < numel(d_vec))
        arclen = distance(interpolated_lats(k),interpolated_longs(k),interpolated_lats(k+1),interpolated_longs(k+1),wgs84);
        d_run(k+1) = d_run(k) + arclen;
    else
        arclen = 0;
    end
    d_vec(k) = arclen;
    
end

%plot(d_vec, SINR_at_waypoint)
%plot(d_run, angle_wrt_cell_at_waypoint)
%plot(d_vec, radial_velocity_wrt_cell)
%geoplot(interpolated_lats, interpolated_longs, lats, longs, 'p') %plot trajectory

%% mapping SINR to code rate and modualtion format

%CAUTION: MCS/CPI table below maps SNR/SINR values to specific modulation
%types and coding rates. It seems to be vendor/manufacturer specific, so
%may also be able to use
%https://www.telecomhall.net/t/mapping-of-sinr-to-mcs/28944 but unsourced
%so i'd rather not for now

%we'll be using the following MCS/QPI table https://www.researchgate.net/figure/Modulation-and-Coding-Schemes-MCS-for-eMBB-and-URLLC-services-with-different-BLERs_tbl1_339635042

% note on modulation: BPSK = 2; QPSK = 4, 16QAM = 16; 64QAM = 64
%threshold SINR
SINR_thresholds_BLER01 = [-6.5, -4.0, -2.6, -1.0, 1.0, 3.0, 6.6, 10.0, 11.4,11.8, 13.0, 13.8, 15.6, 16.8, 17.6];
SINR_thresholds_BLER0001 = SINR_thresholds_BLER01 + 4;

%secon threshold table from
%%https://www.telecomhall.net/uploads/db2683/original/3X/c/f/cfb9dee6c3ae343cd13b9f480a24d36ec662c4f6.jpeg
%%(unsourced)
SINR_thresholds_beta = [-6.7, -4.7,-2.3,0.2,2.4,4.3,5.9,7.2,8.7,10.3,11.7,13.1,14.3,15.8,17.3,18.7,20.0,21.4,24.0,25.3,26.5,27.6,28.7,29.8,30.9,13.9,32.9,33.9];

MOD_thresholds_beta = [4,4,4,4,4,4,4,16,16,16,16,16,64,64,64,64,64,256,256,256,256,256,256,256,256,256,256,256,256];
CR_thresholds_beta = [0,1172,0.1885,0.3008,0.4385,0.5879,0.7402,0.8809,0.332,0.4385,0.5537,0.667,0.7539,0.6016,0.7402,0.8477,0.9258,0.9639,0.7734,0.8701,0.916,0.9482,0.96,0.978,0.9834,0.9869,0.9899,0.9935,0.997,0.9989];


%CODING RATES FOR EACH THRESHOLD, 
CR_thresholds = [1/12, 1/9, 1/6, 1/3, 1/2, 3/5, 1/3, 1/2, 3/5, 1/2, 1/2, 3/5, 3/4, 5/6, 11/12 ];

%spectral efficiency (bits/symb) for each threshold
EFF_thresholds =  [0.15, 0.23, 0.38, 0.6, 0.88, 1.18, 1.48, 1.91, 2.41, 2.73, 3.32, 3.9, 4.52, 5.12, 5.55];
%Modualtion type at each threshold
MOD_thresholds = [4,4,4,4,4,16,16,16,64,64,64,64,64,64];

%fixme: since a lot of the code below just used the BLER01 database and
%names it explicitely, let's just copy to it instead rather than renameing
%everything
SINR_thresholds_BLER01 = SINR_thresholds_beta - 10;
MOD_thresholds = MOD_thresholds_beta;
CR_thresholds = CR_thresholds_beta;


%arrays storing index into MCS table, coding rate, modulation goes here
CR_at_waypoints_BLER01 = zeros(1, length(interpolated_longs));
CR_at_waypoints_BLER0001 = zeros(1, length(interpolated_longs));

MOD_at_waypoints_BLER01 = zeros(1, length(interpolated_longs));
MOD_at_waypoints_BLER0001 = zeros(1, length(interpolated_longs));

EFF_at_waypoints_BLER01 = zeros(1, length(interpolated_longs));
EFF_at_waypoints_BLER0001 = zeros(1, length(interpolated_longs));





for k = 1:1:length(SINR_at_waypoint)
    sinr_wp = SINR_at_waypoint(k);
    %for each datapoint, find index into MCS tables
    idx_BLER01 = 1;
    idx_BLER0001 = 1;
    while ((sinr_wp > SINR_thresholds_BLER01(idx_BLER01))&&(idx_BLER01 < numel(SINR_thresholds_BLER01)))
    %check: if the current sinr at the waypoint is above the threshold sinr
    %in the table, we can look at the next index to see if we exceed that
    %target too
        idx_BLER01 = idx_BLER01 + 1; 
    end
    idx_BLER01 = idx_BLER01 - 1; %take note of last valid index
    
    %do the same for the 0.1% BLR
    while (sinr_wp > SINR_thresholds_BLER0001(idx_BLER0001)) && (idx_BLER0001 < numel(SINR_thresholds_BLER0001))
        idx_BLER0001 = idx_BLER0001 + 1; 
    end
    idx_BLER0001 = idx_BLER0001 - 1; %take note of last valid index

    %check to see if indexes are 0 -> out of range of tower/insufficient
    %SINR
    if idx_BLER01 == 0
        MOD_at_waypoints_BLER01(k) = 0;
        CR_at_waypoints_BLER01(k) = 0;
        %EFF_at_waypoints_BLER01(k) = 0;
    else
        %now find actual coding rates, modualtions and symbol rates for each
        MOD_at_waypoints_BLER01(k) = MOD_thresholds(idx_BLER01);
        CR_at_waypoints_BLER01(k) = CR_thresholds(idx_BLER01);
        %EFF_at_waypoints_BLER01(k) = EFF_thresholds(idx_BLER01);
    end
    if idx_BLER0001 == 0
        MOD_at_waypoints_BLER0001(k) = 0;
        CR_at_waypoints_BLER0001(k) = 0;
        %EFF_at_waypoints_BLER0001(k) = 0;
    else
        MOD_at_waypoints_BLER0001(k) = MOD_thresholds(idx_BLER0001);
        CR_at_waypoints_BLER0001(k) = CR_thresholds(idx_BLER0001);
        %EFF_at_waypoints_BLER0001(k) = EFF_thresholds(idx_BLER0001);
    end

end
% define OFDM format and compute per-subcarrier phase shift
% see https://www.keysight.com/us/en/assets/9921-03326/training-materials/Understanding-the-5G-NR-Physical-Layer.pdf
%-> modulation orders 0,1,2 available for sub-6ghz, MO 2-5 available for
%mmwave (max 480khz)
subcarrier_modulation_order = modulation_order;
subcarrier_width = (2^subcarrier_modulation_order) *  15e3;
subcarrier_channel_count = floor(bw/subcarrier_width); %assumption: we'll use the reciever bandwidth, channel aggregation handled implicitely

%central frequencies for each subcarrier
subcarrier_frequencies = fq - 0.5*bw : subcarrier_width: fq + 0.5*bw - subcarrier_width; 
% we now want to compute the phase shift associated with each subcarrier
% over each waypoint

%define symbol rate data here
% as seen https://www.sharetechnote.com/html/5G/5G_FrameStructure_Candidate.html
%assume 1ms/slot and 14symb/slot for now ->72µs per symbol
symbol_duration = symbol_period;



subcarrier_shifts = zeros(numel(interpolated_longs), numel(subcarrier_frequencies));
subcarrier_doppler = zeros(numel(interpolated_longs), numel(subcarrier_frequencies));
for i = 1:1:numel(interpolated_longs)
    
    radial_v = radial_velocity_wrt_cell(i);
   
    for j = 1:1:numel(subcarrier_frequencies)
        doppler_coef = (c + radial_v)/(c);
        lambda_subc = c/subcarrier_frequencies(j);
        shift = (2*pi*radial_v*symbol_duration)/lambda_subc; %NOTE: this is phase shift over symbol duration as defined above (eq. to average symbol phase shift)
        %note: shift is in radians
        subcarrier_shifts(i,j) = shift;
        %also compute frequency offset
        subcarrier_doppler(i,j) = doppler_coef(1) *  subcarrier_frequencies(j);
    end
end


%% Plot constellation diagrams and shifts here -< not necessary to run this unless you want a visualization

data = (0:63)'; %random data 
QAM64_cd_ref = qammod(data,64); %modulated symbol using 64QAM
md = 5 / max(abs(QAM64_cd_ref)) %normalize radius of constellation here
QAM64_cd_ref = QAM64_cd_ref * md;

pfo = comm.PhaseFrequencyOffset(PhaseOffset=5, FrequencyOffset=0);
QAM_shifted = pfo(QAM64_cd_ref);

%scatterplot(QAM64_cd_ref)
%scatterplot(QAM_shifted)

ref_re = real(QAM64_cd_ref);
ref_im = imag(QAM64_cd_ref);

s_re = real(QAM_shifted);
s_im = imag(QAM_shifted);


%% modulate baseband signal and plot  both original and phase+freq offset constellation diagram
%tgermal noise goes here
noise_ref_db = -204 + rxNoiseFigure; %N0 in dB
noise_ref = 10^(noise_ref_db/10); %convert to watts
interference_spectral_density = ambient_interference - 10*log10(bw); %spectral density of ambient interference

noise_interference_ref_db = sum_dbm_power(interference_spectral_density,noise_ref_db +30) - 30; % (dbm) assumes noise follows gaussian and can be imposed onto noise
noise_interference_ref = 10^(noise_interference_ref_db/10); %convert to watts




unshifed_pairwise_prob = ones(1,numel(interpolated_longs));
unshifed_pairwise_min_dist = ones(1,numel(interpolated_longs));
shifted_pairwise_error_prob = ones(numel(interpolated_longs), numel(subcarrier_frequencies));
shifted_pairwise_error_prob_freq = ones(numel(interpolated_longs), numel(subcarrier_frequencies));
BER_shift_freq_prob =  ones(numel(interpolated_longs), numel(subcarrier_frequencies));
BER_shift_prob = ones(numel(interpolated_longs), numel(subcarrier_frequencies));
BER_static_prob = ones(numel(interpolated_longs), numel(subcarrier_frequencies));



%now for each waypoint and each subcarrrier, compute the minimum pairwise
%distance and max pairwise probability of symbol mismatch
for i = 1:1:numel(interpolated_longs)
    %first: find modulation format that would be used at low speed for this waypoint
    %(assumption: mod stays the same even at higher speed -> not true but
    %we can't yet recompute SNR at higher speeds
    progress_percent = round(i * 100 / numel(interpolated_longs), 2) %for display purposes, this can get quite long so we want to show progress
    mod01 = MOD_at_waypoints_BLER01(i); %only focus on 0.1% metric for now
    
    %define a string of bits to transmit and generate constellation for
    %unmoving case: define modulator object and modulate random sample data
    if mod01 == 0
        %error/no connect case
        continue
    elseif mod01 == 4 %QPSK
        data = (0:mod01-1)';
        unshifted_mod = pskmod(data,mod01);
    elseif mod01 == 16 %16QAM
        data = (0:mod01-1)';
        unshifted_mod = qammod(data,mod01);
    elseif mod01 == 64 %64QAM
        data = (0:mod01-1)';
        unshifted_mod = qammod(data,mod01);
    elseif mod01 == 256 %256QAM
        data = (0:mod01-1)';
        unshifted_mod = qammod(data,mod01);
     
    end
    
    if mod01 ~= 0
        %first find bit energy as (signal power * bits_per_symbol) /  
        datarate_lin = symbol_rate * log2(mod01); %datarate in bits = symbol rates * bits per symbol
        signal_power = SINR_at_waypoint(i) + (sum_dbm_power(rxNoisePower, ambient_interference)); %S/NI => S_db - NI_db -> in dbm here
        energy_per_bit = signal_power- 30 - 10*log10(datarate_lin); %energy per bit in db
        energy_per_bit = 10^(energy_per_bit/10); %convert to joules
        symbol_energy = compute_energy_per_symbol(energy_per_bit, mod01); %commpute ES
        
        %normalize constellation radius to computed symbol energy
        normalized_constellation_radius =  (symbol_energy)^0.5; %note: "radius" here is defined as the magnitude of outward (~diagonal) symbol. for symbols on im or real axis, magnitude is sqrt of radius
        max_cd_mag = normalized_constellation_radius / max(real(unshifted_mod));
        unshifted_mod = unshifted_mod * (max_cd_mag);

        %compute FOR UNSHIFTED CASE the minimum inter-symbol distance and
        %associated maximal pairwise probability
        [d_mins, p1s, p2s] = find_smallest_distance_complex(unshifted_mod);
        %note: first element is always samewise match, should be disgarded
        %by later function computing mismatch probability
        

    
        %SER and BER computations for static case
        prob = prob_overreach(d_mins, noise_interference_ref); %probability that any symbol is mistaken for another 
        unshifed_pairwise_prob(i,:)= prob;
        BER_static_prob(i,:) = unshifed_pairwise_prob(i) / log2(mod01); %assumption: snr is sufficiently high

        
        %SER AND BER for various moving cases
        for j = 1:1:numel(subcarrier_frequencies)
            %get phase and doppler shifts and create modulator objects
            shift_rad = subcarrier_shifts(i,j);
            shift_deg = rad2deg(shift_rad);
            doppler_shift_delta = subcarrier_doppler(i,j) - subcarrier_frequencies(j);

            pfo = comm.PhaseFrequencyOffset(PhaseOffset=shift_deg, FrequencyOffset=0);
            pfo_freq = comm.PhaseFrequencyOffset(PhaseOffset=shift_deg, FrequencyOffset=doppler_shift_delta);
            
            %apply phase shift/freq shift to constellations
            shifted_mod = pfo(unshifted_mod);
            shifted_freq_mod = pfo_freq(unshifted_mod);
        
            %for each case, find odds that symbols in shifted constellation
            %would be misread as something else
            [similarity_mtx, closeidx] = find_smallest_distance_dual(unshifted_mod, shifted_mod);
            [prob_mismatch, prob_match] = prob_overreach_shift(similarity_mtx, noise_interference_ref);
            %numel(prob_mismatch)
            shifted_pairwise_error_prob(i,j) = prob_mismatch; %-> THIS PROBABILITY IS 1 - E SYMBOL ERROR RATE

            [similarity_mtx_f, closeidx_f] = find_smallest_distance_dual(unshifted_mod, shifted_freq_mod);
            [prob_mismatch_f, prob_match_f] = prob_overreach_shift(similarity_mtx_f, noise_interference_ref);
            shifted_pairwise_error_prob_freq(i,j) = prob_mismatch_f;  %-> THIS PROBABILITY IS  SYMBOL ERROR RATE
            
            %now compute BER
            BER_shift_prob(i,j) = shifted_pairwise_error_prob(i,j) / log2(mod01);  %BER phase shift
            BER_shift_freq_prob(i,j) = shifted_pairwise_error_prob_freq(i,j) / log2(mod01); %BER freq + phase shift


        end
    else %if disconnect case (0-ary mod), bit and symbol error rates are 100% 
        
        BER_static_prob(i, :) = 1; %BER static case
        BER_shift_prob(i,:) = 1;  %BER phase shift
        BER_shift_freq_prob(i,:) = 1; %BER freq + phase shift
        shifted_pairwise_error_prob_freq(i,:) = 1; %SER for freq + phase shift computation
        shifted_pairwise_error_prob(i,:) = 1; %SER for phase shift computation
        unshifed_pairwise_prob(i,:) = 1; %SER for static computation
    end
    
end    
save("lastrun.mat"); %saving machine state upon completion as this step is quite long and loaded ram data is easy to accidentally modify without meaning to
%you should be thanking me i am saving you hours of your time perhaps 
% compute average BER for static and moving cases at each waypoint
%Symbol error rates
%static_error_rate = mean(unshifed_pairwise_prob, 2);
moving_error_rate = mean(shifted_pairwise_error_prob, 2);
moving_error_rate_freq = mean(shifted_pairwise_error_prob_freq, 2);
%bit error rates = 
static_ber = BER_static_prob; %not channel dependent
freq_shift_ber = mean(BER_shift_freq_prob, 2);
shift_ber = mean(BER_shift_prob, 2);



%% compute datarates over trajectory here

%we have an SNIR map. S = SNIR * (N + I)

interference_power = ambient_interference; %dbm

NI = sum_dbm_power(interference_power, rxNoisePower); %dbm

signalPower_at_waypoints = SINR_at_waypoint + NI;

datarates_static = zeros(1, numel(interpolated_longs));
datarates_moving = zeros(1, numel(interpolated_longs));
datarates_moving_freq = zeros(1, numel(interpolated_longs));
CR = zeros(subcarrier_channel_count, numel(interpolated_longs)).';
bits_per_symbol = zeros(subcarrier_channel_count, numel(interpolated_longs)).';

%effective datarates after applying coding rate
datarates_eff_static = zeros(subcarrier_channel_count, numel(interpolated_longs));
datarates_eff_moving = zeros(subcarrier_channel_count, numel(interpolated_longs));
datarates_eff_moving_freq = zeros(subcarrier_channel_count, numel(interpolated_longs));

for k = 1:1:numel(interpolated_longs)
    mod_type = MOD_at_waypoints_BLER01(k);
    coding_rate = CR_at_waypoints_BLER01(k);
    if mod_type ~= 0
        bits_per_symbol(k,:) = log2(mod_type);
        CR(k,:) = coding_rate;
    elseif mod_type == 0
        bits_per_symbol(k,:) = 0;
        CR(k,:) = 0;
    end
end
%%
datarates_eff_static = bits_per_symbol.* symbol_rate.* (1 - BER_static_prob) .* CR;
datarates_eff_phase = bits_per_symbol.* symbol_rate.* (1 - BER_shift_prob).* CR;
datarates_eff_freq = bits_per_symbol.* symbol_rate.* (1 - BER_shift_freq_prob).* CR;


%%


%
%difference between static and moving datatrates
delta_static_phase = datarates_eff_static - datarates_eff_phase;

%total areas of datarate * distance (not time dependent -> higher = better
%signal) -> area-under-curve
area_datarate_static = subdivide_length * sum(datarates_eff_static, "all") %should always be highest
area_datarate_phase = subdivide_length * sum(datarates_eff_phase, "all")
area_delta = subdivide_length * sum(delta_static_phase, "all")
%area_datarate_phase_freq = subdivide_length * sum(datarates_eff_moving_freq) %dont use

%fraction of total area-under-curve vs reference (static) area-under-curve)
phase_frac = area_datarate_phase /area_datarate_static
delta_frac = area_delta/area_datarate_static


save(name_of_saved_filestate);


%%
numel(datarates_eff_static)

datarates_static_flat = reshape(datarates_eff_static.',1,[]);
datarates_phase_flat = reshape(datarates_eff_phase.',1,[]);

dr = [datarates_static_flat; datarates_phase_flat].';
% histogram(datarates_eff_static, 10)
% histogram(datarates_eff_phase, 10)
hist(dr, 10)
legend(["Non-shifted/static case", "Phase shifted/moving case"])
xlabel("Per-channel effective datarate (bits/s)")
ylabel("Frequency over trajectory")
titlestr = sprintf(['Distribution of per-channel datarates (fc = %.1f GHz, bw = %.1f MHz). ' ...
    '\n User velocity %.0f km/h, TX power %.1f dBm, ambiant interference. %.1f dBm \n ' ...
    'Ideal-to-moving data transfer ratio %.3f %% '],(fq/(10^9)), (rx_bw / (10^6)) ,(moving_speed * 3.6),tx_power,ambient_interference,(1-delta_frac)*100 )
title(titlestr)










