 %% Match Penguin Foraging Dives to FTLE values
 % based off of 'resample_ftle.m'
 % 29April2024

load '/Users/jveatch/Documents/MATLAB/LCS-Tool-master/demo/ocean_dataset/ftle_LR_binary.mat'
load '/Users/jveatch/Documents/MATLAB/SWARM/Penguin/Penguin_withDist_nonInterp.mat'

load 'LCS_coverage_lat.mat';
load 'LCS_coverage_lon.mat';
%% add in extra days of FTLE coverage

load '/Volumes/T7_Shield/jmv208/LCS-Tool-master/demo/ocean_dataset/ftle_all_jan9jan14_LR.mat'

ftle_ext.ftle = cat(3, ftle_all.ftle, ftle_LR.ftle);
ftle_ext.time = [datenum(ftle_all.time), ftle_LR.time];
ftle_ext.domain = ftle_LR.domain;
ftle_ext.resolution = ftle_LR.resolution;
ftle_ext.x = ftle_LR.x;
ftle_ext.y = ftle_LR.y;


%% loop through ACRO observations, pair to ftle

start_codar = min(ftle_ext.time);
end_codar = max(ftle_ext.time);

polygon = [lon_shrink(:), lat_shrink(:)]; % bounds created with 'create_LCS_bound_random_resample.m'

% ECO_all.rpd_match = NaN(size(ECO_all.midtime));
Penguin.ftle_match = NaN(size(Penguin.datenum));
% create grid

ftle_2d = reshape(ftle_ext.ftle, [4400, 1394]);
[X,Y] = meshgrid(ftle_ext.x, ftle_ext.y);
X = reshape(X, [4400,1]);
Y = reshape(Y, [4400,1]);
ftle_coords = [X,Y];

counter = 1;
for i=start_codar:hours(1):end_codar
    
    start = datenum(i - minutes(30));
    finish = datenum(i + minutes(30));
    ind = find((Penguin.datenum>=datenum(start))&(Penguin.datenum<=datenum(finish)));
    if isempty(ind)
        
    else
        lat = Penguin.lat(ind);
        lon = Penguin.lon(ind);
        point = NaN([length(ind),2]);
        point(:,1) = lon;
        point(:,2) = lat;
%         k = dsearchn(RPD_coords,point);
        [cIdx,cD] = knnsearch(ftle_coords,point,'K',1,'Distance','chebychev');
        for j = 1:(length(ind))
            ind_prof = ind(j);
            lon = Penguin.lon(ind_prof);
            lat = Penguin.lat(ind_prof);
            
            if inpolygon(lon, lat, lon_shrink, lat_shrink)
            
                value = ftle_2d(cIdx(j,:), counter);

                Penguin.ftle_match(ind_prof) = nanmean(value);
                
            else
                Penguin.ftle_match(ind_prof) = NaN;
            end
%             gather_index
        end
        
%         value = rpd_2d(ind,counter);
%         ECO_all.rpd_match(ind) = value;
    end
    counter = counter +1;
end

%% seperate species

indAD = find(strcmp(Penguin.spp, 'ADPE'));
indGE = find(strcmp(Penguin.spp, 'GEPE'));

AD.diveDepth = Penguin.diveDepth(indAD);
AD.ftle_match = Penguin.ftle_match(indAD);
AD.lon = Penguin.lon(indAD);
AD.lat = Penguin.lat(indAD);
AD.time = Penguin.datetime(indAD);
AD.colony = Penguin.colony(indAD);
AD.TagRef = Penguin.TagID(indAD);

GE.diveDepth = Penguin.diveDepth(indGE);
GE.ftle_match = Penguin.ftle_match(indGE);
GE.lon = Penguin.lon(indGE);
GE.lat = Penguin.lat(indGE);
GE.time = Penguin.datetime(indGE);
GE.colony = Penguin.colony(indGE);
GE.TagRef = Penguin.TagID(indGE);




%% repeat for null model (random walk penguin trajectories)

dataTable = readtable('/Users/jveatch/Documents/MATLAB/SWARM/Penguin/sim_ADPE.csv'); % from Matt
% create structure
AD_random.ID = dataTable.sim_penID;
AD_random.hour = dataTable.ADPE_sim_pen_hr;
AD_random.lat = dataTable.ADPE_lat_sim_out;
AD_random.lon = dataTable.ADPE_lon_sim_out;
AD_random.hour_datenum = AD_random.hour./24;

% assign day to each simulated penguin
pen_ID = unique(AD_random.ID);

start_codar = datenum(datetime(2020, 1, 09, 00, 00, 0));
end_codar = datenum(datetime(2020, 1, 19, 00, 00, 0));

sim_days = [start_codar:1:end_codar]; 

AD_random.datenum = NaN(length(AD_random.hour),1);
c = 1;

for i = 1:length(sim_days)
    
    for j = 1:10
        ind = find(strcmp(AD_random.ID, pen_ID(c)));  % for each simulated penguin, assign a day
        day = sim_days(i);
        AD_random.datenum(ind) = AD_random.hour_datenum(ind) + day;
        c = c+1;
    end
    
end

ind = ~isnan(AD_random.datenum);
AD_random.ID = AD_random.ID(ind);
AD_random.hour = AD_random.hour(ind);
AD_random.lat = AD_random.lat(ind);
AD_random.lon = AD_random.lon(ind);
AD_random.hour_datenum = AD_random.hour_datenum(ind);
AD_random.datenum = AD_random.datenum(ind);
%% match FTLE to simulated penguin in space and time --> Adelie

% ECO_all.rpd_match = NaN(size(ECO_all.midtime));
AD_random.ftle_match = NaN(size(AD_random.datenum));
% create grid

counter = 1;
for i=start_codar:hours(1):end_codar
    
    start = datenum(i - minutes(30));
    finish = datenum(i + minutes(30));
    ind = find((AD_random.datenum>=datenum(start))&(AD_random.datenum<=datenum(finish)));
    if isempty(ind)
        
    else
        lat = AD_random.lat(ind);
        lon = AD_random.lon(ind);
        point = NaN([length(ind),2]);
        point(:,1) = lon;
        point(:,2) = lat;
%         k = dsearchn(RPD_coords,point);
        [cIdx,cD] = knnsearch(ftle_coords,point,'K',1,'Distance','chebychev');
        for j = 1:(length(ind))
            ind_prof = ind(j);
            
            lon = AD_random.lon(ind_prof);
            lat = AD_random.lat(ind_prof);
            
            if inpolygon(lon, lat, lon_shrink, lat_shrink)
            
                value = ftle_2d(cIdx(j,:), counter);

                AD_random.ftle_match(ind_prof) = nanmean(value);

                
            else
                 AD_random.ftle_match(ind_prof) = NaN;
            end
%             gather_index
        end
        
%         value = rpd_2d(ind,counter);
%         ECO_all.rpd_match(ind) = value;
    end
    counter = counter +1;
end
%% GENTOO simulated penguins

dataTable = readtable('/Users/jveatch/Documents/MATLAB/SWARM/Penguin/sim_GEPE.csv'); % from Matt
% create structure
GE_random.ID = dataTable.sim_penID;
GE_random.hour = dataTable.GEPE_sim_pen_hr;
GE_random.lat = dataTable.GEPE_lat_sim_out;
GE_random.lon = dataTable.GEPE_lon_sim_out;
GE_random.hour_datenum = GE_random.hour./24;

% assign day to each simulated penguin
pen_ID = unique(GE_random.ID);

start_codar = min(datenum(datetime(2020, 1, 11, 00, 00, 0)));
end_codar = max(datenum(datetime(2020, 1, 25, 00, 00, 0)));

sim_days = [start_codar:1:end_codar]; 

GE_random.datenum = NaN(length(GE_random.hour),1);
c = 1;

for i = 1:length(sim_days)
    
    for j = 1:8
        ind = find(strcmp(GE_random.ID, pen_ID(c)));  % for each simulated penguin, assign a day
        day = sim_days(i);
        GE_random.datenum(ind) = GE_random.hour_datenum(ind) + day;
        c = c+1;
    end
end

ind = ~isnan(GE_random.datenum);
GE_random.ID = GE_random.ID(ind);
GE_random.hour = GE_random.hour(ind);
GE_random.lat = GE_random.lat(ind);
GE_random.lon = GE_random.lon(ind);
GE_random.hour_datenum = GE_random.hour_datenum(ind);
GE_random.datenum = GE_random.datenum(ind);

%% match FTLE to simulated penguin in space and time --> Gentoo

% ECO_all.rpd_match = NaN(size(ECO_all.midtime));
GE_random.ftle_match = NaN(size(GE_random.datenum));
% create grid


counter = 1;
for i=start_codar:hours(1):end_codar
    
    start = datenum(i - minutes(30));
    finish = datenum(i + minutes(30));
    ind = find((GE_random.datenum>=datenum(start))&(GE_random.datenum<=datenum(finish)));
    if isempty(ind)
        
    else
        lat = GE_random.lat(ind);
        lon = GE_random.lon(ind);
        point = NaN([length(ind),2]);
        point(:,1) = lon;
        point(:,2) = lat;
%         k = dsearchn(RPD_coords,point);
        [cIdx,cD] = knnsearch(ftle_coords,point,'K',1,'Distance','chebychev');
        for j = 1:(length(ind))
            ind_prof = ind(j);
            
            lon = GE_random.lon(ind_prof);
            lat = GE_random.lat(ind_prof);
            
            if inpolygon(lon, lat, lon_shrink, lat_shrink)
            
                value = ftle_2d(cIdx(j,:), counter);

                GE_random.ftle_match(ind_prof) = nanmean(value);
            else
                 GE_random.ftle_match(ind_prof) = NaN;
            end
%             gather_index
        end
        
%         value = rpd_2d(ind,counter);
%         ECO_all.rpd_match(ind) = value;
    end
    counter = counter +1;
end




%% plot some histograms


figure(1)
histogram(AD.ftle_match(AD.diveDepth > 5));
hold on;
title('AD > 5');

figure(2)
histogram(AD.ftle_match(AD.diveDepth <= 5));
hold on;
title('AD <= 5');

figure(3)
histogram(GE.ftle_match(GE.diveDepth > 5));
hold on;
title('GE > 5');

figure(4)
histogram(GE.ftle_match(GE.diveDepth <= 5));
hold on;
title('GE <= 5');

ind_lowFTLE_AD = find(AD.ftle_match < 0.05);
ind_lowFTLE_GE = find(GE.ftle_match < 0.05);

ind_highFTLE_AD = find(AD.ftle_match >= 0.05);
ind_highFTLE_GE = find(GE.ftle_match >= 0.05);

figure(5)
scatter(AD.lon(ind_lowFTLE_AD), AD.lat(ind_lowFTLE_AD));
scatter(GE.lon(ind_lowFTLE_GE), GE.lat(ind_lowFTLE_GE));

AD_random.dt = datetime(AD_random.datenum, 'ConvertFrom', 'datenum');
GE_random.dt = datetime(GE_random.datenum, 'ConvertFrom', 'datenum');

dataTable = struct2table(AD);
writetable(dataTable, 'AD_ftle.csv');

dataTable = struct2table(GE);
writetable(dataTable, 'GE_ftle.csv');

dataTable = struct2table(AD_random);
writetable(dataTable, 'AD_random_ftle.csv');

dataTable = struct2table(GE_random);
writetable(dataTable, 'GE_random_ftle.csv');

