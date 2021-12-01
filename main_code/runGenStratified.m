% This code is for generate 'wrapped' stratified atmospheric delay

clear all

SAVEWRAP = 0;
inputRoot = '/Volumes/VolcanicUnrest/Atmosphere/stratify/volcano_2018/';
outputRoot = '/Volumes/VolcanicUnrest/Atmosphere/synthesised_patches/';

if SAVEWRAP == 1
    mkdir([outputRoot, 'set1\wrap\stratified\']);
    mkdir([outputRoot, 'set2\wrap\stratified\']);
end

% areas of volcano
volcanoList = dir(inputRoot);
volcanoList(1:2) = [];

% parameters
totalSamples = 10000; % total sample to generate
imageSize = 500;      % resolution in pixels
incidence=43.7835;
xref  = 115.1;
yref  = -8.4;
dxref = 10;
dyref = dxref;
wavelength = 0.055465;
m2rad = 4.*pi./wavelength;
rad2m = (wavelength./(4.*pi));
zen2los = 1./cos(incidence./180.*pi);
halfcrop = floor(227/2); % as input of Alexnet is 227x227 pixels

gaptime = 16; % days
count = 0;
addind = 1;
for k = 1:length(volcanoList)
    disp([num2str(k),'/',num2str(length(volcanoList))])
    % point to target output directorty
    if SAVEWRAP == 1
        outputDirWrap = [outputRoot, 'set', num2str(2-rem(k,2)),'/wrap/stratified/'];
    end
    outputDirUnwrap = [outputRoot, 'set', num2str(2-rem(k,2)),'/unwrap/stratified/'];
    mkdir(outputDirUnwrap);
    inputDir = [inputRoot, volcanoList(k).name,'/'];
    % get dates
    dateNames = dir([inputDir, '*.ztd']);
    % get stratified atmosphere with 12 days interval
    for n = 1:length(dateNames)-gaptime
        master = dateNames(n).name(1:end-4);
        slave  = dateNames(n+gaptime).name(1:end-4);
        
        [~,~,atmo1] = read_GACOS([inputDir,master]);
        [~,~,atmo2] = read_GACOS([inputDir,slave]);
        
        atmo = (atmo2-atmo1).*zen2los.*m2rad;
        atmo = imresize(atmo, [imageSize imageSize]);
        mask = imerode(atmo~=0,strel('disk',3));
        atmo = atmo(round(size(atmo,1)/2) + (-halfcrop:halfcrop),round(size(atmo,2)/2) + (-halfcrop:halfcrop));
        mask = mask(round(size(atmo,1)/2) + (-halfcrop:halfcrop),round(size(atmo,2)/2) + (-halfcrop:halfcrop));
        mask = mask(end:-1:1,:);
        if (~(((sum(mask(:)==0)/numel(mask)) <= 0.1) || ((sum(mask(:)==0)/numel(mask)) >= 0.5))) && ... % && (range(atmo(:))>pi/2)
                mask(halfcrop,halfcrop)>0
            
            save([outputDirUnwrap, 'S', sprintf('%05d',count + addind), '.mat'], 'atmo');
            if SAVEWRAP == 1
                atmo = wrapTo2Pi(atmo)-pi;
                atmo = (atmo-min(atmo(:)))/range(atmo(:)).*mask;
                imwrite(atmo, [outputDirWrap, 'S', sprintf('%05d',count + addind), '.png']);
            end
            count = count + 1;
        else
            break;
        end
        if count >= totalSamples
            break;
        end
    end
end

%% check
% 
% checkDir = 'G:\VolcanicUnrest\Atmosphere\synthesised_patches\set2\wrap\stratified\';
% namefiles = dir([checkDir, '*.png']);
% 
% for k = 1:length(namefiles)
%     img = im2double(imread([checkDir, namefiles(k).name]));
%     if (sum(img(:)==0)/numel(img)) > 0.75
%         delete([checkDir, namefiles(k).name]);
%     end
% end
