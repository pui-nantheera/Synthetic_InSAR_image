clear all

rootDir = 'G:/VolcanicUnrest/Atmosphere/synthesised_patches/';

% parameters
samplesPerClass = 10000;

% run for each set
for setnum = 1:2
    % input directories
    patchDirWrap = [rootDir, 'set', num2str(setnum),'/wrap/'];
    patchDirUnwrap = [rootDir, 'set', num2str(setnum),'/unwrap/'];
    deformList    = dir([patchDirUnwrap,'deform/*.mat']);
    turbulentList = dir([patchDirUnwrap,'turbulent/*.mat']);
    stratifiedList = dir([patchDirUnwrap,'stratified/*.mat']);
    % get shuffled index for combination
    indDeform = randperm(length(deformList),samplesPerClass);
    indTurbulent = randperm(length(turbulentList),samplesPerClass);
    indStratified = randperm(length(stratifiedList),samplesPerClass);
    % output directories
    outputDirWrap = [patchDirWrap,'combine/'];
    outputDirUnwrap = [patchDirUnwrap, 'combine/'];
    mkdir(outputDirWrap);
    mkdir(outputDirUnwrap);
    % marging process
    for k = 1:samplesPerClass
        % get deformation
        load([patchDirUnwrap, 'deform/', deformList(indDeform(k)).name(1:end-3),'mat']);
        load([patchDirUnwrap, 'turbulent/', turbulentList(indTurbulent(k)).name(1:end-3),'mat']);
        load([patchDirUnwrap, 'stratified/', stratifiedList(indStratified(k)).name(1:end-3),'mat']);
        if range(los_grid(:))<=15
            los_grid = los_grid*18/range(los_grid(:));
            save([patchDirUnwrap, 'deform/', deformList(indDeform(k)).name(1:end-3),'mat'],'los_grid');
        elseif range(los_grid(:))>=50
            los_grid = los_grid*40/range(los_grid(:));
            save([patchDirUnwrap, 'deform/', deformList(indDeform(k)).name(1:end-3),'mat'],'los_grid');
        end
        insarImg = los_grid + curTur + atmo;
        mask = imerode(atmo~=0,strel('disk',3));
        insarWrap = (wrapTo2Pi(insarImg)-pi);
       
        
        insarWrap = (insarWrap-min(insarWrap(:)))/range(insarWrap(:)).*mask;
        imwrite(insarWrap, [outputDirWrap, outputName, '.png']);

    end
end

