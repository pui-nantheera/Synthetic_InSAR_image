% This code is for generate turbu atmospheric delay

clear all

addpath('..\turbulent\');

outputRoot = 'G:\VolcanicUnrest\Atmosphere\synthesised_patches\';
SAVEWRAP = 0;
mkdir([outputRoot, 'set1\unwrap\turbulent\']);
mkdir([outputRoot, 'set2\unwrap\turbulent\']);
if SAVEWRAP == 1
    mkdir([outputRoot, 'set1\wrap\turbulent\']);
    mkdir([outputRoot, 'set2\wrap\turbulent\']);
end
% parameters
rows = 100;
cols = 100;
psizex = 1;
psizey = 1;
covmodel_type = 0;
N = 20000/25;
halfcrop = floor(227/2); % size of input of Alexnet = 227x227 pixels
imageSize = 500;         % resolution in pixels

for maxvar = 7.5+[-2 -1 0 0.75 1.5]
    for alpha = 0.008*[0.5 0.75 1 1.5 2]
        disp([sprintf('%0.2f',maxvar),'_', sprintf('%0.4f',alpha)]);
        % generate turbulent atmosphere
        atm_pets = pcmc_atm(rows,cols,maxvar,alpha,covmodel_type,N,psizex,psizey);
        for k = 1:N
            curTur = imresize(atm_pets(:,:,k),[imageSize imageSize]);
            curTur = curTur(round(size(curTur,1)/2) + (-halfcrop:halfcrop),round(size(curTur,2)/2) + (-halfcrop:halfcrop));
            
            outputDir = [outputRoot, 'set', num2str(2-rem(k,2)),'\unwrap\turbulent\'];
            save([outputDir, 'turb_', sprintf('%0.2f',maxvar),'_', sprintf('%0.4f',alpha), '_', sprintf('%03d', k), '.mat'],'curTur');
            if SAVEWRAP == 1
                curTur = wrapTo2Pi(curTur)-pi;
                curTur = (curTur-min(curTur(:)))/range(curTur(:));
                outputDir = [outputRoot, 'set', num2str(2-rem(k,2)),'\wrap\turbulent\'];
                imwrite(curTur, [outputDir, 'turb_', sprintf('%0.2f',maxvar),'_', sprintf('%0.4f',alpha), '_', sprintf('%03d', k), '.png']);
            end
        end
    end
end

