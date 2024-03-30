
function [] = generate_channelKey(saveDir, nCh)
names = [ "Nucleus"; repmat("ch name", nCh-1, 1)];

isLocalizable = vertcat([0], repmat([1], nCh-1, 1));

uM_FISH_ch = (1:nCh)';

uM_Bead_ch = (0:nCh-1)';


scope = repmat("Nikon", nCh, 1);

dye = ["DAPI"; repmat("fluor", nCh-1,1)];


ckT = table(names, isLocalizable, uM_FISH_ch, uM_Bead_ch, scope, dye)

save_path = fullfile(saveDir, 'channelKey.txt'); % save destianion

writetable(ckT, save_path);

disp('Channel Key Table saved to:')
disp(save_path);
disp('Please open and modify accordingly')

end