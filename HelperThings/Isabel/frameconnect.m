[M,txt,raw] = xlsread('C:\Users\schwsama\Documents\Data\Isabel\TE20_AF647_3D_7_new.csv');
timewindow = 40;
Mcopy = M;
for frame = 1:max(M(:,2)) %step through each frame
    currtimelocs = M(:,2) == frame;
    currtimelocIDs = find(currtimelocs);
    if isempty(currtimelocIDs)
        continue;
    end
    for locsATframe = 1:size(currtimelocIDs,1) %step through each localization at frame tt
        neartimelocs = (M(:,2) > frame) & (M(:,2) <= frame+timewindow);
        neartimelocIDs = find(neartimelocs);
        % note need to do above step inside the loops because
        % M is updated each step
        badids = [];
        for ntlocs = 1:size(neartimelocIDs,1) %step through each localization that is close in time
            if abs(M(currtimelocIDs(locsATframe),5) - M(neartimelocIDs(ntlocs),5)) < 3*M(currtimelocIDs(locsATframe),8) ||...
               abs(M(currtimelocIDs(locsATframe),6) - M(neartimelocIDs(ntlocs),6)) < 3*M(currtimelocIDs(locsATframe),8) ||...
               abs(M(currtimelocIDs(locsATframe),7) - M(neartimelocIDs(ntlocs),7)) < 3*M(currtimelocIDs(locsATframe),9)
                % remove it if it's too close to 
                badids = [badids neartimelocIDs(ntlocs)];
            end
        end
        M(badids,:) = [];    
    end
end


%% plotting

figure;scatter3(M(:,5)',M(:,6)',M(:,7)',4,M(:,2)'); title('Time Collapsed Fits');
xlabel('XDistance'); ylabel('YDistance'); zlabel('ZDistance');

figure;scatter3(Mcopy(:,5)',Mcopy(:,6)',Mcopy(:,7)',4,Mcopy(:,2)'); title('All Fits');
xlabel('XDistance'); ylabel('YDistance'); zlabel('ZDistance');

