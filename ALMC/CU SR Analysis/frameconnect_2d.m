%[M,txt,raw] = xlsread('C:\Users\schwsama\Documents\Data\Isabel\TE20_AF647_3D_7_new.csv');
timewindow = 10;
ti%Mcopy = M;
M = [data.SML_data.ch2.frame,data.SML_data.ch2.position_x,data.SML_data.ch2.position_y,data.SML_data.ch2.precision];
c
for frame = 1:200%max(M(:,1)) %step through each frame
    currtimelocs = M(:,1) == frame;
    currtimelocIDs = find(currtimelocs);
    if isempty(currtimelocIDs)
        continue;
    end
    for locsATframe = 1:size(currtimelocIDs,1) %step through each localization at frame tt
        neartimelocs = (M(:,1) > frame) & (M(:,1) <= frame+timewindow);
        neartimelocIDs = find(neartimelocs);
        % note need to do above step inside the loops because
        % M is updated each step
        badids = [];
        for ntlocs = 1:size(neartimelocIDs,1) %step through each localization that is close in time
            if abs(M(currtimelocIDs(locsATframe),2) - M(neartimelocIDs(ntlocs),2)) < 3*M(currtimelocIDs(locsATframe),4) ||...
               abs(M(currtimelocIDs(locsATframe),3) - M(neartimelocIDs(ntlocs),3)) < 3*M(currtimelocIDs(locsATframe),4)
                % remove it if it's too close to 
                badids = [badids neartimelocIDs(ntlocs)];
            end
        end
        M(badids,:) = [];    
    end
end
toc

%% plotting

figure;scatter3(M(:,2)',M(:,3)',M(:,4)',4,M(:,2)'); title('Time Collapsed Fits');
xlabel('XDistance'); ylabel('YDistance'); zlabel('ZDistance');

figure;scatter3(Mcopy(:,5)',Mcopy(:,6)',Mcopy(:,7)',4,Mcopy(:,2)'); title('All Fits');
xlabel('XDistance'); ylabel('YDistance'); zlabel('ZDistance');

