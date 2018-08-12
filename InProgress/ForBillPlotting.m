%% 
% b = [];
% now open variable 'b' in 'Workspace Tab';
%%
clear xs;
base = ones(1,size(b,1));
for ii = 1:size(b,2)     
xs(ii,:) = base + ii-1;
end

figure; hold on;
% font = 'Times';
for ii = 1:size(b,2)    
scatter(b(:,ii),xs(ii,:),'marker','none');
text(b(:,ii),xs(ii,:),'|');
end
ylim([0 size(b,2)+3]);

