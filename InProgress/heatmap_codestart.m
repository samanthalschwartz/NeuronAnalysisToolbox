test{1} = floor(rand(5,1)*1000);
test{2} = floor(rand(5,1)*1000);
test{3} = floor(rand(5,1)*1000);
times = 1:max([max(test{1}),max(test{2}),max(test{3})]);
inzeros = times*0;
vals = cell(numel(test),1);
for nn = 1:numel(vals)
   vals{nn} = inzeros;
   vals{nn}(test{nn}) = 1;
end

figure;
for ii = 1:numel(vals)
    for jj = 1:numel(vals{ii})
        if vals{ii}(jj)
            FaceColor = 'k';
        else
            FaceColor = 'w';
        end
        p = patch([times(jj),times(jj)+1,times(jj)+1,times(jj)],[ii-1, ii-1, ii, ii],vals{ii}(jj));
   set(p,'FaceColor',FaceColor,'EdgeColor','none');
    end
end