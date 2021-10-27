%% Hill function parameters 
SAT = 0.005; %This was found from crossvalidation
c = 0.001; %The normalization coeff c
n = 2; %The power to which x is raised
th = -75;
make_coch = 0;
% SAT = 0.16; %This was found from crossvalidation
% c = 0.01; %The normalization coeff c
% n = 1.77; %The power to which x is raised
%% Input
dx = 0.1;
x = [0:dx:130];
%% Make and load data
if make_coch
    save_distributions(th);
end
load('/mnt/40086D4C086D41D0/Reverb_normative/Cochleagrams/Troubleshoot/coch.mat');
edges = x;
for j = 1:4
counts.(coch(j).reverb_cond) = histcounts(coch(j).X_ft(:),edges);
end
%% Run and plot
row = 2;
col = 2;
per = 0.03;
edgel = per; edger = 0.02; edgeh = per; edgeb = per; space_h = per; space_v = 0.075;
[pos]=subplot_pos(row,col,edgel,edger,edgeh,edgeb,space_h,space_v);

lw = 3;
x_to_n = (c.*x).^n;
y = (x_to_n)./(SAT + x_to_n);
max_y = max(y);


for j = 1:4
    subplot('position',pos{j});
    hold on;
    plot(x,y,'Color','r','LineWidth',lw);
    xlabel('Input');
    ylabel('Output');
    title(['Hill function ',coch(j).reverb_cond]);
    histogram('BinEdges', edges,'BinCounts',(counts.(coch(j).reverb_cond)/max(counts.(coch(j).reverb_cond)))* max_y,'FaceColor','none','EdgeColor','b','DisplayStyle','stairs','LineWidth', lw);
    hold off;
end

