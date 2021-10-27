kernel_dir = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/DNet_fits';
save_dir = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Fits';
tau_type = 'sq';
model_name = 'dnet';
nh = 6;
nf = 34;
dt_ms = 5;
f_min = 200;
f_max = 22000;
NPSP_th = 40;
lambda_name = 'lambda_ind';
batch_name = 'batch_36';
coch_name = 'spechill';
r_type{1} = 'big';
r_type{2} = 'small';
r_type{3} = 'anech';
n_rooms = length(r_type);
%% Criteria

%% Load the data
load(fullfile(kernel_dir,'noise_ratios'),'NR','neuronlist');
load(fullfile(kernel_dir,'freqs'),'freqs'); % Get the freqs that were used
%% Select the data
%Sort in increasing NPSP
[NPSPs,ix_select] = sort(NR,'ascend');
n_clust = length(ix_select);

for k = 1:n_clust
    ix = ix_select(k);
    temp = split(neuronlist{ix},'_');
    animal_names = temp{1};
    pen_names = temp{2};
    cluster_ids = temp{3}(1:end-4);
end

%% First load all the selected kernels
fprintf('== Extracting DNet features ==\n');tic;
for r = 1:n_rooms
    c_name = fullfile(kernel_dir,['run_dn_mlp_Reverb_Datacond_',r_type{r},'_30ms']);
    load(c_name,'model');
    for k = 1:n_clust
        ix = ix_select(k);
        theta = model(ix).theta;
        effective_HU{k}.(r_type{r}) = get_STRF_from_DNet(theta,dt_ms,nf,nh,tau_type);
    end
end
fprintf('== Done! This took %0.fs ==\n',toc);
%% Get things out which will be useful for plotting
for r = 1:n_rooms
    room = r_type{r};
    taus.(room).all = [];
    ie_score.(room).all = [];
    bf.(room).all = [];
end

for k = 1:n_clust
    for r = 1:n_rooms
        room = r_type{r};
        n_hu = length(effective_HU{k}.(room));
        temp_tau = [effective_HU{k}.(room).STRF_tau]';
        temp_ie = [effective_HU{k}.(room).STRF_IE_score]';
        taus.(room).mean(k) = mean(temp_tau);
        taus.(room).std(k) = std(temp_tau);
        taus.(room).all = [taus.(room).all;temp_tau];
        ie_score.(room).all = [ie_score.(room).all;temp_ie];
        for hu = 1:n_hu
            k_fh = effective_HU{k}.(room)(hu).STRF_weights;
            [~,ix_bf] = max(abs(k_fh(:)));
            [ix_bf,~] = ind2sub(size(k_fh),ix_bf);
            bf_hu = freqs(ix_bf);
            bf.(room).all = [bf.(room).all;bf_hu];
        end
    end
end

%% Plot a histogram with all the data 
row = 2;
col = 2;
per = 0.01;
edgel = 0.05; edger = per; edgeh = 0.03; edgeb = 0.06; space_h = 0.03; space_v = 0.1;
[pos]=subplot_pos(row,col,edgel,edger,edgeh,edgeb,space_h,space_v);

c{1} = 'r';
c{2} = 'g';
c{3} = 'b';
axis_sz = 15;
font_sz = axis_sz;

model_dir = fullfile(save_dir,model_name);
if ~exist(model_dir,'dir')
    mkdir(model_dir);
end

coch_dir = fullfile(model_dir,coch_name);
if ~exist(coch_dir, 'dir')
    mkdir(coch_dir);
end

lambda_dir = fullfile(coch_dir,lambda_name);
if ~exist(lambda_dir, 'dir')
    mkdir(lambda_dir);
end

save_dir_full = fullfile(lambda_dir,batch_name);
if ~exist(save_dir_full, 'dir')
    mkdir(save_dir_full);
end
            
figure('units','normalized','outerposition',[0 0 1 1]);
edges = logspace(log10(5),log10(500),20);
for r = 1:n_rooms
    room = r_type{r};
    subplot('position',pos{1});
    hold on;
    histogram(taus.(room).all,edges,'FaceColor','none','EdgeColor',c{r},'DisplayStyle','stairs','LineWidth',2,'Normalization','probability');
    
    subplot('position',pos{3});hold on;       
    plot(taus.(room).all,ie_score.(room).all,'.','Color',c{r});
    
    subplot('position',pos{4});hold on;
    plot(bf.(room).all,taus.(room).all,'.','Color',c{r});
    
end

subplot('position',pos{2});
x = [1 2 3];
data = [mean(taus.big.mean) mean(taus.small.mean) mean(taus.anech.mean)]';
errhigh = [std(taus.big.mean) std(taus.small.mean) std(taus.anech.mean)];
errlow  = [0 0 0];
bar(x,data);
set(gca,'xticklabel',{'Big','Small','Anech'})
hold on;
er = errorbar(x,data,errlow,errhigh);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

set(gcf,'color','w');

subplot('position',pos{1});
set(gca, 'xscale','log');
xlabel('\tau [ms]','FontSize',font_sz,'FontWeight','bold');
ylabel('Probability','FontSize',font_sz,'FontWeight','bold');
title('Histograms of \tau');
legend('Big','Small','Anech');
set(gca,'FontName','Arial','FontSize',axis_sz,'FontWeight','Bold');

subplot('position',pos{2});
xlabel('Room','FontSize',font_sz,'FontWeight','bold');
ylabel('\tau [ms]','FontSize',font_sz,'FontWeight','bold');
set(gca,'FontName','Arial','FontSize',axis_sz,'FontWeight','Bold');

subplot('position',pos{3});
xlabel('\tau [ms]','FontSize',font_sz,'FontWeight','bold');
ylabel('IE score','FontSize',font_sz,'FontWeight','bold');
title('IE score vs \tau');
legend('Big','Small','Anech');
set(gca,'FontName','Arial','FontSize',axis_sz,'FontWeight','Bold');
xlim([0 500]);

subplot('position',pos{4});
xlabel('Freqeuncy [kHz]','FontSize',font_sz,'FontWeight','bold');
ylabel('\tau [ms]','FontSize',font_sz,'FontWeight','bold');
title('\tau vs Frequency');
legend('Big','Small','Anech');
set(gca, 'XScale', 'log');
set(gca,'FontName','Arial','FontSize',axis_sz,'FontWeight','Bold');
ylim([0 500]);

save_name = fullfile(save_dir_full,['DNet plots of tau for all neurons NPSP<',num2str(NPSP_th),' cochlea ',coch_name,' ',lambda_name,' ',batch_name,'.png']);
export_fig(save_name);
close all;
%% Plot histograms as subplots to see better