kernel_dir = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Kernel_fits/dnet/spechill/lambda_ind/batch_4';
save_dir = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Fits';
tau_type = 'sq';
model = 'dnet';

r_type{1} = 'big';
r_type{2} = 'small';
r_type{3} = 'anech';
n_rooms = length(r_type);
%% Criteria
NPSP_th = 40;

%% Get the names 
[temp,batch_name] = fileparts(kernel_dir);
[temp,lambda_name] = fileparts(temp);
[~,coch_name] = fileparts(temp);
%% Load the data
load(fullfile(kernel_dir,'info'),'info');

%% Select the data
ix = info.NPSP<NPSP_th;
NPSPs = info.NPSP(ix);
cluster_ids = info.cluster_id(ix);
animal_names = info.animal_name(ix);
pen_names = info.pen_name(ix);
% qualities = info.quality(ix);

%Sort in increasing NPSP
[NPSPs,ix_select] = sort(NPSPs,'ascend');
animal_names = animal_names(ix_select);
pen_names = pen_names(ix_select);
cluster_ids = cluster_ids(ix_select);
% qualities = qualities(ix_select);
n_clust = length(cluster_ids);
%% First load all the selected kernels
fprintf('== Extracting DNet features ==\n');tic;
for k = 1:n_clust
    c_name = fullfile(kernel_dir,strjoin({animal_names{k},pen_names{k},num2str(cluster_ids(k))},'_'));
    load(c_name,'kernel');
    dt_ms = round(kernel.dt_ms);
    nh = kernel.n_h;
    nf = length(kernel.freqs);
    for r = 1:n_rooms
        theta = kernel.(r_type{r}).theta;
        effective_HU{k}.(r_type{r}) = get_STRF_from_DNet(theta,dt_ms,nf,nh,tau_type);
    end
end
freqs = fliplr(kernel.freqs);
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

model_dir = fullfile(save_dir,model);
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