function ASD_responder_subtype
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
% Written by 
% Qi liu 
% Weihua Zhao zarazhao@uestc.edu.cn
% Last edited July 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear

load('data_ori.mat')
data=data_ori(:,1:23);  % features
response=data_ori(:,24); % reaponder label
[nsub,nfea]=size(data);
data_norm = zeros(nsub,nfea);

% normalization
for j=1:nfea
    norm_Min = min(data(:,j));
    norm_Max = max(data(:,j));
    data_norm(:,j) = (data(:,j)-norm_Min)./(norm_Max-norm_Min);
    
end



%% Cluster analysis with all features
data_index = [1:23];
data_new = data_norm(:,data_index);
maxk=10;
Kmeans_results ={};
% k-means cluster (Matlab2015b)
for k = 2 : maxk 
    disp(['Calculating for ' num2str(k) 'clusters'])
    [IDX, C, SUMD, D]=kmeans(data_new,k,'Distance','cityblock','Start','plus','MaxIter',10000,'Replicates',20,'Display','final'); %,'Options',opt);   
    Kmeans_results{k}.IDX=IDX;   
    Kmeans_results{k}.C=C;        
    Kmeans_results{k}.SUMD=SUMD; 
    Kmeans_results{k}.D=D; 
end

save Kmeans_results_allfea.mat Kmeans_results

load('Kmeans_results_allfea.mat')

% calculated variance ratio score (VRS)
Tol_CSS=zeros(10,1);      
SSB = zeros(10,1);
VRS= zeros(10,1);
for k = 2 : 10 
    
     for clu_num = 1:k
         for win_num = 1:length(data_new)
             if Kmeans_results{1,k}.IDX(win_num) == clu_num
                single_CSS = (data_new(win_num,:) - Kmeans_results{1,k}.C(clu_num,:))*(data_new(win_num,:) - Kmeans_results{1,k}.C(clu_num,:))';
                Tol_CSS(k)=Tol_CSS(k)+single_CSS;
             end
         end
         clu_index = find(Kmeans_results{1,k}.IDX == clu_num);
         SSB(k) = SSB(k)+length(clu_index)*(mean(data_new) - Kmeans_results{1,k}.C(clu_num,:))*(mean(data_new) - Kmeans_results{1,k}.C(clu_num,:))';
     end
     VRS(k) = SSB(k)/Tol_CSS(k)*(length(data_new)-k)/(k-1);
end
x=[2:1:10]';
plot(x,VRS(2:10),'o')
hold on
p= plot(2:10,VRS(2:10),'LineWidth',3);

% subtype_result
ori_res_sub1=mean(response(find(Kmeans_results{1,2}.IDX==1)));
ori_res_sub2=mean(response(find(Kmeans_results{1,2}.IDX==2)));
subtype_result_allfea(1,1)=length(find(Kmeans_results{1,2}.IDX==1)); %number of subtype 1
subtype_result_allfea(1,3)=length(find(Kmeans_results{1,2}.IDX==2)); %number of subtype 1
subtype_result_allfea(1,2)=ori_res_sub1; %response rate of subtype 1
subtype_result_allfea(1,4)=ori_res_sub2; %response rate of subtype 2

% permutation test for the response rate
per_res_sub1 = zeros(10000,1);
per_res_sub2 = zeros(10000,1);
nper = 10000;
for per_num = 1 : nper
    per_index = randperm(41);
    per_res = response(per_index);
    per_res_sub1(per_num,1) = mean(per_res(find(Kmeans_results{1,2}.IDX==1)));
    per_res_sub2(per_num,1) = mean(per_res(find(Kmeans_results{1,2}.IDX==2)));
end
subtype_result_allfea(1,6)=length(find(per_res_sub1>ori_res_sub1));  %permutation test result of subtype 1
subtype_result_allfea(1,7)=length(find(per_res_sub2>ori_res_sub2));  %permutation test result of subtype 2

% group differences between two subtypes on each feature
mean_value = zeros(23,6);
for i = 1 :23
    [h(i),p(i)]=ttest2(data(find(Kmeans_results{1,2}.IDX==1),i),data(find(Kmeans_results{1,2}.IDX==2),i));
    mean_value(i,1) = mean(data(find(Kmeans_results{1,2}.IDX==1),i));
    mean_value(i,2) = std(data(find(Kmeans_results{1,2}.IDX==1),i));
    mean_value(i,3) = length(find(Kmeans_results{1,2}.IDX==1));
    mean_value(i,4) = mean(data(find(Kmeans_results{1,2}.IDX==2),i));
    mean_value(i,5) = std(data(find(Kmeans_results{1,2}.IDX==2),i));
    mean_value(i,6) = length(find(Kmeans_results{1,2}.IDX==2));
end
P_fdr = mafdr(p,'BHFDR', true);  % FDR corrected



%% ablation analysis
for i = 1:23
    
    data_index = [1:i-1 i+1:23];
    data_new = data_norm(:,data_index);
    maxk=10;
    
    for k = 2 : maxk
        disp(['Calculating for ' num2str(k) 'clusters'])
        [IDX, C, SUMD, D]=kmeans(data_new,k,'Distance','cityblock','Start','plus','MaxIter',10000,'Replicates',20,'Display','final'); %,'Options',opt);
        Kmeans_results{k}.IDX=IDX;    
        Kmeans_results{k}.C=C;        
        Kmeans_results{k}.SUMD=SUMD;
        Kmeans_results{k}.D=D;
    end
    
    % calculated variance ratio score (VRS)
    Tol_CSS=zeros(10,1);      
    SSB = zeros(10,1);
    VRS= zeros(10,1);
    for k = 2 : 10
          
        for clu_num = 1:k
            for win_num = 1:length(data_new)
                if Kmeans_results{1,k}.IDX(win_num) == clu_num
                    single_CSS = (data_new(win_num,:) - Kmeans_results{1,k}.C(clu_num,:))*(data_new(win_num,:) - Kmeans_results{1,k}.C(clu_num,:))';
                    Tol_CSS(k)=Tol_CSS(k)+single_CSS;
                end
            end
            clu_index = find(Kmeans_results{1,k}.IDX == clu_num);
            SSB(k) = SSB(k)+length(clu_index)*(mean(data_new) - Kmeans_results{1,k}.C(clu_num,:))*(mean(data_new) - Kmeans_results{1,k}.C(clu_num,:))';
        end
        VRS(k) = SSB(k)/Tol_CSS(k)*(length(data_new)-k)/(k-1);
    end

    ori_res_sub1=mean(response(find(Kmeans_results{1,2}.IDX==1)));
    ori_res_sub2=mean(response(find(Kmeans_results{1,2}.IDX==2)));
    subtype_result_abla(i,1)=length(find(Kmeans_results{1,2}.IDX==1));
    subtype_result_abla(i,3)=length(find(Kmeans_results{1,2}.IDX==2));
    subtype_result_abla(i,2)=ori_res_sub1;
    subtype_result_abla(i,4)=ori_res_sub2;
    subtype_result_abla(i,5)=0.615-max(ori_res_sub1,ori_res_sub2); %difference from the response rate of above
    
    per_res_sub1 = zeros(10000,1);
    per_res_sub2 = zeros(10000,1);
    nper = 10000;
    for per_num = 1 : nper
        per_index = randperm(41);
        per_res = response(per_index);
        per_res_sub1(per_num,1) = mean(per_res(find(Kmeans_results{1,2}.IDX==1)));
        per_res_sub2(per_num,1) = mean(per_res(find(Kmeans_results{1,2}.IDX==2)));
    end
    subtype_result_abla(i,6)=length(find(per_res_sub1>ori_res_sub1));
    subtype_result_abla(i,7)=length(find(per_res_sub2>ori_res_sub2));
    
end


%% 7-index cluster analysis
data_index = [4 5 12 13 14 18 19]; 
data_new = data_norm(:,data_index);
maxk=10;
Kmeans_results_7index ={};
% k-means cluster (Matlab2015b)
for k = 2 : maxk 
    disp(['Calculating for ' num2str(k) 'clusters'])
    [IDX, C, SUMD, D]=kmeans(data_new,k,'Distance','cityblock','Start','plus','MaxIter',10000,'Replicates',20,'Display','final'); %,'Options',opt);   
    Kmeans_results_7index{k}.IDX=IDX;    
    Kmeans_results_7index{k}.C=C;      
    Kmeans_results_7index{k}.SUMD=SUMD; 
    Kmeans_results_7index{k}.D=D; 
end


% calculated variance ratio score (VRS)
Tol_CSS=zeros(10,1);      
SSB = zeros(10,1);
VRS= zeros(10,1);
for k = 2 : 10
    
     
     for clu_num = 1:k
         for win_num = 1:length(data_new)
             if Kmeans_results_7index{1,k}.IDX(win_num) == clu_num
                single_CSS = (data_new(win_num,:) - Kmeans_results_7index{1,k}.C(clu_num,:))*(data_new(win_num,:) - Kmeans_results_7index{1,k}.C(clu_num,:))';
                Tol_CSS(k)=Tol_CSS(k)+single_CSS;
             end
         end
         clu_index = find(Kmeans_results_7index{1,k}.IDX == clu_num);
         SSB(k) = SSB(k)+length(clu_index)*(mean(data_new) - Kmeans_results_7index{1,k}.C(clu_num,:))*(mean(data_new) - Kmeans_results_7index{1,k}.C(clu_num,:))';
     end
     VRS(k) = SSB(k)/Tol_CSS(k)*(length(data_new)-k)/(k-1);
end

% subtype_result
ori_res_sub1=mean(response(find(Kmeans_results_7index{1,2}.IDX==1)));
ori_res_sub2=mean(response(find(Kmeans_results_7index{1,2}.IDX==2)));
subtype_result_7index(1,1)=length(find(Kmeans_results_7index{1,2}.IDX==1)); %number of subtype 1
subtype_result_7index(1,3)=length(find(Kmeans_results_7index{1,2}.IDX==2)); %number of subtype 1
subtype_result_7index(1,2)=ori_res_sub1; %response rate of subtype 1
subtype_result_7index(1,4)=ori_res_sub2; %response rate of subtype 2
