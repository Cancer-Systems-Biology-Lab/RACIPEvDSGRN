params = {'Inh_of_GRHL2ToZEB', 'Num_of_GRHL2ToZEB', 'Trd_of_GRHL2ToZEB', 'Inh_of_ZEBToGRHL2', 'Num_of_ZEBToGRHL2', 'Trd_of_ZEBToGRHL2', 'Inh_of_miR200ToZEB', 'Num_of_miR200ToZEB', 'Trd_of_miR200ToZEB', 'Inh_of_ZEBTomiR200', 'Num_of_ZEBTomiR200', 'Trd_of_ZEBTomiR200', 'Inh_of_SLUGTomiR200', 'Num_of_SLUGTomiR200', 'Trd_of_SLUGTomiR200', 'Act_of_SLUGToZEB', 'Num_of_SLUGToZEB', 'Trd_of_SLUGToZEB', 'Act_of_ZEBToZEB', 'Num_of_ZEBToZEB', 'Trd_of_ZEBToZEB', 'Act_of_SLUGToSLUG', 'Num_of_SLUGToSLUG', 'Trd_of_SLUGToSLUG', 'Prod_of_GRHL2', 'Deg_of_GRHL2', 'Prod_of_ZEB', 'Deg_of_ZEB', 'Prod_of_miR200', 'Deg_of_miR200', 'Prod_of_SLUG', 'Deg_of_SLUG'};
prs_file = 'GRHL_52To61.prs';
prs_new = 'GRHL_52To61.dat';
copyfile(prs_file, prs_new);
par_list = readtable(prs_new);
par_list = string(par_list.Parameter);
par_list = arrayfun(@(x) replace(x, '-', ''),par_list);
par_order = zeros(length(par_list),1);
for i = 1:length(params)
    par_order(i) = find(par_list == params(i));
end
n_nodes = sum(contains(par_list, 'Prod'));
par_file = 'GRHL_52To61_parameters.dat';
parameters = table2array(readtable(par_file));
parameters =  parameters(:, 3:end);
initFile = 'GRHL_52To61_parameters_eqval.dat';
fid = fopen(initFile);
initConds = strings;
tline = fgetl(fid);
while ischar(tline)
    initConds = [initConds,tline];
    tline = fgetl(fid);
end
fclose(fid);
tspan = 1:100;
ss_tot = zeros(1, n_nodes+2);
parfor i = 1:size(parameters,1)
    p = parameters(i,par_order);
    p1 = parameters(i,:);
    ss=zeros(1,n_nodes);
    initDat = initConds(i+1);
    if initDat == '-1'
       continue;
    end
    initDat = split(initDat);
    l = length(initDat);
    for j = 1:l
        y0 = str2num(initDat(j));
        [t y] = ode23s(@(t,y)GRHL_52To61(t,y,p), tspan, y0);
        while sum(round((y(end, :) - y((end-1), :)), 3)) ~= 0
            y0 = y(end,:);
            [t y] = ode23s(@(t,y)GRHL_52To61(t,y,p), tspan, y0);
        end
        %if all(abs(sum(ss - y(end, :),2)) > 0.01)
        ss = [ss;round(y(end, :),3)];
        %end
    end
    ss = ss(2:end, :);
    [aa,~,c] = unique(ss, 'rows');
    ss = [aa, histcounts(c,1:max(c)+1)'];
    ss = [repelem(i, size(ss,1))' ss];
    ss_tot = [ss_tot; ss];
end
ss_tot = ss_tot(2:end,:);
ss_tot = array2table(ss_tot);
ss_tot.Properties.VariableNames = {'ParIndex' 'GRHL2' 'ZEB' 'miR200' 'SLUG' 'Count'};
writetable(ss_tot, 'GRHL_52To61_ss.csv');
