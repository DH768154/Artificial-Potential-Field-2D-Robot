clc;clear;close all
addpath('myfunc','-begin')

ag = linspace(0,360,361);
C_obs{1} = 0.6*[cosd(ag);sind(ag)]+[1;3.5];
%C_obs{1} = 0.6*[cosd(ag);sind(ag)]+[4;3];

L = [1,1,1,1,1];
n = length(L);
links0 = [0,L;zeros(1,n+1)];
%%
n_ctrlpt = 4;
ctrlpts0 = reshape(L.*linspace(0,1,n_ctrlpt+1)',n_ctrlpt+1,1,[]);
ctrlpts0 = pagectranspose([ctrlpts0,zeros(size(ctrlpts0))]);
ctrlpts0 = ctrlpts0(:,2:end,:);
ctrlpts0_all = reshape(ctrlpts0,2,[],1);

[frames,~] = KineSimpleRob(links0,[0,0,0,0,0]');
ctrlptsH = NaN(2,n_ctrlpt,n);
for i = 1:n
    ctrlptsH(:,:,i) = applySE2(frames(:,:,i),ctrlpts0(:,:,i));
end
w = 1./ctrlptsH(1,:,:).^2;
%w = ones(size(w));
%%
qI = [10,20,30,30,20]'/180*pi;
qG = [120,8,7,6,5]'/180*pi;

[frames,linksG] = KineSimpleRob(links0,qG);
ctrlptsG = NaN(2,n_ctrlpt,n);
for i = 1:n
    ctrlptsG(:,:,i) = applySE2(frames(:,:,i),ctrlpts0(:,:,i));
end
%%
f = figure;
patch(C_obs{1}(1,:),C_obs{1}(2,:),'y'); hold on
plot(linksG(1,:),linksG(2,:),'.--r','LineWidth',1); hold on
pl = plot(NaN,NaN,'.-k','LineWidth',1,'MarkerSize',8); hold on
pc = plot(NaN,NaN,'.r','MarkerSize',10); hold on
pr = quiver(NaN,NaN,NaN,NaN,'Color','b','LineWidth',1); hold on
pa = quiver(NaN,NaN,NaN,NaN,'Color',[0.4660 0.6740 0.1880],'LineWidth',1); hold on
plot([-1,1]/2,[0,0],'-k','LineWidth',3); hold on
legend([pc,pr,pa],{'ctrl pts','repulse','attract'})
grid on; axis equal
xlim([-5,5]); ylim([-0.2,5]);
title('Attractive Force applied on Control Points in W-Space')
set(f,'Units','normalized','Position',[0.175,0.175,0.65,0.65]);

coeff_att = 10;%[20,20,5,1,1]';
maxdist = 10;
coeff_rep = 1;
rad_rep = 1;
stepsize = 0.001;
max_step = [1,1,1,1,1]'*0.5/180*pi.*[-1,1];
max_iter = 1000;

qC = qI;
for kk = 1:max_iter
    [frames,links] = KineSimpleRob(links0,qC);
    ctrlpts = NaN(2,n_ctrlpt,n);
    for i = 1:n
        ctrlpts(:,:,i) = applySE2(frames(:,:,i),ctrlpts0(:,:,i));
    end
    if mean(sqrt(sum((ctrlpts-ctrlptsG).^2)),'all')<=0.2
        break
    end
    ctrlpts_all = reshape(ctrlpts,2,[],1);
%% 
% 

    f_rep = NaN(size(ctrlpts));
    f_att = NaN(size(ctrlpts));

    tq_rep = zeros(n,n_ctrlpt,n);
    tq_att = zeros(n,n_ctrlpt,n);
    for i = 1:n % linknum
        for j = 1:n_ctrlpt % j-th control point
            ccpt = ctrlpts(:,j,i);
            J = NaN(2,i);
            for ii = 1:i
                P = ctrlpts(:,j,i)-links(:,ii);
                J(:,ii) = [-P(2);P(1)];
            end

            [closestpt,closestdist] = closestpt2Shape2d(ctrlpts(:,j,i),C_obs);
            f_rep(:,j,i) = dUrepFunc(ctrlpts(:,j,i),coeff_rep,rad_rep,closestdist,closestpt);
            tq_rep(1:i,j,i) = J' * f_rep(:,j,i);
            f_att(:,j,i) = dUattFunc(ctrlpts(:,j,i),ctrlptsG(:,j,i),coeff_att,maxdist);

            %if mean(sqrt(sum((ctrlpts-ctrlptsG).^2)),'all')>0.5
            f_att(:,j,i) = f_att(:,j,i) .* w(:,j,i);
%             else
%                 stepsize=0.001;
% %                 if j~=n_ctrlpt || i~=n
% %                     f_att(:,j,i) = 0;
% %                 end
%             end

            tq_att(1:i,j,i) = J' * f_att(:,j,i);
        end
    end
    tq_all = sum(tq_rep+tq_att,2);
    tq_all = sum(tq_all,3);

     f_rep = reshape(f_rep,2,[],1);
%     f_rep_all = sum(f_rep,2);
    

    %tq_att = dUattFunc(qC,qG,coeff_att,maxdist);


    dq = stepsize * tq_all;
    dq = bound2range(dq,max_step);
    qC = qC-dq;
%% 
% 

    pl.XData = links(1,:);
    pl.YData = links(2,:);
    pr.XData = ctrlpts_all(1,:);
    pr.YData = ctrlpts_all(2,:);
    pr.UData = -f_rep(1,:);
    pr.VData = -f_rep(2,:);
    pa.XData = ctrlpts_all(1,:);
    pa.YData = ctrlpts_all(2,:);
    pa.UData = -f_att(1,:);
    pa.VData = -f_att(2,:);

    pc.XData = ctrlpts_all(1,:);
    pc.YData = ctrlpts_all(2,:);
    drawnow;
end
%%


%%
function [frames,links] = KineSimpleRob(links0,ag)
n = length(ag);
T = NaN(3,3,n+1);
for i = 1:n
    T(:,:,i) = [eye(2),links0(:,i);0,0,1]*Rot2d(ag(i),[0;0]);
end
T(:,:,end) = [eye(2),links0(:,i);0,0,1];
frames = cummult(T);

links = reshape(frames(1:2,3,:),2,[],1);
end
% toc
%% 
% 

function pt = applySE2(T,pt)
pt = T*[pt;ones(1,size(pt,2))];
pt = pt(1:2,:);
end
%% 
% 

function F = dUrepFunc(q,coeff_rep,rad_rep,closestdist,closestpt)
ind = closestdist<rad_rep;
F = zeros(length(q),length(closestdist));
if any(ind)
    F(:,ind) = coeff_rep*(1/rad_rep-1./closestdist(ind))./closestdist(ind).^3.*(q-closestpt(:,ind));
end
end
%% 
% 

function F = dUattFunc(q,qG,coeff_att,maxdist)
dq = q-qG;
dist = norm(dq);
if dist<=maxdist
    F = coeff_att.*dq;
else
    F = maxdist*coeff_att.*dq/dist;
end
end