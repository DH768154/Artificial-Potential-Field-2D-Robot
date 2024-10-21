clc;clear;close all
addpath('myfunc','-begin')
ag = linspace(0,360,361);
C_obs{1} = 0.6*[cosd(ag);sind(ag)]+[1;3.5];

L = [1,1,1,1,1];
n = length(L);
links0 = [0,L;zeros(1,n+1)];
%%
n_ctrlpt = 4;
ctrlpts0 = reshape(L.*linspace(0,1,n_ctrlpt+1)',n_ctrlpt+1,1,[]);
ctrlpts0 = pagectranspose([ctrlpts0,zeros(size(ctrlpts0))]);
ctrlpts0 = ctrlpts0(:,2:end,:);
ctrlpts0_all = reshape(ctrlpts0,2,[],1);
%%
qI = [10,20,30,30,20]'/180*pi;
qG = [120,8,7,6,5]'/180*pi;

[frames,linksG] = KineSimpleRob(links0,qG);
%%
f = figure;
patch(C_obs{1}(1,:),C_obs{1}(2,:),'y'); hold on
patch(C_obs{1}(1,:),C_obs{1}(2,:),'y'); hold on
plot(linksG(1,:),linksG(2,:),'.--r','LineWidth',1); hold on
pl = plot(NaN,NaN,'.-k','LineWidth',1,'MarkerSize',14); hold on
pc = plot(NaN,NaN,'.r','MarkerSize',10); hold on
pr = quiver(NaN,NaN,NaN,NaN,'Color','b','LineWidth',1); hold on
plot([-1,1]/2,[0,0],'-k','LineWidth',3); hold on
legend([pc,pr],{'ctrl pts','repulse'})
grid on; axis equal
xlim([-5,5]); ylim([-0.2,5]);
title('Attractive Force applied in C-Space')
set(f,'Units','normalized','Position',[0.175,0.175,0.65,0.65]);

coeff_att = [50,30,3,1,1]';
maxdist = 10;
coeff_rep = 0.01;
rad_rep = 1.2;
stepsize = 0.01;
max_step = [1,1,1,1,1]'*0.5/180*pi.*[-1,1];

qC = qI;
for kk = 1:1000
    [frames,links] = KineSimpleRob(links0,qC);
    ctrlpts = NaN(2,n_ctrlpt,n);
    for i = 1:n
        ctrlpts(:,:,i) = applySE2(frames(:,:,i),ctrlpts0(:,:,i));
    end
    ctrlpts_all = reshape(ctrlpts,2,[],1);
%% 
% 

    f = NaN(size(ctrlpts));

    tq_rep = zeros(n,n_ctrlpt,n);
    for i = 1:n % linknum
        for j = 1:n_ctrlpt % j-th control point
            ccpt = ctrlpts(:,j,i);
            J = NaN(2,i);
            for ii = 1:i
                P = ctrlpts(:,j,i)-links(:,ii);
                J(:,ii) = [-P(2);P(1)];
            end

            [closestpt,closestdist] = closestpt2Shape2d(ctrlpts(:,j,i),C_obs);
            f(:,j,i) = dUrepFunc(ctrlpts(:,j,i),coeff_rep,rad_rep,closestdist,closestpt);
            tq_rep(1:i,j,i) = J' * f(:,j,i);
        end
    end
    tq_rep_all = sum(tq_rep,2);
    tq_rep_all = sum(tq_rep_all,3);

     f_rep = reshape(f,2,[],1);
%     f_rep_all = sum(f_rep,2);
    
    dist45 = norm(qC(4:end)-qG(4:end));
    if all(abs(qC(1:3)-qG(1:3))<=5/180*pi)
        coeff_att_new = coeff_att;
        coeff_att_new(4:end) = coeff_att(4:end)*5;
    else
        coeff_att_new = coeff_att;
    end
    tq_att = dUattFunc(qC,qG,coeff_att_new,maxdist);


    dq = stepsize * (tq_att+tq_rep_all);
    dq = bound2range(dq,max_step);
    qC = qC-dq;
    if norm(qC-qG)<=norm([1,1,1,1,1]/5/180*pi)
        break
    end
%% 
% 

    pl.XData = links(1,:);
    pl.YData = links(2,:);
    pr.XData = ctrlpts_all(1,:);
    pr.YData = ctrlpts_all(2,:);
    pr.UData = -f_rep(1,:);
    pr.VData = -f_rep(2,:);
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