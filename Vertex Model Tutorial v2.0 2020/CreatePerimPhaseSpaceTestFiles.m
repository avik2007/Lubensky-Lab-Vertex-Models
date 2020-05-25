%%MAKE SURE YOU PUT DETAILED DESCRIPTION OF WHAT YOU ARE DOING HERE 
%%
K = 10;
%%
xlim = 3; % number of x coordinates to check 
ylim = 4; % number of y coordinates to check

xcoords = zeros(xlim,1); %contains x coordinates of checked points in phase spcae
ycoords = zeros(ylim,1); %contains y coordinates of checked points in phase space

xint = 0.1; % precision with which space is searched in x-direction
yint = 0.5; % precision with which space is searched in y-direction

% x is LAMBDA, from Meryl's vertex model
% y is GAMMA, the one from LISA MANNINGS paper (see GenerationVoronoi_Lloyd
%Iteration 2/26/20 comment
%% create run files for testing vertex model
% this is the 2/26 tested region of phase space -COPY AND PASTE THIS CODE TO EXPAND THE REGION 
xp = zeros(xlim*ylim,1); %for plotting phase points
yp = zeros(xlim*ylim,1);

for indexLA = 1:xlim
    for indexGA = 1:ylim
        xcoords(indexLA,1) = 0.1 + xint*indexLA;
        ycoords(indexGA,1) = -2  + yint*indexGA;
        xc = xcoords(indexLA,1);
        yc = ycoords(indexGA,1);
        nf = strcat('LNbC_LA_', num2str(xc,'%1.1f'), '_GA_', num2str(yc,'%1.1f'),  '.txt');
        %GenerationVoronoi_LloydIteration_PerimPhaseTests(2000, 20, 30,nf , xcoords(indexLA,1), ycoords(indexGA,1));
        a = (indexLA-1) * (ylim)+ indexGA;
        b = (indexLA-1) * (ylim) + indexGA;
        xp(a,1) = xc;
        yp(b,1) = yc;
    end 
end
%ALL OF THESE COMPILED AND RAN WITHOUT ERRORS
%% do phase space comparisons between Meryl's and Manning's model
xprel1 = xp/K;
yprel1 = yp/K;
plot(xprel1,yprel1,'o', 'MarkerFaceColor', 'b');% Compiled points in phase space can be solid, others hollow
xlabel('$$\frac{\Lambda}{K}$$', 'interpreter', 'latex');
ylabel('$$\frac{\gamma}{K}$$','interpreter', 'latex');
title('$(\frac{\gamma}{K}, \frac{\Lambda}{K})$ phase plot ', 'interpreter', 'latex');

%% I want to try a significantly larger section of phase space
%% create run files for testing vertex model
% this is the 12/27 tested region of phase space -COPY AND PASTE THIS CODE TO EXPAND THE REGION 

%x refers to LAMBDA/K
%y refers to GAMMA/K
xlim = 4; % number of x coordinates to check 
ylim = 4; % number of y coordinates to check

xcoords = zeros(xlim,1); %contains x coordinates of checked points in phase spcae
ycoords = zeros(ylim,1); %contains y coordinates of checked points in phase space

xmax = 1;
ymax = 1;
xmin = 0;
ymin = -1;

xint = (xmax - xmin) / (xlim); % precision with which space is searched in x-direction
yint = (ymax - ymin) / (ylim); % precision with which space is searched in y-direction


xprel2 = zeros((xlim+1)*(ylim+1),1); %for plotting phase points
yprel2 = zeros((xlim+1)*(ylim+1),1);

samplelist = struct('name',{});
for indexLA = 1:xlim+1
    for indexGA = 1:ylim+1
        xcoords(indexLA,1) = xmin + xint*(indexLA-1);
        ycoords(indexGA,1) = ymin  + yint*(indexGA-1);
        xc = xcoords(indexLA,1);
        yc = ycoords(indexGA,1);
        nf = strcat('LNbC_LA_', num2str(xc,'%1.2f'), '_GA_', num2str(yc,'%1.2f'),  '.txt');
        %GenerationVoronoi_LloydIteration_PerimPhaseTests(2000, 20, 30,nf , K*xcoords(indexLA,1), K*ycoords(indexGA,1));
        a = (indexLA-1) * (ylim + 1)+ indexGA;
        samplelist(end+1) = struct('name',{nf});
        xprel2(a,1) = xc;
        yprel2(a,1) = yc;
    end 
end
%% Classify phase space points from 12/27 test region

%You have 4 types of phase points - nonshrinking, nonbreaking; shrinking,
%breaking; shrinking, nonbreaking; nonshrinking, breaking
%We will call them nn, sb, sn, nb

nnindex = find(~[samplelist.shrink] & ~[samplelist.errors]); %these will be labelled as blue circles
sbindex = find([samplelist.shrink] & [samplelist.errors]); %these will be green squares's
snindex = find([samplelist.shrink] & ~[samplelist.errors]); % these will be blue squares's
nbindex = find(~[samplelist.shrink] & [samplelist.errors]); % these will be green circles
%new fields - shrink, errors

plot(xprel2(nnindex,1), yprel2(nnindex,1),'o', 'MarkerFaceColor', 'b');
hold on;
plot(xprel2(sbindex,1), yprel2(sbindex,1),'s', 'MarkerFaceColor', 'g');
plot(xprel2(snindex,1), yprel2(snindex,1),'s', 'MarkerFaceColor', 'b');
plot(xprel2(nbindex,1), yprel2(nbindex,1),'o', 'MarkerFaceColor', 'g');
plot(xprel1,yprel1,'s', 'MarkerFaceColor', 'b');% Compiled points in phase space can be solid, others hollow
xlabel('$$\frac{\Lambda}{K}$$', 'interpreter', 'latex');
ylabel('$$\frac{\gamma}{K}$$','interpreter', 'latex');
title('$(\frac{\gamma}{K}, \frac{\Lambda}{K})$ phase plot (ignoring units), Lloyd iteration 30', 'interpreter', 'latex');
legend('no shrink, no break','shrinks, breaks', 'shrinks, no break', 'no shrink, break');
hold off;


