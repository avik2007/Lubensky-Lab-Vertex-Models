    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script produces a txt file with the initial conditions needed to run the 
% vertex model cpp exicutable. The Global values are set by the user. The
% cell positions are a honeycomb tilling of the plane. The resulting txt file
% contians all the global parameter values and cell, edge, and vertex
% locatons and neighbors. 
%
% To use this file specify the name of the output txt file and the number
% of cells you would like in your tissue. You may wish to modify some
% global parameters to fit your needs. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

%% Inputs needed
INITIALSTEPSIZE = 0.0001;    %Initial stepsize for the numberic integrator 
nomfichier=strcat('shearmodulus_check0.95','.txt');  %the file name
n=1;                       %1/2 of the Cells per row
m=2;                       %Cells per collumn 
r=89542;                    %Seed for the random number generators

%% Global Parameters 
GAMMA =10;                  %Tension on an edge
l= sqrt(2) * 25;            %Average edge length    
t0 = 1;                     %Standard unit of time 
RLXLY = 1;                  %Ratio of LX to LY (used to make rectangular boxes) 
Ratio = 0.3;                %ratio of pressure to tension forces (standard value 0.3)
Ralpha=1000000000000000000000000000;                   %Ralpha=10 means the box evolves 10 times slower then the vertices
TOTALTIME =	80;             %How long to run the simulation 
CURRENTTIME = 0.0;          %The initial time
TIMEINTERVAL =	0.1;        %How often to output the current configuration of the tissue
%INITIALSTEPSIZE = 0.002;    %Initial stepsize for the numberic integrator 
error_tolerance_absolute_mechanical = 0.01;  %Error tolerence for use with the adaptive stepsize integrators 
LX=n*3*l+1;
LY=m*sqrt(3)*l+1; 
NCells=n*m*2;
A0 = LX * LY / NCells;%4033.1139     ;                    %LX*LY/NCells;                   %Cells want to have area A0 if box size can change
%Here, I'm making a strain test for the hexagonal lattice by calculating
%what the preferred area should be for the initial area of all the cells to
% the equilibrium area.
MINAREA = A0/100;                    %Minimum allowed cell area (smaller cells will be removed)
MINLENGTH =	l/40;                    %Minimum edge length (smaller edges will be merged to form a fourfold vertex
NEWEDGELENGTH =	l/20;                %Length of a new edge created in the T1 process%
%K = Ratio * GAMMA / l^3 ;           %for pressure=K(A-A0)
K = Ratio * 10 * GAMMA / l;          %for pressure=K(1-A/A0)
Alpha_M = GAMMA*t0/l;                % Drag coefficent
TISSUE_ALPHA = 10000000000000000000000000000000 * Alpha_M;%10000000*Ralpha*Alpha_M;       % How fast the box size updates relative to the vertex motion
TISSUE_EXT_STRESS_Y = 0;             % External Stress on the tissue (positive values expancile) 
TISSUE_EXT_STRESS_X = 0;
P0 = 0.*(sqrt(24 / sqrt(3))) * sqrt(A0); % sqrt(24 / sqrt(3)) = regular hexagon shape parameter
LAMBDA = 0.01 * K * A0;%K / 0.05; % mimics Manning Model 
VERBOSEON = 1; %boolean, tells code to output extra data for cells and edges when true
equilibriumtolerance = 0.00000000001;

%% Create the Hexagonal tillining 
SEED = r;
rng(SEED); 
f2=figure();
offset=0;%0.01*l; 
Points1=zeros(NCells,2); 
k=1; 
strainx = 0;%.5;

%Grid One
for i=0:n-1
    for j=0:m-1
        Points1(k,:)=[(3*i*l +rand()*offset)*(1 +strainx) ,sqrt(3)*j*l+rand()*offset]; 
        k=k+1; 
    end
end
% Grid two

for i=0:n-1
    for j=0:m-1
        Points1(k,:)=[((3*i+3/2)*l +rand()*offset)* (1 + strainx) , (sqrt(3)*j+sqrt(3)/2)*l+rand()*offset]; 
        k=k+1; 
    end
end

LX = LX * (1 + strainx);
Pointsa=Points1+[LX *ones(NCells,1) LY*zeros(NCells,1)];
Pointsb=Points1-[LX*zeros(NCells,1) LY*ones(NCells,1)];
Pointsc=Points1-[LX*ones(NCells,1) LY*zeros(NCells,1)];
Pointsd=Points1+[LX*zeros(NCells,1) LY*ones(NCells,1)];
Pointse=Points1+[LX*ones(NCells,1) LY*ones(NCells,1)];
Pointsf=Points1+[LX*ones(NCells,1) -LY*ones(NCells,1)];
Pointsg=Points1+[-LX*ones(NCells,1) LY*ones(NCells,1)];
Pointsh=Points1+[-LX*ones(NCells,1) -LY*ones(NCells,1)];
Pointst=[Points1; Pointsa; Pointsb; Pointsc; Pointsd; Pointse; Pointsf; Pointsg; Pointsh];
[v,c]=voronoin(Pointst); 


%% Plot a Representation of the result
set(0,'CurrentFigure',f2);
plot(Points1(:,1),Points1(:,2),'r+');
for i = 1:length(c) 
    %ca are the cells that are in the interior only
    IN=(all(v(c{i},1)>0)&all(v(c{i},1)<LX)&all(v(c{i},2)>0)&all(v(c{i},2)<LY));
    %ca have cells that are at the border
    IN2=(any(v(c{i},1)>0)&any(v(c{i},1)<LX)&any(v(c{i},2)>0)&any(v(c{i},2)<LY))&(~IN);
    if (all(c{i}~=1)&&IN)  
    set(0,'CurrentFigure',f2)
    patch(v(c{i},1),v(c{i},2),'y'); % use color yellow.
    end
    % cells crossing the border
    if (all(c{i}~=1)&&(IN2))  
    set(0,'CurrentFigure',f2)
    patch(v(c{i},1),v(c{i},2),'g'); % use color green.
    end
end
set(0,'CurrentFigure',f2)
axis equal
xlim([0 LX]);
ylim([0 LY]);


for i = 1:length(c)
   INGLOBAL=any((v(c{i},1)>0)&v(c{i},1)<LX&v(c{i},2)>0&v(c{i},2)<LY);
   if INGLOBAL
       ListVertex=c{i};
       area=0;
       for k=1:length(c{i})
           kplus=mod(k,length(c{i}))+1;
           area=area+v(c{i}(k),1)*v(c{i}(kplus),2)-v(c{i}(kplus),1)*v(c{i}(k),2);
       end
       if(area<0) c{i}=fliplr(c{i}); end
   end
end

%% Export to a txt file
fid=fopen(nomfichier,'wt');
fprintf(fid, 'Voronoi Tesselation Data: \n');
fprintf(fid, 'NumbCells %d \t rng %d \n', NCells, r);
fprintf(fid, 'ratio pressure to tension\t%f\n', Ratio);
fprintf(fid, 'Global:\n');
fprintf(fid, 'LX\t%f\nLY\t%f\n', LX, LY);
fprintf(fid, 'A0\t%f\n', A0);
fprintf(fid, 'SEED\t%d\n', SEED);
fprintf(fid, 'TIMEINTERVAL\t%f\n', TIMEINTERVAL);
fprintf(fid, 'MINAREA\t%f\n', MINAREA);
fprintf(fid, 'MINLENGTH\t%f\n', MINLENGTH);
fprintf(fid, 'NEWEDGELENGTH\t%f\n', NEWEDGELENGTH);
fprintf(fid, 'K\t%f\n', K);
fprintf(fid, 'GAMMA\t%f\n', GAMMA);
fprintf(fid, 'LAMBDA\t%f\n', LAMBDA);
fprintf(fid, 'P0\t%f\n', P0);
fprintf(fid, 'TOTALTIME\t%f\n', TOTALTIME);
fprintf(fid, 'CURRENTTIME\t%f\n', CURRENTTIME);
fprintf(fid, 'VERBOSEON\t%d\n', VERBOSEON);
fprintf(fid, 'INITIALSTEPSIZE\t%f\n', INITIALSTEPSIZE);
fprintf(fid, 'error_tolerance_absolute_mechanical\t%f\n', error_tolerance_absolute_mechanical);
fprintf(fid, 'TISSUE_ALPHA\t%f\n', TISSUE_ALPHA);
fprintf(fid, 'ALPHAM\t%f\n', Alpha_M);
fprintf(fid, 'TISSUE_EXT_STRESS_X\t%f\n', TISSUE_EXT_STRESS_X);
fprintf(fid, 'TISSUE_EXT_STRESS_Y\t%f\n', TISSUE_EXT_STRESS_Y);
fprintf(fid, 'equilibriumtolerance\t%d\n', equilibriumtolerance);

fprintf(fid, 'endGlobal\n'); 
fprintf(fid, 'Vertices\n');
fprintf(fid, 'id\tx\ty\n');

NbVertex=0;
CentralVertices=zeros(length(v),2);
for i = 2:length(v)   % v(1) is always a point at infinity, so skip (NOTED ADDED BY DKL 7/5/10)
    % the vertices don't necessarily have to be in order because the code
    % will skip over vertices that are outside of the window of interest
    % (NOTE ADDED BY ABM 10/14/23
   if(v(i,1)>0&&v(i,1)<LX&&v(i,2)>0&&v(i,2)<LY)
       fprintf(fid,'%d',i); 
       fprintf(fid,'\t %2.4f',v(i,1));
       fprintf(fid,'\t %2.4f',v(i,2));
       fprintf(fid,'\n');
       NbVertex=NbVertex+1;
       CentralVertices(NbVertex,1)=v(i,1);
       CentralVertices(NbVertex,2)=v(i,2);
   end
end

CentralVertices=CentralVertices(1:NbVertex,1:2);
EqVertex=zeros(NbVertex,1);

for i = 2:length(v)
   vcorrige=zeros(1,2);
   vtest=v(i,:);

   if (vtest(1)>0&&vtest(1)<LX&&vtest(2)>0&&vtest(2)<LY) vcorrige=vtest;end
   vtest=v(i,:)-LX*[1 0];
   if (vtest(1)>0&&vtest(1)<LX&&vtest(2)>0&&vtest(2)<LY) vcorrige=vtest;  end
   vtest=v(i,:)-LY*[0 1];
   if (vtest(1)>0&&vtest(1)<LX&&vtest(2)>0&&vtest(2)<LY) vcorrige=vtest;  end
   vtest=v(i,:)+LX*[1 0];
   if (vtest(1)>0&&vtest(1)<LX&&vtest(2)>0&&vtest(2)<LY) vcorrige=vtest;  end
   vtest=v(i,:)+LY*[0 1];
   if (vtest(1)>0&&vtest(1)<LX&&vtest(2)>0&&vtest(2)<LY) vcorrige=vtest;  end
   vtest=v(i,:)-[LX LY];
   if (vtest(1)>0&&vtest(1)<LX&&vtest(2)>0&&vtest(2)<LY) vcorrige=vtest;  end
   vtest=v(i,:)-[LX -LY];
   if (vtest(1)>0&&vtest(1)<LX&&vtest(2)>0&&vtest(2)<LY) vcorrige=vtest;  end
   vtest=v(i,:)+[LX LY];
   if (vtest(1)>0&&vtest(1)<LX&&vtest(2)>0&&vtest(2)<LY) vcorrige=vtest;  end
   vtest=v(i,:)+[LX -LY];
   if (vtest(1)>0&&vtest(1)<LX&&vtest(2)>0&&vtest(2)<LY) vcorrige=vtest;  end
   column=find(all(abs((repmat(vcorrige,size(v,1),1))'-v')<0.000001));%strmatch(vcorrige, v);
   if(~isempty(column))
       EqVertex(i)=column;
   end

end
 
fprintf(fid, 'Edges\n');
fprintf(fid, 'id\tvertex1\tvertex2\txflag\tyflag\n');
NbEdgeMax=3*NbVertex;
Edges.Vertices=zeros(NbEdgeMax,2);
Edges.Crossx12=zeros(NbEdgeMax,1);
Edges.Crossx21=zeros(NbEdgeMax,1);
Edges.Crossy12=zeros(NbEdgeMax,1);
Edges.Crossy21=zeros(NbEdgeMax,1);

EdgeIncrement=0; %the number of the edge so you get 1,2,3,4 as edge lables
for i = 1:length(c) %For each vertex
   INGLOBAL=any( (v(c{i},1)>0) & v(c{i},1)<LX & v(c{i},2)>0 & v(c{i},2)<LY);
   if INGLOBAL %If it is in the reigon of intrest
       ListVertex=c{i}; %list of edges in each cell in the reigon 
       for j=1:length(ListVertex) %for each edge in the cell
           testdejainscrit=0;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
           % if we look at the current edge is in the list
           for k=1:EdgeIncrement
               if ((EqVertex((ListVertex(j)))==Edges.Vertices(k,1)...
                       &&EqVertex(ListVertex(mod(j,length(ListVertex))+1))==Edges.Vertices(k,2))...%the modulus is to be on that loop
                       ||...
                       (EqVertex(ListVertex(j))==Edges.Vertices(k,2)...
                       &&EqVertex(ListVertex(mod(j,length(ListVertex))+1))==Edges.Vertices(k,1)))
                   testdejainscrit=1;
               else
               end
           end
           % if it is not there we added at the end
           if(~testdejainscrit)
               EdgeIncrement=EdgeIncrement+1;
               Edges.Vertices(EdgeIncrement,1)=EqVertex(ListVertex(j)); %the first vertex
               Edges.Vertices(EdgeIncrement,2)=EqVertex(ListVertex(mod(j,length(ListVertex))+1)); %the second vertex
               Vertex1encours=ListVertex(j);
               Vertex1encourseq=EqVertex(ListVertex(j));
               Vertex2encours=ListVertex(mod(j,length(ListVertex))+1);
               Vertex2encourseq=EqVertex(ListVertex(mod(j,length(ListVertex))+1));
               % and is written in the file
               fprintf(fid,'%d',EdgeIncrement);
               fprintf(fid,'\t %d', Edges.Vertices(EdgeIncrement,1));
               fprintf(fid,'\t %d', Edges.Vertices(EdgeIncrement,2));
               % must be a test to see if the one-vertex through a
               % rating (a vertex in, one out) 2-if there are two vertices
               % out but through the side anyway (happens on the corner)
               VertexCross21=(EqVertex(ListVertex(j))~=ListVertex(j))&&(EqVertex(ListVertex(mod(j,length(ListVertex))+1))==ListVertex(mod(j,length(ListVertex))+1));
               VertexCross12=(EqVertex(ListVertex(j))==ListVertex(j))&&(EqVertex(ListVertex(mod(j,length(ListVertex))+1))~=ListVertex(mod(j,length(ListVertex))+1));
               % cas 1

               %  x12=x21=y12=y21=0;

               %Periodic boundary condition testing
               if((EqVertex(ListVertex(j))~=ListVertex(j))&&(EqVertex(ListVertex(mod(j,length(ListVertex))+1))==ListVertex(mod(j,length(ListVertex))+1)))
                   if(v(ListVertex(j),1)<0) Edges.Crossx12(EdgeIncrement)=1;%fprintf(fid,' xCrossBoundary12');
                       %    x12=1;
                   end
                   if(v(ListVertex(j),1)>LX) Edges.Crossx21(EdgeIncrement)=1;%fprintf(fid,' xCrossBoundary21');
                       %   x21=1;
                   end
                   if(v(ListVertex(j),2)<0) Edges.Crossy12(EdgeIncrement)=1;%fprintf(fid,' yCrossBoundary12');
                       %   y12=1;
                   end
                   if(v(ListVertex(j),2)>LY) Edges.Crossy21(EdgeIncrement)=1;%fprintf(fid,' yCrossBoundary21');
                       %    y21=1;
                   end
               end
               if((EqVertex(ListVertex(j))==ListVertex(j))&&(EqVertex(ListVertex(mod(j,length(ListVertex))+1))~=ListVertex(mod(j,length(ListVertex))+1)))
                   if(v(ListVertex(mod(j,length(ListVertex))+1),1)<0) 
                       Edges.Crossx21(EdgeIncrement)=1;
                   end
                   if(v(ListVertex(mod(j,length(ListVertex))+1),1)>LX) 
                       Edges.Crossx12(EdgeIncrement)=1;
                   end
                   if(v(ListVertex(mod(j,length(ListVertex))+1),2)<0) 
                       Edges.Crossy21(EdgeIncrement)=1;
                   end
                   if(v(ListVertex(mod(j,length(ListVertex))+1),2)>LY) 
                       Edges.Crossy12(EdgeIncrement)=1;
                   end
               end
               %cas 2
               if((EqVertex(ListVertex(j))~=ListVertex(j))&&(EqVertex(ListVertex(mod(j,length(ListVertex))+1))~=ListVertex(mod(j,length(ListVertex))+1)))
                   % if the vertex is in the range which is not
                   % moving in the same direction, is that both
                   % ends of the edge are not contained in the same
                   % in the time interval.
                   v1=v(ListVertex(j),:);
                   v2=v(ListVertex(mod(j,length(ListVertex))+1),:);
                   eqv1=v(EqVertex(ListVertex(j)),:);
                   eqv2=v(EqVertex(ListVertex(mod(j,length(ListVertex))+1)),:);
                   %Periodic boundary condition chehck again
                   if(v2(1)>v1(1))&&(eqv2(1)<eqv1(1)) 
                       Edges.Crossx12(EdgeIncrement)=1;
                   end
                   if(v2(1)<v1(1))&&(eqv2(1)>eqv1(1)) 
                       Edges.Crossx21(EdgeIncrement)=1;
                   end
                   if (v2(2)>v1(2))&&(eqv2(2)<eqv1(2))
                       Edges.Crossy12(EdgeIncrement)=1;
                   end
                   if (v2(2)<v1(2))&&(eqv2(2)>eqv1(2)) 
                       Edges.Crossy21(EdgeIncrement)=1;
                   end
               end
              
               %New output added by MS 4/17/2014
               if ( (Edges.Crossx12(EdgeIncrement)==1) || (Edges.Crossx21(EdgeIncrement)==1) )
                   newxflag=1;
               end
               if ( (Edges.Crossx12(EdgeIncrement)==0) && (Edges.Crossx21(EdgeIncrement)==0) )
                   newxflag=0;
               end
               if ( (Edges.Crossy12(EdgeIncrement)==1) || (Edges.Crossy21(EdgeIncrement)==1) )
                   newyflag=1;
               end
               if ( (Edges.Crossy12(EdgeIncrement)==0) && (Edges.Crossy21(EdgeIncrement)==0) )
                   newyflag=0;
               end
               
               fprintf(fid,'\t %d \t %d \n', newxflag, newyflag);

           end
       end

   end
end

cinside=cell(length(c),1);
CellIncrement=0;
for i=1:length(c)
   INGLOBAL=any((v(c{i},1)>0)&v(c{i},1)<LX&v(c{i},2)>0&v(c{i},2)<LY);
   if INGLOBAL
       CellIncrement=CellIncrement+1;
       cinside{CellIncrement}=EqVertex(c{i});

   end
end
NbCellsInside=CellIncrement;
cinside=cinside(1:NbCellsInside);
cinside2=cell(NbCellsInside,1);
CellIncrement=0;
for i=1:NbCellsInside
   dejavu=false;
   for j=1:i-1
       if ((length(cinside{i})==length(cinside{j}))&&all(sort(cinside{i})==sort(cinside{j}))) dejavu=true; end
   end
   if (~dejavu) CellIncrement=CellIncrement+1;cinside2{CellIncrement}=cinside{i}; end
end
cinside2=cinside2(1:CellIncrement);
fprintf(fid, 'Cells\n');
fprintf(fid, 'id\tA0\tEdgeNumb\tEdgeList\n');
for i=1:length(cinside2)
   ListVertex=c{i};
   fprintf(fid,'%d\t%f\t%d\t',i, A0, length(ListVertex));
   for j=1:length(ListVertex)
       NumeroEdge=0;
       % if we look at the current edge is in the list
       for k=1:EdgeIncrement
           if (EqVertex((ListVertex(j)))==Edges.Vertices(k,1)...
                   &&EqVertex(ListVertex(mod(j,length(ListVertex))+1))==Edges.Vertices(k,2))
               NumeroEdge=k;
           end % modulo is to be on that loop
           if(EqVertex(ListVertex(j))==Edges.Vertices(k,2)...
                   &&EqVertex(ListVertex(mod(j,length(ListVertex))+1))==Edges.Vertices(k,1))
               NumeroEdge=k; %MS change to always be positive
           end
       end
       % and is written in the file
       fprintf(fid,'%d,',NumeroEdge);
   end
   fprintf(fid, '\n');
end
fprintf(fid, 'End\n');
fclose(fid);
