function [Points1] = GenerationVoronoi_FloydIteration()


r = 504; %seed
NbPoints = 2000; %number of cells 
limit = 3; %number of iterations


%% Generation of a series of points
close all;

%% Parameters
rng(r);
Mesh = 250 * sqrt(2) * sqrt(NbPoints / 100);
A0 = Mesh * Mesh / NbPoints;

%Global parameters
SEED = 48269; %originally 48263
TIMEINTERVAL =	1.0;
MINAREA = 0.1;
MINLENGTH =	1.0;
NEWEDGELENGTH =	2.0;
K =	0.003264;  
GAMMA =	17.4859;
TOTALTIME = 500;
CURRENTTIME = 0.0;
INITIALSTEPSIZE = .0001;
error_tolerance_absolute_mechanical =	0.001;
TISSUE_ALPHA = 100;
TISSUE_EXT_STRESS_X = 0.4946;
TISSUE_EXT_STRESS_Y = 0.4946;

%MS: Output should be labled as 'VoronoiCells NbC (number of cells) RNG
%(the rng number).txt'
nomfichier=strcat('NbC_', num2str(NbPoints,'%04d'), '_rng_', num2str(r,'%4d'),  '.txt');

%%
f1=figure();
f2=figure();
Points1=Mesh*rand(NbPoints,2); %Matrix of x,y locatios of cell centers
% copy on neighboring points
% take to points you have and coppy them into the top bottom sides and
% corrner reigons sorounding the 400 by 400 middle protion thus you end
% with a huge space of 1200 by 1200  that is repeated 9 times in squares of
% 400
it = 0;
while it < limit
    Pointsa=Points1+Mesh*[ones(NbPoints,1) zeros(NbPoints,1)];
    Pointsb=Points1-Mesh*[zeros(NbPoints,1) ones(NbPoints,1)];
    Pointsc=Points1-Mesh*[ones(NbPoints,1) zeros(NbPoints,1)];
    Pointsd=Points1+Mesh*[zeros(NbPoints,1) ones(NbPoints,1)];
    Pointse=Points1+Mesh*[ones(NbPoints,1) ones(NbPoints,1)];
    Pointsf=Points1+Mesh*[ones(NbPoints,1) -ones(NbPoints,1)];
    Pointsg=Points1+Mesh*[-ones(NbPoints,1) ones(NbPoints,1)];
    Pointsh=Points1+Mesh*[-ones(NbPoints,1) -ones(NbPoints,1)];
    Pointst=[Points1; Pointsa; Pointsb; Pointsc; Pointsd; Pointse; Pointsf; Pointsg; Pointsh];
    
    
    %% Calculate the Voronoi diagram
    % c contains a list of each cell
    % v(c{5},1) and v(c{5},2) contains a list of the x and y positions of each
    % vertex sourounding a cell
    [v,c]=voronoin(Pointst);
    
    centroids = zeros(NbPoints,2);
    centroidcount = 0;
    for i = 1:length(c)
        if (all(c{i}~=1))
            centroidx = mean(v(c{i},1));
            centroidy = mean(v(c{i},2));
            if centroidx > 0 && centroidx < Mesh && centroidy > 0 && centroidy < Mesh
                centroidcount = centroidcount + 1;
                centroids(centroidcount,:) = [centroidx centroidy];
            end
        end
    end
    
    it = it + 1;
    Points1 = centroids;
end

 %% Representation of the result
    set(0,'CurrentFigure',f1);
    plot(Points1(:,1),Points1(:,2),'r+');
    set(0,'CurrentFigure',f2);
    plot(Points1(:,1),Points1(:,2),'r+');
    hold on
    for i = 1:length(c)
        
        % all traceable cells
        if (all(c{i}~=1))
            set(0,'CurrentFigure',f1)
            patch(v(c{i},1),v(c{i},2),i)
        end
        %ca are the cells that are in the interior only
        IN=(all(v(c{i},1)>0)&all(v(c{i},1)<Mesh)&all(v(c{i},2)>0)&all(v(c{i},2)<Mesh));
        %ca have cells that are at the border
        IN2=(any(v(c{i},1)>0)&any(v(c{i},1)<Mesh)&any(v(c{i},2)>0)&any(v(c{i},2)<Mesh))&(~IN);
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
    xlim([0 Mesh]);
    ylim([0 Mesh]);
    plot(Points1(:,1),Points1(:,2),'bo');
    scatter(centroids(:,1),centroids(:,2),'r+')


for i = 1:length(c)
    INGLOBAL=any((v(c{i},1)>0)&v(c{i},1)<Mesh&v(c{i},2)>0&v(c{i},2)<Mesh);
    
    if INGLOBAL
        ListVertex=c{i};
        area=0;
        for k=1:length(c{i})
            kplus=mod(k,length(c{i}))+1;
            area=area+v(c{i}(k),1)*v(c{i}(kplus),2)-v(c{i}(kplus),1)*v(c{i}(k),2);
        end
        if(area<0)    c{i}=fliplr(c{i}); end
    end
end

%% Export to a file
fid=fopen(nomfichier,'wt');
fprintf(fid, 'Voronoi Tesselation Data: \n');
fprintf(fid, 'NumbCells %d \t rng %d \n', NbPoints, r);
fprintf(fid, 'Global:\n');
fprintf(fid, 'LX\t%f\nLY\t%f\n', Mesh, Mesh);
fprintf(fid, 'A0\t%f\n', A0);
fprintf(fid, 'SEED\t%d\n', SEED);
fprintf(fid, 'TIMEINTERVAL\t%f\n', TIMEINTERVAL);
fprintf(fid, 'MINAREA\t%f\n', MINAREA);
fprintf(fid, 'MINLENGTH\t%f\n', MINLENGTH);
fprintf(fid, 'NEWEDGELENGTH\t%f\n', NEWEDGELENGTH);
fprintf(fid, 'K\t%f\n', K);
fprintf(fid, 'GAMMA\t%f\n', GAMMA);
fprintf(fid, 'TOTALTIME\t%f\n', TOTALTIME);
fprintf(fid, 'CURRENTTIME\t%f\n', CURRENTTIME);
fprintf(fid, 'INITIALSTEPSIZE\t%f\n', INITIALSTEPSIZE);
fprintf(fid, 'error_tolerance_absolute_mechanical\t%f\n', error_tolerance_absolute_mechanical);
fprintf(fid, 'TISSUE_ALPHA\t%f\n', TISSUE_ALPHA);
fprintf(fid, 'TISSUE_EXT_STRESS_X\t%f\n', TISSUE_EXT_STRESS_X);
fprintf(fid, 'TISSUE_EXT_STRESS_Y\t%f\n', TISSUE_EXT_STRESS_Y);
fprintf(fid, 'endGlobal\n');
fprintf(fid, 'Vertices\n');
fprintf(fid, 'id\tx\ty\n');
% all those in the interior must be included
NbVertex=0;
CentralVertices=zeros(length(v),2);
for i = 2:length(v)   % v(1) is always a point at infinity, so skip (NOTED ADDED BY DKL 7/5/10)
    if(v(i,1)>0&&v(i,1)<Mesh&&v(i,2)>0&&v(i,2)<Mesh)
        fprintf(fid,'%d',i);
        fprintf(fid,'\t %2.4f',v(i,1));
        fprintf(fid,'\t %2.4f',v(i,2));
        fprintf(fid,'\n');
        NbVertex=NbVertex+1;
        CentralVertices(NbVertex,1)=v(i,1);
        CentralVertices(NbVertex,2)=v(i,2);
    end
end

% keep only the interesting portion
CentralVertices=CentralVertices(1:NbVertex,1:2);
% I want to create a file of equivalence, which gives for each vertex
% its equivalent in the main mesh (ca will be useful later)
EqVertex=zeros(NbVertex,1);
% the first is less inf in any case

for i = 2:length(v)
    vcorrige=zeros(1,2);
    % vcorrige find the vertex position equivalent
    vtest=v(i,:);
    
    if (vtest(1)>0&&vtest(1)<Mesh&&vtest(2)>0&&vtest(2)<Mesh) vcorrige=vtest;end
    vtest=v(i,:)-Mesh*[1 0];
    if (vtest(1)>0&&vtest(1)<Mesh&&vtest(2)>0&&vtest(2)<Mesh) vcorrige=vtest;  end
    vtest=v(i,:)-Mesh*[0 1];
    if (vtest(1)>0&&vtest(1)<Mesh&&vtest(2)>0&&vtest(2)<Mesh) vcorrige=vtest;  end
    vtest=v(i,:)+Mesh*[1 0];
    if (vtest(1)>0&&vtest(1)<Mesh&&vtest(2)>0&&vtest(2)<Mesh) vcorrige=vtest;  end
    vtest=v(i,:)+Mesh*[0 1];
    if (vtest(1)>0&&vtest(1)<Mesh&&vtest(2)>0&&vtest(2)<Mesh) vcorrige=vtest;  end
    vtest=v(i,:)-Mesh*[1 1];
    if (vtest(1)>0&&vtest(1)<Mesh&&vtest(2)>0&&vtest(2)<Mesh) vcorrige=vtest;  end
    vtest=v(i,:)-Mesh*[1 -1];
    if (vtest(1)>0&&vtest(1)<Mesh&&vtest(2)>0&&vtest(2)<Mesh) vcorrige=vtest;  end
    vtest=v(i,:)+Mesh*[1 1];
    if (vtest(1)>0&&vtest(1)<Mesh&&vtest(2)>0&&vtest(2)<Mesh) vcorrige=vtest;  end
    vtest=v(i,:)+Mesh*[1 -1];
    if (vtest(1)>0&&vtest(1)<Mesh&&vtest(2)>0&&vtest(2)<Mesh) vcorrige=vtest;  end
    % should give the issue of what has been found
    column=find(all(abs((repmat(vcorrige,size(v,1),1))'-v')<0.000001));%strmatch(vcorrige, v);
    if(~isempty(column))
        EqVertex(i)=column;
    end
    
end

% then you must create the list of Edges
fprintf(fid, 'Edges\n');
fprintf(fid, 'id\tvertex1\tvertex2\txflag\tyflag\n');
NbEdgeMax=3*NbVertex; %NbVertex=number of verticies in center 400 by 400 grid
% contains the number of both sides
Edges.Vertices=zeros(NbEdgeMax,2);
Edges.Crossx12=zeros(NbEdgeMax,1);
Edges.Crossx21=zeros(NbEdgeMax,1);
Edges.Crossy12=zeros(NbEdgeMax,1);
Edges.Crossy21=zeros(NbEdgeMax,1);

EdgeIncrement=0; %the number of the edge so you get 1,2,3,4 as edge lables
for i = 1:length(c) %For each vertex
    %IN=(all(v(c{i},1)>0)&all(v(c{i},1)<1)&all(v(c{i},2)>0)&all(v(c{i},2)<1));
    %INGLOBAL=(any(v(c{i},1)>0)&any(v(c{i},1)<1)&any(v(c{i},2)>0)&any(v(c{i},2)<1));
    % condition INGLOBAL takes only cells that have at least one vertex
    % in the central square
    INGLOBAL=any( (v(c{i},1)>0) & v(c{i},1)<Mesh & v(c{i},2)>0 & v(c{i},2)<Mesh);
    % if cells that are at least partially in the interior
    
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
                    if(v(ListVertex(j),1)>Mesh) Edges.Crossx21(EdgeIncrement)=1;%fprintf(fid,' xCrossBoundary21');
                        %   x21=1;
                    end
                    if(v(ListVertex(j),2)<0) Edges.Crossy12(EdgeIncrement)=1;%fprintf(fid,' yCrossBoundary12');
                        %   y12=1;
                    end
                    if(v(ListVertex(j),2)>Mesh) Edges.Crossy21(EdgeIncrement)=1;%fprintf(fid,' yCrossBoundary21');
                        %    y21=1;
                    end
                end
                if((EqVertex(ListVertex(j))==ListVertex(j))&&(EqVertex(ListVertex(mod(j,length(ListVertex))+1))~=ListVertex(mod(j,length(ListVertex))+1)))
                    if(v(ListVertex(mod(j,length(ListVertex))+1),1)<0) Edges.Crossx21(EdgeIncrement)=1;%fprintf(fid,' xCrossBoundary21');
                        % x21=1;
                    end
                    if(v(ListVertex(mod(j,length(ListVertex))+1),1)>Mesh) Edges.Crossx12(EdgeIncrement)=1;%fprintf(fid,' xCrossBoundary12');
                        % x12=1;
                    end
                    if(v(ListVertex(mod(j,length(ListVertex))+1),2)<0) Edges.Crossy21(EdgeIncrement)=1;%fprintf(fid,' yCrossBoundary21');
                        % y21=1;
                    end
                    if(v(ListVertex(mod(j,length(ListVertex))+1),2)>Mesh) Edges.Crossy12(EdgeIncrement)=1;%fprintf(fid,' yCrossBoundary12');
                        % y12=1;
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
                    if(v2(1)>v1(1))&&(eqv2(1)<eqv1(1)) Edges.Crossx12(EdgeIncrement)=1;%fprintf(fid,' xCrossBoundary12');
                    end
                    if(v2(1)<v1(1))&&(eqv2(1)>eqv1(1)) Edges.Crossx21(EdgeIncrement)=1;%fprintf(fid,' xCrossBoundary21');
                    end
                    if (v2(2)>v1(2))&&(eqv2(2)<eqv1(2))Edges.Crossy12(EdgeIncrement)=1;%fprintf(fid,' yCrossBoundary12');
                    end
                    if (v2(2)<v1(2))&&(eqv2(2)>eqv1(2)) Edges.Crossy21(EdgeIncrement)=1;%fprintf(fid,' yCrossBoundary21');
                    end
                end
                %  fprintf(fid,'\n');
                %         if(Edges.Crossx12(EdgeIncrement)==1||Edges.Crossx21(EdgeIncrement)==1)
                %            fprintf(fid,' 1 ');
                %       else fprintf(fid,' 0 ');
                %      end
                %     if(Edges.Crossy12(EdgeIncrement)==1||Edges.Crossy21(EdgeIncrement)==1)
                %        fprintf(fid,' 1\n');
                %   else fprintf(fid,' 0\n');
                %  end
                
                
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

%% Then we go to the cells
% I want to create a shortened list of cells that point to the
% of green inside
cinside=cell(length(c),1);
CellIncrement=0;
for i=1:length(c)
    INGLOBAL=any((v(c{i},1)>0)&v(c{i},1)<Mesh&v(c{i},2)>0&v(c{i},2)<Mesh);
    if INGLOBAL
        CellIncrement=CellIncrement+1;
        cinside{CellIncrement}=EqVertex(c{i});
    end
end
NbCellsInside=CellIncrement;
% to remove the useful lines
cinside=cinside(1:NbCellsInside);
cinside2=cell(NbCellsInside,1);
% loop to remove repetitions
CellIncrement=0;
for i=1:NbCellsInside
    dejavu=false;
    for j=1:i-1
        if ((length(cinside{i})==length(cinside{j}))&&all(sort(cinside{i})==sort(cinside{j}))) dejavu=true; end
    end
    if (~dejavu) CellIncrement=CellIncrement+1;cinside2{CellIncrement}=cinside{i}; end
end
cinside2=cinside2(1:CellIncrement);
% from this list we can write cells
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


