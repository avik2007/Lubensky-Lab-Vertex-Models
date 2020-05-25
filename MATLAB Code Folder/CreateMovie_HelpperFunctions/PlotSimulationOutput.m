function PlotSimulationOutput(V,E,C,BoxSize,vertexMarker,isCellIncluded,xMax, yMax)
%If xMax and yMax are specified, set the maximum x and y axis to be xMax
%and yMax, respectively. Otherwise, the axes are determined from BoxSize.



%%
hold on;



%% Plot Cells (only if isCellIncluded = true)

if isCellIncluded
%if isCellIncluded = false, we won't try to draw the cells and this entire
%section WILL BE SKIPPED.


    %Find cells that do not cross any boundary
    cind = find(C.xflag==0 & C.yflag==0);

    %for every cells, go through its vertices in CWVertices and find the vertex
    %index
    for i = 1:length(cind)
        cwverts = C.CWVertices{cind(i)};
        indices = zeros(length(cwverts), 1);
        for j = 1:length(cwverts)
            indices(j) = find(V.id == cwverts(j));
        end

        %tag containing keyword 'Cell' followed by cell id
        TagStr = ['Cell ', num2str(C.id( cind(i)) )];

        %Determine color based on cell type
        ctype = C.Type(cind(i));
        switch ctype
            case 0
                Color = 'w';
            case 1
                Color = 'm';
            case 2
                Color = 'm';
            case 3
                Color = 'c';
            case 4
                Color = 'w';
            case 5 
                Color = 'k';
            case 6 
                Color = 'w';
            case 7 
                Color = 'y';
            otherwise
                error(['Cell Type ', ctype, ' is not associated with any color. Add the color for that cell type in PlotSimulationOutput.m'])    
        end


        fill(V.x(indices), V.y(indices), Color, 'Tag', TagStr)
    end



end %if isCellIncluded
    


%% Plot Edges

N = size(E.id,1);
for i = 1:N
    vind1 = find(V.id==E.vertex1(i));
    vind2 = find(V.id==E.vertex2(i));
    TagStr = ['Edge ', num2str(E.id(i))];
    line([V.x(vind1),V.x(vind2)], [V.y(vind1),V.y(vind2)], 'Color', 'k', 'LineWidth', 2.0, 'Tag', TagStr);
end




%% Plot Vertices

if vertexMarker == true
N = size(V.id, 1);
    for i = 1:N
        TagStr = ['Vertex ', num2str(V.id(i))];
        if V.id(i) > 0          %Plot 'real' vertices with circle
            plot(V.x(i),V.y(i), 'o','MarkerFaceColor','b', 'MarkerSize',1.0,'Tag', TagStr);
        else                    %Plot boundary vertices (with negative ID) without circle
            plot(V.x(i),V.y(i), 'Tag', TagStr);
        end
    end
end



%% Draw box at BoxSize
LX = BoxSize(1);
LY = BoxSize(2);
line([0, LX], [0, 0], 'HitTest', 'off', 'Color', 'k');
line([0, LX], [LY, LY], 'HitTest', 'off', 'Color', 'k');
line([0, 0], [0, LY], 'HitTest', 'off', 'Color', 'k');
line([LX, LX], [0, LY], 'HitTest', 'off', 'Color', 'k');



%% Set axis properties

%If xMax and yMax are specified, set the maximum x and y axis to be xMax
%and yMax, respectively. Otherwise the axes limits are taken from BoxSize
ha = gca;
if nargin == 8
    set(ha, 'XLim', [0 xMax], 'YLim', [0 yMax], 'DataAspectRatio', [1 1 1]);
else
    set(ha, 'XLim', [0 BoxSize(1)], 'YLim', [0 BoxSize(2)], 'DataAspectRatio', [1 1 1]);
end


hold off

