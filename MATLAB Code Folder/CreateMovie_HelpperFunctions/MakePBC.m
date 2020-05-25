function [V_pbc,E_pbc,C_pbc] = MakePBC(V,E,C,BoxSize,isCellIncluded)

%For edges that cross boundaries, do the following:
%   1. Create vertices at crossing points (use negative id to distinguish
%   them from 'real' vertices)
%   2. For ease of drawing, split the crossing edges into segments defined
%   by the 2 original vertices and the vertices created above. Give the
%   edge segments the same id as the original edge, then delete the
%   original edge.
%   3. Split cells that cross boundaries into multiple cells with the same
%   id (but bounded by different edges and vertices). (TO BE IMPLEMENTED)


V_pbc = V;
E_pbc = E;
C_pbc = C;

vcount = -1;     %counter for newly created vertex id


%% PBC for edges

%For each edges that cross a boundary, find 2 crossing points on the
%boundary and split the edge into 2 segments: vertex1 to the 1st crossing
%point, and the 2nd crossing point to vertex2.
%Add these edge segments to the list (with the same id, length, etc except
%for vertex and flag data), then delete the original crossing edge.

%Also, keep a 'cell' list of ordered vertices like:
%   NewEdgeList{} = [edgeId], [vertexId1], [vboundaryId1, vboundaryId2  , ..., vboundaryIdN]
%                   [edgeId], [vertexId2], [vboundaryIdN, vboundaryIdN-1, ..., vboundaryId1] 
%useful for doing PBC on cells later.


N = size(E.id, 1);
edelete = [];               %for keeping track crossing edges to be deleted later
NewEdgeList{1,3} = [];    %useful for doing PBC on cells later
count = 1;

for i = 1:N     %go through all edges
    vertex = [E.vertex1(i), E.vertex2(i)];
    ind1 = find(V.id == vertex(1));
    ind2 = find(V.id == vertex(2));
    x = [V.x(ind1), V.x(ind2)];
    y = [V.y(ind1), V.y(ind2)];
    
    if (E.xflag(i)==0 && E.yflag(i)==0)         %If edge does not cross any boundary
        continue
        
    elseif (E.xflag(i)==1 && E.yflag(i)==0)     %If edge crosses vertical boundary
        [xsort,isort] = sort(x, 'ascend');      %sort x ascending
        vstart = vertex(isort(2));
        vend = vertex(isort(1));              %starting and ending vertices
        slope = (y(isort(1))-y(isort(2))) / (xsort(1)+BoxSize(1)-xsort(2));
        xboundary = zeros(2,1);
        yboundary = zeros(2,1);
        xboundary(1) = BoxSize(1);
        yboundary(1) = y(isort(2)) + (BoxSize(1) - xsort(2))*slope;  %point right before crossing boundary
        xboundary(2) = 0;
        yboundary(2) = yboundary(1);                                 %point right after crossing boundary
        
    elseif (E.xflag(i)==0 && E.yflag(i)==1)      %If edge crosses horizontal boundary
        [ysort,isort] = sort(y, 'ascend');       %sort y ascending
        vstart = vertex(isort(2));
        vend = vertex(isort(1));              %starting and ending vertices
        slope = (x(isort(1))-x(isort(2))) / (ysort(1)+BoxSize(2)-ysort(2));
        xboundary = zeros(2,1);
        yboundary = zeros(2,1);
        xboundary(1)= x(isort(2)) + (BoxSize(2) - ysort(2))*slope;
        yboundary(1) = BoxSize(2);                                   %point right before crossing boundary
        xboundary(2) = xboundary(1);
        yboundary(2) = 0;                                            %point right after crossing boundary
        
     elseif (E.xflag(i)==1 && E.yflag(i)==1)      %If edge crosses both boundaries
        [xsort,isortx] = sort(x, 'ascend');
        [ysort,isorty] = sort(y, 'ascend');
        xboundary = zeros(4,1);
        yboundary = zeros(4,1);
        
        if (isortx(1)==isorty(1))
            %crossing from top right to bottom left
            vstart = vertex(isortx(2));
            vend = vertex(isortx(1));              %starting and ending vertices
            slope = ((ysort(1)+BoxSize(2))-ysort(2)) / ((xsort(1)+BoxSize(1))-xsort(2));
            ycross = ysort(2) + (BoxSize(1)-xsort(2))*slope;   %y when crossing vertical boundary
            xcross = xsort(2) + (BoxSize(2)-ysort(2))/slope;   %x when crossing horizontal boundary
            
            if (xcross<=BoxSize(1) && ycross>=BoxSize(2))
                %crossing horizontal axis first, then vertical axis (or
                %crossing through a corner)
                xboundary(1) = xcross;
                yboundary(1) = BoxSize(2);                        %point right before horizontal axis crossing
                xboundary(2) = xboundary(1);
                yboundary(2) = 0;                                 %point right after horizontal axis crossing
                xboundary(3) = BoxSize(1);
                yboundary(3) = ycross - BoxSize(2);               %point right before vertical axis crossing
                xboundary(4) = 0;
                yboundary(4) = yboundary(3);                      %point right after vertical axis crossing
            elseif (xcross>BoxSize(1) && ycross<BoxSize(2))
                %crossing vertical axis first, then horizontal axis
                xboundary(1) = BoxSize(1);
                yboundary(1) = ycross;                       %point right before vertical axis crossing
                xboundary(2) = 0;
                yboundary(2) = yboundary(1);                   %point right after vertical axis crossing
                xboundary(3) = xcross - BoxSize(1);
                yboundary(3) = BoxSize(2);                   %point right before horizontal axis crossing
                xboundary(4) = xboundary(3);
                yboundary(4) = 0;                            %point right after horizontal axis crossing
            else
                error('??!');
            end
            
            
        elseif (isortx(1)==isorty(2))
            %crossing from bottom right to top left
            vstart = vertex(isortx(2));
            vend = vertex(isortx(1));              %starting and ending vertices
            slope = ((ysort(2)-BoxSize(2))-ysort(1)) / ((xsort(1)+BoxSize(1))-xsort(2));
            ycross = ysort(1) + (BoxSize(1)-xsort(2))*slope;   %y when crossing vertical boundary
            xcross = xsort(2) - ysort(1)/slope;                %x when crossing horizontal boundary
            
            if (xcross<BoxSize(1) && ycross<0)
                %crossing horizontal axis first, then vertical axis
                xboundary(1) = xcross;
                yboundary(1) = 0;                                 %point right before horizontal axis crossing
                xboundary(2) = xboundary(1);
                yboundary(2) = BoxSize(2);                        %point right after horizontal axis crossing
                xboundary(3) = BoxSize(1);
                yboundary(3) = BoxSize(2) + ycross;               %point right before vertical axis crossing
                xboundary(4) = 0;
                yboundary(4) = yboundary(3);                      %point right after vertical axis crossing
            elseif (xcross>BoxSize(1) && ycross>0)
                %crossing vertical axis first, then horizontal axis
                xboundary(1) = BoxSize(1);
                yboundary(1) = ycross;                            %point right before vertical axis crossing
                xboundary(2) = 0;
                yboundary(2) = yboundary(1);                        %point right after vertical axis crossing
                xboundary(3) = xcross - BoxSize(1);
                yboundary(3) = 0;                                 %point right before horizontal axis crossing
                xboundary(4) = xboundary(3);
                yboundary(4) = BoxSize(2);                        %point right after horizontal axis crossing
            else
                error('??');
            end
            
        else
            error('?!?');
        end
        
    else
        error('Flag must be either 0 or 1!');
    end
    
    
    orderedBoundaryVerts = [];
    
    %Add boundary points to list of vertices
    n = length(xboundary);
    vnewidTemp = zeros(n,1);        %list of ID of the newly-created boundary vertices
    for j = 1:n
        V_pbc.id(end+1) = vcount;
        V_pbc.x(end+1) = xboundary(j);
        V_pbc.y(end+1) = yboundary(j);
        V_pbc.EdgesNumber(end+1) = NaN;
        vnewidTemp(j) = vcount;
        vcount = vcount - 1;        %give boundary vertices negative IDs
    end
    
    %Add edge segments to list of edges
    vnewid = [vstart; vnewidTemp; vend];
    for j = 1:2:length(vnewid)-1        %stepsize is 2
        E_pbc.id(end+1) = E.id(i);
        E_pbc.vertex1(end+1) = vnewid(j);
        E_pbc.vertex2(end+1) = vnewid(j+1);
        E_pbc.xflag(end+1) = 0;
        E_pbc.yflag(end+1) = 0;
        E_pbc.length(end+1) = E.length(i);
        E_pbc.cell1(end+1) = E.cell1(i);
        E_pbc.cell2(end+1) = E.cell2(i);
    end
    
    %Record index for deleting original (crossing) edges later
    edelete = [edelete, i];
    
    
    %Create entry of crossedgelist
    NewEdgeList{count,1} = E.id(i);
    NewEdgeList{count,2} = vstart;
    NewEdgeList{count,3} = vnewidTemp;
    
    NewEdgeList{count+1, 1} = E.id(i);
    NewEdgeList{count+1, 2} = vend;
    NewEdgeList{count+1, 3} = vnewidTemp(end:-1:1);   %flip the order
    
    count = count + 2;
        
end


%Remove the original edge from the list
E_pbc.id(edelete) = [];
E_pbc.vertex1(edelete) = [];
E_pbc.vertex2(edelete) = [];
E_pbc.xflag(edelete) = [];
E_pbc.yflag(edelete) = [];
E_pbc.length(edelete) = [];
E_pbc.cell1(edelete) = [];
E_pbc.cell2(edelete) = [];



%Add field E.Midpoint specifying coordinates of the midpoint of an edge.
%Useful for being able to click on an edge and display properties
N = length(E_pbc.id);
E_pbc(1).midpoint = zeros(N,2);
for i = 1:N
    vind1 = find(V_pbc.id == E_pbc.vertex1(i));
    vind2 = find(V_pbc.id == E_pbc.vertex2(i));
    xmid = (V_pbc.x(vind1) + V_pbc.x(vind2))/2;
    ymid = (V_pbc.y(vind1) + V_pbc.y(vind2))/2;
    E_pbc.midpoint(i,:) = [xmid, ymid];
end


    
    
%% PBC for cells (only if isCellIncluded = true)

%if isCellIncluded = false, we won't try to draw the cells and this entire
%section WILL BE SKIPPED.
if isCellIncluded


    % For every cells that cross some boundary, break the cells into multiple
    % connected segments by adding boundary points (i.e. vertex with negative
    % id created above when splitting crossing edges)

    % Outline: create 4 empty lists, corresponding to cell segments in each
    % quadrant. (Some lists may be empty, e.g. if the cell only cross one
    % boundary)
    %   NewCell{i,1} = list of vertices where no boundary has been crossed from the starting point (flag = [0 0])
    %   NewCell{i,2} = list of vertices where only xboundary has been crossed (flag = [1 0])
    %   NewCell{i,3} = list of vertices where only yboundary has been crossed (flag = [0 1])
    %   NewCell{i,4} = list of vertices where both boundaries has been crossed (flag = [1 1])
    % Walk through CW-ordered list of vertices, and put each vertex in the
    % appropriate list. Finally, if the cell crosses both boundaries, determine
    % from the order of boundaries crossed whether we need to include the
    % corner points.


    % put id and vertex1 from NewEdgeList into a vector for ease of use
    NewEdgeId = [NewEdgeList{:,1}];
    NewEdgeV1 = [NewEdgeList{:,2}];


    % indices of cells that cross one or more boundary
    cind = find(C.xflag==1 | C.yflag==1);
    N = length(cind);
    NewCell = cell(N, 4);

    for i = 1:N
        cwverts = C.CWVertices{cind(i)};
        cwverts(end+1) = cwverts(1);            %make the list 'circular'

        currentList = 1;
        NewCell{i,currentList} = cwverts(1);              %put 1st vertex in 1st list

        for j = 1:length(cwverts)-1
            v1 = cwverts(j);
            v2 = cwverts(j+1);
            %find index of edge between v1 and v2 (be careful to use the 
            %non-pbc list E, instead of Epbc)
            eindex = find( (E.vertex1==v1 & E.vertex2==v2) | (E.vertex1==v2 & E.vertex2==v1) );
            eid = E.id(eindex);
            xflag = E.xflag(eindex);
            yflag = E.yflag(eindex);


            if (xflag==0 && yflag==0);
                %if no boundary is crossed, put v2 in list and continue
                NewCell{i,currentList} = [NewCell{i,currentList}, v2];
                continue;
            elseif (xflag==1 && yflag==0)
                %if xboundary is crossed: find vboundary1 & vboundary2, put
                %vboundary1 into current list, switch to list #2, put
                %vboundary2 and v2 into that list.
                index = find(NewEdgeId==eid & NewEdgeV1==v1, 1);
                vboundary1 = NewEdgeList{index,3}(1);
                vboundary2 = NewEdgeList{index,3}(2);

                NewCell{i,currentList} = [NewCell{i,currentList}, vboundary1];

                if currentList == 1
                    currentList = 2;
                elseif currentList == 2
                    currentList = 1;
                elseif currentList == 3
                    currentList = 4;
                elseif currentList == 4
                    currentList = 3;
                end

                NewCell{i,currentList} = [NewCell{i,currentList}, vboundary2, v2];
                continue;

            elseif (xflag==0 && yflag==1)
                %if xboundary is crossed: find vboundary1 & vboundary2, put
                %vboundary1 into current list, switch to list #3, put
                %vboundary2 and v2 into that list.
                index = find(NewEdgeId==eid & NewEdgeV1==v1, 1);
                vboundary1 = NewEdgeList{index,3}(1);
                vboundary2 = NewEdgeList{index,3}(2);

                NewCell{i,currentList} = [NewCell{i,currentList}, vboundary1];

                if currentList == 1
                    currentList = 3;
                elseif currentList == 2
                    currentList = 4;
                elseif currentList == 3
                    currentList = 1;
                elseif currentList == 4
                    currentList = 2;
                end

                NewCell{i,currentList} = [NewCell{i,currentList}, vboundary2, v2];
                continue;

            elseif (xflag==1 && yflag==1)
                %if both boundary is crossed, find out which boundary is first
                %crossed. Then do something similar to the above.
                index = find(NewEdgeId==eid & NewEdgeV1==v1, 1);
                vboundary1 = NewEdgeList{index,3}(1);
                vboundary2 = NewEdgeList{index,3}(2);
                vboundary3 = NewEdgeList{index,3}(3);
                vboundary4 = NewEdgeList{index,3}(4);

                %get x position of vboundary1 and vboundary2
                indboundary1 = (V_pbc.id == vboundary1);
                indboundary2 = (V_pbc.id == vboundary2);
                xboundary1 = V_pbc.x(indboundary1);
                xboundary2 = V_pbc.x(indboundary2);

                NewCell{i,currentList} = [NewCell{i,currentList}, vboundary1];

                if (xboundary1==0) || (xboundary2==0)
                    %if xboundary is crossed first
                    if currentList == 1             %cross xboundary
                        currentList = 2;
                    elseif currentList == 2
                        currentList = 1;
                    elseif currentList == 3
                        currentList = 4;
                    elseif currentList == 4
                        currentList = 3;
                    end
                    NewCell{i,currentList} = [NewCell{i,currentList}, vboundary2, vboundary3];
                    if currentList == 1             %cross yboundary
                        currentList = 3;
                    elseif currentList == 2
                        currentList = 4;
                    elseif currentList == 3
                        currentList = 1;
                    elseif currentList == 4
                        currentList = 2;
                    end
                    NewCell{i,currentList} = [NewCell{i,currentList}, vboundary4, v2];

                else
                    %if yboundary is crossed first
                    if currentList == 1             %cross yboundary
                        currentList = 3;
                    elseif currentList == 2
                        currentList = 4;
                    elseif currentList == 3
                        currentList = 1;
                    elseif currentList == 4
                        currentList = 2;
                    end
                    NewCell{i,currentList} = [NewCell{i,currentList}, vboundary2, vboundary3];
                    if currentList == 1             %cross xboundary
                        currentList = 2;
                    elseif currentList == 2
                        currentList = 1;
                    elseif currentList == 3
                        currentList = 4;
                    elseif currentList == 4
                        currentList = 3;
                    end
                    NewCell{i,currentList} = [NewCell{i,currentList}, vboundary4, v2];
                end
            end

        end


        %for cell that has segments in all four quadrants, we need to include
        %the corner points
        if ( ~isempty(NewCell{i,1}) && ~isempty(NewCell{i,2}) ...
                && ~isempty(NewCell{i,3}) && ~isempty(NewCell{i,4}) )
            %create the four corner points
            minId = min(V_pbc.id);
            V_pbc.id(end+1) = minId - 1;        %bottom-left
            V_pbc.x(end+1) = 0;
            V_pbc.y(end+1) = 0;
            V_pbc.id(end+1) = minId - 2;        %bottom-right
            V_pbc.x(end+1) = BoxSize(1);
            V_pbc.y(end+1) = 0;
            V_pbc.id(end+1) = minId - 3;        %top-left
            V_pbc.x(end+1) = 0;
            V_pbc.y(end+1) = BoxSize(2);
            V_pbc.id(end+1) = minId - 4;        %top-right
            V_pbc.x(end+1) = BoxSize(1);
            V_pbc.y(end+1) = BoxSize(2);

            for k = 1:4
                %loop over all segments to insert the appropriate corner points
                temp = NewCell{i,k};
                pos = find(temp < 0);  %find boundary vertices
                boundaryverts = temp(pos);
                %get position of the boundary vertices
                ind1 = find(V_pbc.id == boundaryverts(1), 1);
                ind2 = find(V_pbc.id == boundaryverts(2), 1);
                x1 = V_pbc.x(ind1);
                y1 = V_pbc.y(ind1);
                x2 = V_pbc.x(ind2);
                y2 = V_pbc.y(ind2);

                if ( (x1==0 || x2==0) && (y1==0 || y2==0) )
                    %if the cell segment contains the bottom-left corner
                    CornerVertId = minId - 1;
                elseif (x1==BoxSize(1) || x2==BoxSize(1)) && (y1==0 || y2==0)
                    %if the cell segment contains the bottom-right corner
                    CornerVertId = minId - 2;
                elseif (x1==0 || x2==0) && (y1==BoxSize(2) || y2==BoxSize(2))
                    %if the cell segment contains the top-left corner
                    CornerVertId = minId - 3;
                elseif (x1==BoxSize(1) || x2==BoxSize(1)) && (y1==BoxSize(2) || y2==BoxSize(2))
                    %if the cell segment contains the top-right corner
                    CornerVertId = minId - 4;
                else
                    error('err..?')
                end

                %now insert the corner point between the two boundary vertices
                %(which have negative id's)
                if pos(1) == 1
                    NewCell{i,k} = [CornerVertId, temp];
                else
                    NewCell{i,k} = [temp(1:pos(1)), CornerVertId, temp(pos(2):end)];
                end


            end
        end
    end

    Nold = size(C.id,1);
    count = 0;

    %Add 'new' cell segments into Cpbc
    for i = 1:N
        for j = 1:4
            if isempty(NewCell{i,j})
                continue;
            end
            count = count + 1;

            %copy values from original cell
            C_pbc.id(Nold+count) = C.id(cind(i));
            C_pbc.EdgeNumb(Nold+count) = C.EdgeNumb(cind(i));
            C_pbc.EdgeList(Nold+count) = C.EdgeList(cind(i));
            C_pbc.Type(Nold+count) = C.Type(cind(i));
            C_pbc.A0(Nold+count) = C.A0(cind(i));
            C_pbc.Area(Nold+count) = C.Area(cind(i));
            C_pbc.Pressure(Nold+count) = C.Pressure(cind(i));
            C_pbc.Center_x(Nold+count) = C.Center_x(cind(i));
            C_pbc.Center_y(Nold+count) = C.Center_y(cind(i));


            %new values for cell segment
            C_pbc.CWVertices(Nold+count) = NewCell(i,j);
            C_pbc.xflag(Nold+count) = 0;
            C_pbc.yflag(Nold+count) = 0;

        end
    end


    %delete original cells
    C_pbc.id(cind) = [];
    C_pbc.EdgeNumb(cind) = [];
    C_pbc.EdgeList(cind) = [];
    C_pbc.Type(cind) = [];
    C_pbc.CWVertices(cind) = [];
    C_pbc.A0(cind) = [];
    C_pbc.Area(cind) = [];
    C_pbc.Pressure(cind) = [];
    C_pbc.xflag(cind) = [];
    C_pbc.yflag(cind) = [];
    C_pbc.Center_x(cind) = [];
    C_pbc.Center_y(cind) = [];



end %if isCellIncluded







