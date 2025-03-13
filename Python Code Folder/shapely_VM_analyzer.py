import os
import pandas as pd
import numpy as np
from scipy.spatial import Voronoi, voronoi_plot_2d, ConvexHull, Delaunay
import matplotlib.pyplot as plt
import seaborn as sns
from PIL import Image
from process_voronoi import MyVoronoi, getAreas, PolyArea, transform_array_to_coords, length, doIntersect, point_hash, hash_based_edge_check
import pickle
from scipy.ndimage import gaussian_filter, uniform_filter, uniform_filter1d,laplace, convolve1d, convolve, binary_dilation, zoom
import sys
sys.path.insert(0, '../raw_data_and_preprocessing')
from preprocessing import get_fly_index
import xarray as xr
import geopandas as gpd
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['savefig.facecolor'] = 'white'
import image_smoothing_toolbox
import timeit
from skimage import measure
from shapely import Polygon, MultiPolygon, LineString
#from shapely_fly_analyzer import maskChecker
import network_polygons as tp
import geopandas as gpd
import matplotlib.cm as cm
#import shapely_fly_analyzer as SFA
#from shapely_fly_analyzer import maskChecker
#from shapely_fly_analyzer import maskChecker, GetPandasDFOfBristleCoordinates, maskChecker
#from network_polygons import getBoundaryCorrectedGraphLaplacianFromA, getLaplacianMatrixFromA
#import wing_ensemble_analyzer as WEA
import image_smoothing_toolbox as ist
#from matplotlib_scalebar.scalebar import ScaleBar
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.patches as mpatches
from skimage.draw import polygon

class ShapelyVMAnalyzer():
    def __init__(self, file = None):
        self.file = file
        self.SVM = None
        self.celldf = None
        self.edgedf = None
        self.vertexdf = None
        self.Lx = 0
        self.Ly = 0
        self.isValid = False;
        self.P0 = 0;
        self.K = 0;
        self.Lambda = 0;

    def readFileMakeVars_HOMOGENEOUS(self):
        # read the file 
        filetext= open(self.file, 'r')
        content = filetext.read()
        filetext.close()

        Lx_index = content.find("LX\t")
        #print(Lx_index)
        Lx_endindex = content[3+Lx_index:].find("\n")

        Ly_index = content.find("LY\t")
        Ly_endindex = content[3+Ly_index:].find("\n")
        self.Lx = float(content[3+Lx_index:3+ Lx_index + Lx_endindex])
        self.Ly = float(content[3+Ly_index:3+ Ly_index + Ly_endindex])
        P0_index = content.find("P0\t")
        P0_endindex = content[3+P0_index:].find("\n")
        self.P0 = float(content[3+P0_index:3+P0_index+P0_endindex])
        K_index = content.find("K\t")
        K_endindex = content[2+K_index:].find("\n")
        self.K = float(content[2+K_index:2+K_index+K_endindex])
        lambda_index = content.find("LAMBDA\t")
        lambda_endindex = content[len('LAMBDA\t')+lambda_index:].find("\n")
        self.Lambda = float(content[len('LAMBDA\t') + lambda_index: len('LAMBDA\t') + lambda_index + lambda_endindex])

        
        # find where the various sections go
        vertex_index = content.find("Vertices")
        edge_index = content.find("Edges\n")
        cell_index = content.find("Cells")
        end_index = content.find("End\n")

        # read the vertex bits
        vertex_substr = content[vertex_index:edge_index]
        vertex_list = vertex_substr.split("\n")
        #vertex_list[2:] # list of vertex data
        vert_var_list = np.array(vertex_list[1].split("\t"))
        vertex_data_array = np.array([subl.split("\t") for subl in vertex_list[2:-1]])

        # read the edge bits
        edge_substr = content[edge_index:cell_index]
        edge_list = edge_substr.split("\n")
        edge_var_list = np.array(edge_list[1].split("\t"))
        edge_data_array = np.array([subl.split("\t") for subl in edge_list[2:-1]])

        # read the cell bits
        cell_substr = content[cell_index:end_index]
        cell_list = cell_substr.split("\n")
        cell_var_list = np.array(cell_list[1].split("\t"))
        cell_data_array = np.array([subl.split("\t") for subl in cell_list[2:-1]])

        # make pandas arrays 
        vertindexarray = [int(v) for v in vertex_data_array[:, 0]]
        edgeindexarray = [int(e) for e in edge_data_array[:,0]]
        cellindexarray = [int(c) for c in cell_data_array[:,0]]
        self.vertexdf = pd.DataFrame(columns = vert_var_list[1:], index = vertindexarray, data = vertex_data_array[:,1:])
        self.edgedf = pd.DataFrame(columns = edge_var_list[1:], index = edgeindexarray, data = edge_data_array[:,1:])
        self.celldf = pd.DataFrame(columns = cell_var_list[1:], index = cellindexarray, data = cell_data_array[:,1:])

        # here we extract the locations of all the things we need
        vertexindexlist = self.vertexdf.columns.tolist()
        xindex = vertexindexlist.index('x')
        yindex = vertexindexlist.index('y')


        edgeindexlist = self.edgedf.columns.tolist()
        exflagindex = edgeindexlist.index('xflag')
        eyflagindex = edgeindexlist.index('yflag')
        vertex1index = edgeindexlist.index('vertex1')
        vertex2index = edgeindexlist.index('vertex2')
        lengthindex = edgeindexlist.index('length')
        cell1index = edgeindexlist.index('_cell1')
        cell2index = edgeindexlist.index('_cell2')

        cellindexlist = self.celldf.columns.tolist()
        A0index = cellindexlist.index('A0')
        #Kindex = cellindexlist.index('K')
        xflagindex = cellindexlist.index('xflag')
        yflagindex = cellindexlist.index('yflag')
        areaindex = cellindexlist.index('Area')
        CWVerticesindex = cellindexlist.index('CWVertices')
        EdgeListindex = cellindexlist.index('EdgeList')
        #try:
        #    PMindex = cellindexlist.index('Gamma')
        #except ValueError:
        #    PMindex = cellindexlist.index('Lambda')
        


        # assemble shapely cell objects
        cellPolygons = []
        cellkeepindexList = []
        #Klist = []
        A0list = []
        #PMlist = []
        cwunitveclist = []
        boundaryCellboolean = []
        boundaryCells = []
        for cellindex in range(0,len(self.celldf)):
            celldata = self.celldf.iloc[cellindex]
            if ( (celldata[xflagindex] + celldata[yflagindex]) == '00'):
                cellkeepindexList.append(self.celldf.index[cellindex])
                stringcwvertices = celldata[CWVerticesindex]
                cwverticesList = [int(j) for j in stringcwvertices.split(',')[:-1]]
                cellVertices = []
                for vert in cwverticesList:
                    pointx = float(self.vertexdf.loc[vert][xindex])
                    pointy = float(self.vertexdf.loc[vert][yindex])
                    point = [pointx, pointy]
                    cellVertices.append(point)
                polygon = Polygon(cellVertices)
                cellPolygons.append(polygon)
                #Klist.append(celldata[Kindex])
                A0list.append(celldata[A0index])
                #PMlist.append(celldata[PMindex])
                boundaryCellboolean.append(False)
            else:
                boundaryCells.append(cellindex)
        areas = [p.area for p in cellPolygons]
        for bc in boundaryCells:
            celldata = self.celldf.iloc[bc]
            stringcwvertices = celldata[CWVerticesindex]
            stringedges = celldata[EdgeListindex]
            cwverticesList = [int(j) for j in stringcwvertices.split(',')[:-1]]
            cwVertcoords = [[float(self.vertexdf.loc[firstvertex][xindex]), float(self.vertexdf.loc[firstvertex][yindex])] for firstvertex in cwverticesList]
            cwVertcoords = np.array(cwVertcoords + [cwVertcoords[0]])
            edgevectors = np.diff(cwVertcoords, axis=0)


            edgeList = [int(j) for j in stringedges.split(',')[:-1]]
            # let's align the edgeList with the cwverticesList
            ev1 = [int(self.edgedf.loc[e][vertex1index]) for e in edgeList]
            ev2 = [int(self.edgedf.loc[e][vertex2index]) for e in edgeList]
            reorderingedges = []
            cwverticesList_round = cwverticesList + [cwverticesList[0]]
            for ind in range(len(cwverticesList)):
                v1 = cwverticesList_round[ind]
                v2 = cwverticesList_round[ind+1]
                en = 0
                while (en < len(ev1)):
                    eeev1 = ev1[en]
                    eeev2 = ev2[en]
                    bool1 = ( (eeev1 == v1) and (eeev2 == v2))
                    bool2 = ( (eeev1 == v2) and (eeev2 == v1))
                    test = (bool1 or bool2)
                    if (test):
                        reorderingedges.append(en)
                        en = len(ev1)
                    en+=1

            #reorderingedges = [ev1.index(cwverticesList[i]) for i in range(len(cwverticesList))]
            cwedgeList = [edgeList[i] for i in reorderingedges]
            cwedgeXflagindicator = [int(self.edgedf.loc[e][exflagindex]) for e in cwedgeList]
            cwedgeYflagindicator = [int(self.edgedf.loc[e][eyflagindex]) for e in cwedgeList]
            
            XedgeCrossingDetermination = [ edgevectors[j][0] / np.abs(edgevectors[j][0])*cwedgeXflagindicator[j] for j in range(len(cwedgeXflagindicator))]
            YedgeCrossingDetermination = [ edgevectors[j][1] / np.abs(edgevectors[j][1])*cwedgeYflagindicator[j] for j in range(len(cwedgeXflagindicator))]

            if (np.sum(cwedgeXflagindicator) > 0):
                if (XedgeCrossingDetermination.index(-1) < XedgeCrossingDetermination.index(1)):
                    rightorleft = 1
                else:
                    rightorleft = -1
            else:
                rightorleft = 0

            if (np.sum(cwedgeYflagindicator) > 0):
                if (YedgeCrossingDetermination.index(-1) < YedgeCrossingDetermination.index(1)):
                    upordown = 1
                else:
                    upordown = -1
            else:
                upordown = 0
                    
            xtot = 0
            ytot = 0
            cwedgeXcrossdirec = []
            cwedgeYcrossdirec = []
            for i in range(len(cwedgeXflagindicator)):
                xtot = xtot + cwedgeXflagindicator[i]
                if ((cwedgeXflagindicator[i] == 1) and xtot > 1):
                    cwedgeXcrossdirec.append(-1)
                elif (xtot < 2):
                    cwedgeXcrossdirec.append(cwedgeXflagindicator[i])
            
            for i in range(len(cwedgeYflagindicator)):
                ytot = ytot + cwedgeYflagindicator[i]
                if ((cwedgeYflagindicator[i] == 1) and ytot > 1):
                    cwedgeYcrossdirec.append(-1)
                elif (ytot < 2):
                    cwedgeYcrossdirec.append(cwedgeYflagindicator[i])
            
            cwedgeXa = [(f % 2) for f in np.cumsum(cwedgeXflagindicator)]
            cwedgeYa = [(f % 2) for f in np.cumsum(cwedgeYflagindicator)]
                
            # YOU NEED SOME WAY TO TELL THE THE THING THAT THE SECOND CROSSING EDGE HAS TO CROSS IN THE OPPOSITE DIRECTION
            cell1Vertices = []
            cell2Vertices = []
            pointx = cwVertcoords[0][0]
            pointy = cwVertcoords[0][1]
            point = [pointx, pointy]
            cell1Vertices.append(point)
            maxx = 0
            maxy = 0
            minx = 2*self.Lx
            miny = 2*self.Ly
            for vv in range(0, len(cwverticesList)-1):        
                pointx = cwVertcoords[vv][0] + edgevectors[vv][0] + rightorleft*cwedgeXa[vv]*self.Lx
                pointy = cwVertcoords[vv][1] + edgevectors[vv][1] +  upordown*cwedgeYa[vv]*self.Ly
                point = [pointx, pointy]
                cell1Vertices.append(point)

            for vv in range(0, len(cwverticesList)):
                pointx = cell1Vertices[vv][0]
                pointy = cell1Vertices[vv][1]
                if (pointx > maxx):
                    maxx = pointx
                if (pointx < minx):
                    minx = pointx
                if (pointy > maxy):
                    maxy = pointy
                if (pointy < miny):
                    miny = pointy
            npcell1Vertices = np.array(cell1Vertices)
            xoffset = 0
            yoffset = 0
            if (int(celldata[xflagindex]) > 0):
                if (maxx > 0 and minx < 0):
                    xoffset = 1
                else:
                    xoffset = -1
            if (int(celldata[yflagindex]) > 0):
                if (maxy > 0 and miny < 0):
                    yoffset = 1
                else:
                    yoffset = -1
            npcell2Vertices = npcell1Vertices + np.array([xoffset*self.Lx, yoffset*self.Ly])
            polygon1 = Polygon(cell1Vertices)
            polygon2 = Polygon(npcell2Vertices)
            cellPolygons.append(polygon1)
            cellkeepindexList.append(self.celldf.index[bc])
            #Klist.append(celldata[Kindex])
            A0list.append(celldata[A0index])
            #PMlist.append(celldata[PMindex])
            areas.append(polygon1.area)
            boundaryCellboolean.append(True)
            cellPolygons.append(polygon2)
            cellkeepindexList.append(self.celldf.index[bc])
            #Klist.append(celldata[Kindex])
            A0list.append(celldata[A0index])
            #PMlist.append(celldata[PMindex])
            areas.append(polygon2.area)
            boundaryCellboolean.append(True)
            if (xoffset*yoffset != 0):
                npcell3Vertices = npcell1Vertices + np.array([xoffset*self.Lx, 0])
                polygon3 = Polygon(npcell3Vertices)
                cellPolygons.append(polygon3)
                cellkeepindexList.append(self.celldf.index[bc])
                #Klist.append(celldata[Kindex])
                A0list.append(celldata[A0index])
                #PMlist.append(celldata[PMindex])
                areas.append(polygon3.area)
                boundaryCellboolean.append(True)
                npcell4Vertices = npcell1Vertices + np.array([0, yoffset*self.Ly])
                polygon4 = Polygon(npcell4Vertices)
                cellPolygons.append(polygon4)
                cellkeepindexList.append(self.celldf.index[bc])
                #Klist.append(celldata[Kindex])
                A0list.append(celldata[A0index])
                #PMlist.append(celldata[PMindex])
                areas.append(polygon4.area)
                boundaryCellboolean.append(True)

        scaled_areas = [a / np.mean(areas) for a in areas]

        #make adjacency list
        neighborsList = [[] for i in range(len(cellPolygons))]

        for edgeindex in range(0, len(self.edgedf)):
            edgedata = self.edgedf.iloc[edgeindex]
            v1 = int(edgedata[vertex1index])
            v2 = int(edgedata[vertex2index])
            onecell = int(edgedata[cell1index])
            twocell = int(edgedata[cell2index])
            if ( (onecell in cellkeepindexList) and (twocell in cellkeepindexList) ):
                neighborsList[onecell].append(twocell)
                neighborsList[twocell].append(onecell)
        #neighborsListtrunc = [nonemptylist for nonemptylist in neighborsList if (len(nonemptylist) > 0)]
        # I'VE ONLY CALCULATED NEIGHBORS FOR THE NON BOUNDARY CROSSING CELLS - FOR NETWORKS
        #geopandas list (SVM = "shapely vertex model")
        areas = [float(a) for a in areas]
        scaled_areas = [float(a) for a in scaled_areas]
        #Klist = [float(k) for k in Klist]
        A0list = [float(a) for a in A0list]
        #PMlist = [float(p) for p in PMlist]
        first = gpd.GeoDataFrame(data = {'geometry':cellPolygons, 'VMindex':cellkeepindexList, 'neighbors':neighborsList, 'isitboundary': boundaryCellboolean, 'areas':areas, 'scaled areas':scaled_areas, 'A0':A0list})
        valid = []

        touching = []
        for cell in first.geometry:
            valid.append(cell.is_valid)
        second = first[valid]

        for k in second.geometry:
            touchcomb = False
            for j in second.geometry:
                touchcomb = touchcomb or k.touches(j)
            touching.append(touchcomb)
        if ( len(first) == len(second) ):
            self.isValid = True
        self.SVM = second[touching]


    def readFileMakeVars(self):
        # read the file 
        filetext= open(self.file, 'r')
        content = filetext.read()
        filetext.close()

        Lx_index = content.find("LX\t")
        #print(Lx_index)
        Lx_endindex = content[3+Lx_index:].find("\n")

        Ly_index = content.find("LY\t")
        Ly_endindex = content[3+Ly_index:].find("\n")
        self.Lx = float(content[3+Lx_index:3+ Lx_index + Lx_endindex])
        self.Ly = float(content[3+Ly_index:3+ Ly_index + Ly_endindex])
        P0_index = content.find("P0\t")
        P0_endindex = content[3+P0_index:].find("\n")
        self.P0 = float(content[3+P0_index:3+P0_index+P0_endindex])
        
        # find where the various sections go
        vertex_index = content.find("Vertices")
        edge_index = content.find("Edges\n")
        cell_index = content.find("Cells")
        end_index = content.find("End\n")

        # read the vertex bits
        vertex_substr = content[vertex_index:edge_index]
        vertex_list = vertex_substr.split("\n")
        #vertex_list[2:] # list of vertex data
        vert_var_list = np.array(vertex_list[1].split("\t"))
        vertex_data_array = np.array([subl.split("\t") for subl in vertex_list[2:-1]])

        # read the edge bits
        edge_substr = content[edge_index:cell_index]
        edge_list = edge_substr.split("\n")
        edge_var_list = np.array(edge_list[1].split("\t"))
        edge_data_array = np.array([subl.split("\t") for subl in edge_list[2:-1]])

        # read the cell bits
        cell_substr = content[cell_index:end_index]
        cell_list = cell_substr.split("\n")
        cell_var_list = np.array(cell_list[1].split("\t"))
        cell_data_array = np.array([subl.split("\t") for subl in cell_list[2:-1]])

        # make pandas arrays 
        vertindexarray = [int(v) for v in vertex_data_array[:, 0]]
        edgeindexarray = [int(e) for e in edge_data_array[:,0]]
        cellindexarray = [int(c) for c in cell_data_array[:,0]]
        self.vertexdf = pd.DataFrame(columns = vert_var_list[1:], index = vertindexarray, data = vertex_data_array[:,1:])
        self.edgedf = pd.DataFrame(columns = edge_var_list[1:], index = edgeindexarray, data = edge_data_array[:,1:])
        self.celldf = pd.DataFrame(columns = cell_var_list[1:], index = cellindexarray, data = cell_data_array[:,1:])

        # here we extract the locations of all the things we need
        vertexindexlist = self.vertexdf.columns.tolist()
        xindex = vertexindexlist.index('x')
        yindex = vertexindexlist.index('y')


        edgeindexlist = self.edgedf.columns.tolist()
        exflagindex = edgeindexlist.index('xflag')
        eyflagindex = edgeindexlist.index('yflag')
        vertex1index = edgeindexlist.index('vertex1')
        vertex2index = edgeindexlist.index('vertex2')
        lengthindex = edgeindexlist.index('length')
        cell1index = edgeindexlist.index('_cell1')
        cell2index = edgeindexlist.index('_cell2')

        cellindexlist = self.celldf.columns.tolist()
        A0index = cellindexlist.index('A0')
        Kindex = cellindexlist.index('K')
        xflagindex = cellindexlist.index('xflag')
        yflagindex = cellindexlist.index('yflag')
        areaindex = cellindexlist.index('Area')
        CWVerticesindex = cellindexlist.index('CWVertices')
        EdgeListindex = cellindexlist.index('EdgeList')
        try:
            PMindex = cellindexlist.index('Gamma')
        except ValueError:
            PMindex = cellindexlist.index('Lambda')
        


        # assemble shapely cell objects
        cellPolygons = []
        cellkeepindexList = []
        Klist = []
        A0list = []
        PMlist = []
        cwunitveclist = []
        boundaryCellboolean = []
        boundaryCells = []
        for cellindex in range(0,len(self.celldf)):
            celldata = self.celldf.iloc[cellindex]
            if ( (celldata[xflagindex] + celldata[yflagindex]) == '00'):
                cellkeepindexList.append(self.celldf.index[cellindex])
                stringcwvertices = celldata[CWVerticesindex]
                cwverticesList = [int(j) for j in stringcwvertices.split(',')[:-1]]
                cellVertices = []
                for vert in cwverticesList:
                    pointx = float(self.vertexdf.loc[vert][xindex])
                    pointy = float(self.vertexdf.loc[vert][yindex])
                    point = [pointx, pointy]
                    cellVertices.append(point)
                polygon = Polygon(cellVertices)
                cellPolygons.append(polygon)
                Klist.append(celldata[Kindex])
                A0list.append(celldata[A0index])
                PMlist.append(celldata[PMindex])
                boundaryCellboolean.append(False)
            else:
                boundaryCells.append(cellindex)
        areas = [p.area for p in cellPolygons]
        for bc in boundaryCells:
            celldata = self.celldf.iloc[bc]
            stringcwvertices = celldata[CWVerticesindex]
            stringedges = celldata[EdgeListindex]
            cwverticesList = [int(j) for j in stringcwvertices.split(',')[:-1]]
            cwVertcoords = [[float(self.vertexdf.loc[firstvertex][xindex]), float(self.vertexdf.loc[firstvertex][yindex])] for firstvertex in cwverticesList]
            cwVertcoords = np.array(cwVertcoords + [cwVertcoords[0]])
            edgevectors = np.diff(cwVertcoords, axis=0)


            edgeList = [int(j) for j in stringedges.split(',')[:-1]]
            # let's align the edgeList with the cwverticesList
            ev1 = [int(self.edgedf.loc[e][vertex1index]) for e in edgeList]
            ev2 = [int(self.edgedf.loc[e][vertex2index]) for e in edgeList]
            reorderingedges = []
            cwverticesList_round = cwverticesList + [cwverticesList[0]]
            for ind in range(len(cwverticesList)):
                v1 = cwverticesList_round[ind]
                v2 = cwverticesList_round[ind+1]
                en = 0
                while (en < len(ev1)):
                    eeev1 = ev1[en]
                    eeev2 = ev2[en]
                    bool1 = ( (eeev1 == v1) and (eeev2 == v2))
                    bool2 = ( (eeev1 == v2) and (eeev2 == v1))
                    test = (bool1 or bool2)
                    if (test):
                        reorderingedges.append(en)
                        en = len(ev1)
                    en+=1

            #reorderingedges = [ev1.index(cwverticesList[i]) for i in range(len(cwverticesList))]
            cwedgeList = [edgeList[i] for i in reorderingedges]
            cwedgeXflagindicator = [int(self.edgedf.loc[e][exflagindex]) for e in cwedgeList]
            cwedgeYflagindicator = [int(self.edgedf.loc[e][eyflagindex]) for e in cwedgeList]
            
            XedgeCrossingDetermination = [ edgevectors[j][0] / np.abs(edgevectors[j][0])*cwedgeXflagindicator[j] for j in range(len(cwedgeXflagindicator))]
            YedgeCrossingDetermination = [ edgevectors[j][1] / np.abs(edgevectors[j][1])*cwedgeYflagindicator[j] for j in range(len(cwedgeXflagindicator))]

            if (np.sum(cwedgeXflagindicator) > 0):
                if (XedgeCrossingDetermination.index(-1) < XedgeCrossingDetermination.index(1)):
                    rightorleft = 1
                else:
                    rightorleft = -1
            else:
                rightorleft = 0

            if (np.sum(cwedgeYflagindicator) > 0):
                if (YedgeCrossingDetermination.index(-1) < YedgeCrossingDetermination.index(1)):
                    upordown = 1
                else:
                    upordown = -1
            else:
                upordown = 0
                    
            xtot = 0
            ytot = 0
            cwedgeXcrossdirec = []
            cwedgeYcrossdirec = []
            for i in range(len(cwedgeXflagindicator)):
                xtot = xtot + cwedgeXflagindicator[i]
                if ((cwedgeXflagindicator[i] == 1) and xtot > 1):
                    cwedgeXcrossdirec.append(-1)
                elif (xtot < 2):
                    cwedgeXcrossdirec.append(cwedgeXflagindicator[i])
            
            for i in range(len(cwedgeYflagindicator)):
                ytot = ytot + cwedgeYflagindicator[i]
                if ((cwedgeYflagindicator[i] == 1) and ytot > 1):
                    cwedgeYcrossdirec.append(-1)
                elif (ytot < 2):
                    cwedgeYcrossdirec.append(cwedgeYflagindicator[i])
            
            cwedgeXa = [(f % 2) for f in np.cumsum(cwedgeXflagindicator)]
            cwedgeYa = [(f % 2) for f in np.cumsum(cwedgeYflagindicator)]
                
            # YOU NEED SOME WAY TO TELL THE THE THING THAT THE SECOND CROSSING EDGE HAS TO CROSS IN THE OPPOSITE DIRECTION
            cell1Vertices = []
            cell2Vertices = []
            pointx = cwVertcoords[0][0]
            pointy = cwVertcoords[0][1]
            point = [pointx, pointy]
            cell1Vertices.append(point)
            maxx = 0
            maxy = 0
            minx = 2*self.Lx
            miny = 2*self.Ly
            for vv in range(0, len(cwverticesList)-1):        
                pointx = cwVertcoords[vv][0] + edgevectors[vv][0] + rightorleft*cwedgeXa[vv]*self.Lx
                pointy = cwVertcoords[vv][1] + edgevectors[vv][1] +  upordown*cwedgeYa[vv]*self.Ly
                point = [pointx, pointy]
                cell1Vertices.append(point)

            for vv in range(0, len(cwverticesList)):
                pointx = cell1Vertices[vv][0]
                pointy = cell1Vertices[vv][1]
                if (pointx > maxx):
                    maxx = pointx
                if (pointx < minx):
                    minx = pointx
                if (pointy > maxy):
                    maxy = pointy
                if (pointy < miny):
                    miny = pointy
            npcell1Vertices = np.array(cell1Vertices)
            xoffset = 0
            yoffset = 0
            if (int(celldata[xflagindex]) > 0):
                if (maxx > 0 and minx < 0):
                    xoffset = 1
                else:
                    xoffset = -1
            if (int(celldata[yflagindex]) > 0):
                if (maxy > 0 and miny < 0):
                    yoffset = 1
                else:
                    yoffset = -1
            npcell2Vertices = npcell1Vertices + np.array([xoffset*self.Lx, yoffset*self.Ly])
            polygon1 = Polygon(cell1Vertices)
            polygon2 = Polygon(npcell2Vertices)
            cellPolygons.append(polygon1)
            cellkeepindexList.append(self.celldf.index[bc])
            Klist.append(celldata[Kindex])
            A0list.append(celldata[A0index])
            PMlist.append(celldata[PMindex])
            areas.append(polygon1.area)
            boundaryCellboolean.append(True)
            cellPolygons.append(polygon2)
            cellkeepindexList.append(self.celldf.index[bc])
            Klist.append(celldata[Kindex])
            A0list.append(celldata[A0index])
            PMlist.append(celldata[PMindex])
            areas.append(polygon2.area)
            boundaryCellboolean.append(True)
            if (xoffset*yoffset != 0):
                npcell3Vertices = npcell1Vertices + np.array([xoffset*self.Lx, 0])
                polygon3 = Polygon(npcell3Vertices)
                cellPolygons.append(polygon3)
                cellkeepindexList.append(self.celldf.index[bc])
                Klist.append(celldata[Kindex])
                A0list.append(celldata[A0index])
                PMlist.append(celldata[PMindex])
                areas.append(polygon3.area)
                boundaryCellboolean.append(True)
                npcell4Vertices = npcell1Vertices + np.array([0, yoffset*self.Ly])
                polygon4 = Polygon(npcell4Vertices)
                cellPolygons.append(polygon4)
                cellkeepindexList.append(self.celldf.index[bc])
                Klist.append(celldata[Kindex])
                A0list.append(celldata[A0index])
                PMlist.append(celldata[PMindex])
                areas.append(polygon4.area)
                boundaryCellboolean.append(True)

        scaled_areas = [a / np.mean(areas) for a in areas]

        #make adjacency list
        neighborsList = [[] for i in range(len(cellPolygons))]

        for edgeindex in range(0, len(self.edgedf)):
            edgedata = self.edgedf.iloc[edgeindex]
            v1 = int(edgedata[vertex1index])
            v2 = int(edgedata[vertex2index])
            onecell = int(edgedata[cell1index])
            twocell = int(edgedata[cell2index])
            if ( (onecell in cellkeepindexList) and (twocell in cellkeepindexList) ):
                neighborsList[onecell].append(twocell)
                neighborsList[twocell].append(onecell)
        #neighborsListtrunc = [nonemptylist for nonemptylist in neighborsList if (len(nonemptylist) > 0)]
        # I'VE ONLY CALCULATED NEIGHBORS FOR THE NON BOUNDARY CROSSING CELLS - FOR NETWORKS
        #geopandas list (SVM = "shapely vertex model")
        areas = [float(a) for a in areas]
        scaled_areas = [float(a) for a in scaled_areas]
        Klist = [float(k) for k in Klist]
        A0list = [float(a) for a in A0list]
        PMlist = [float(p) for p in PMlist]
        first = gpd.GeoDataFrame(data = {'geometry':cellPolygons, 'VMindex':cellkeepindexList, 'neighbors':neighborsList, 'isitboundary': boundaryCellboolean, 'areas':areas, 'scaled areas':scaled_areas, 'K':Klist, 'A0':A0list, 'PM': PMlist})
        valid = []

        touching = []
        for cell in first.geometry:
            valid.append(cell.is_valid)
        second = first[valid]

        for k in second.geometry:
            touchcomb = False
            for j in second.geometry:
                touchcomb = touchcomb or k.touches(j)
            touching.append(touchcomb)
        if ( len(first) == len(second) ):
            self.isValid = True
        self.SVM = second[touching]
        tensions = (np.array(self.SVM.geometry.length) - self.P0)*np.array(self.SVM.PM, dtype = 'float')
        pressures = (self.SVM.A0 - self.SVM.areas)*self.SVM.K
        self.SVM['pressures'] = pressures
        self.SVM['tensions'] = tensions
        normaltensionstress = []
        normalstress = []
        for index in range(0, len(self.SVM)):
            vmindex = self.SVM.iloc[index]["VMindex"]
            area = self.SVM.iloc[index]["areas"]
            stredgelist = self.celldf.loc[vmindex]["EdgeList"]
            edited = stredgelist.split(',')
            edgelist = [int(e) for e in edited[:-1]]
            tension = np.array([float(self.edgedf.loc[ed]["tension"]) for ed in edgelist])
            edgelength = np.array([float(self.edgedf.loc[ed]["length"]) for ed in edgelist])
            normaltensionstress.append((0.5 / area) * (np.dot(tension, edgelength))) 
            normalstress.append((0.5 / area) * (np.dot(tension, edgelength)) - self.SVM.iloc[index]["pressures"])
        self.SVM['normaltensionstress'] = np.array(normaltensionstress)
        self.SVM['normalstress'] = np.array(normalstress)

    def calculateHomogeneousNondimEnergy(self): # requires use of readFileMakeVars_HOMOGENEOUS() first
        Ncells = np.max(self.SVM.VMindex.tolist())
        isitused = np.zeros((Ncells), dtype='bool')
        energy = 0
        area_diff = 0
        perim_diff = 0
        for i in range(len(self.SVM)):
            cellindex = self.SVM.VMindex.iloc[i]
            if (not isitused[cellindex-1]):
                area = self.SVM.iloc[i].geometry.area
                perimeter = self.SVM.iloc[i].geometry.length
                A0 = self.SVM.iloc[i].A0
                P0 = self.P0
                K = self.K
                Lambda = self.Lambda
                Lx = self.Lx
                Ly = self.Ly
                area_diff += (area/A0 - 1)**2  
                perim_diff += (perimeter - P0)**2 
                isitused[cellindex-1] = True
        #SVtest.SVM.iloc[0].geometry.length
        energy = (K*A0*area_diff*.5 + Lambda*perim_diff*.5)
        normfactor = (K * Lx**2 * Ly**2) / Ncells
        energy_SC = (K*A0*area_diff*.5 + Lambda*perim_diff*.5) / normfactor
        return energy_SC

    def calculateHomogeneousSQUEEZEDNondimEnergy(self, xscale): # requires use of readFileMakeVars_HOMOGENEOUS() first
        Ncells = np.max(self.SVM.VMindex.tolist())
        isitused = np.zeros((Ncells), dtype='bool')
        energy = 0
        area_diff = 0
        perim_diff = 0
        newgs = self.SVM.affine_transform([xscale, 0, 0,  1/xscale, 0, 0])
        for i in range(len(self.SVM)):
            cellindex = self.SVM.VMindex.iloc[i]
            if (not isitused[cellindex-1]):
                area = newgs[i].area
                perimeter = newgs[i].length
                A0 = self.SVM.iloc[i].A0
                P0 = self.P0
                K = self.K
                Lambda = self.Lambda
                Lx = self.Lx
                Ly = self.Ly
                area_diff += (area/A0 - 1)**2  
                perim_diff += (perimeter - P0)**2 
                isitused[cellindex-1] = True
        #SVtest.SVM.iloc[0].geometry.length
        energy = (K*A0*area_diff*.5 + Lambda*perim_diff*.5)
        normfactor = (K * Lx**2 * Ly**2) / Ncells
        energy_SC = (K*A0*area_diff*.5 + Lambda*perim_diff*.5) / normfactor
        return energy_SC
    
    def getPiecewise(self, axmax , aymax , fill_val = 'scaled areas'):
        unionshape = gpd.GeoDataFrame(geometry = gpd.GeoSeries(self.SVM.buffer(0).unary_union))
        xmin = 0#int( np.round( unionshape.bounds['minx'].min() * ( 1 - bounds_excess ) ) )
        ymin = 0#int( np.round( unionshape.bounds['miny'].min() * ( 1 - bounds_excess ) ) )
        xmax = int(axmax)#int( np.round( unionshape.bounds['maxx'].max() * ( 1 + bounds_excess ) ) )
        ymax = int(aymax)#int( np.round( unionshape.bounds['maxy'].max() * ( 1 + bounds_excess ) ) )
        funcVals = np.zeros(shape=(ymax - ymin, xmax - xmin)) 
        gpdlist = self.SVM
        polylist = self.SVM.affine_transform([(xmax-xmin) / self.Lx, 0, 0, (ymax-ymin) / self.Ly, -xmin, -ymin]).geometry.values#translate(xoff = -xmin, yoff = -ymin).geometry.values
        for i in range(len(polylist)):
            if (polylist[i].geom_type != 'LineString' and polylist[i].geom_type != 'Point'):    
            #vertices = np.array([self.vertices[v] for v in vertList])
                if (polylist[i].geom_type == 'MultiPolygon'):
                # THIS SECTION CAN HELP CLEAN UP THE NEIGHBORPOLYGON LIST OF POLYGONS - USE THIS IN FUTURE!!!
                    MP = polylist[i]
                    area = MP.area
                    mlist = list(MP.geoms)
                    arearatiolist = [p.area / area for p in mlist]
                    polynew = mlist[np.argmax(arearatiolist)]
                    A = polynew.exterior.coords.xy
                elif (polylist[i].geom_type != 'LineString'):
                    A = polylist[i].exterior.coords.xy
                rr, cc = polygon(A[1], A[0] , funcVals.shape)
                if fill_val in self.SVM.columns:
                    funcVals[rr, cc] = gpdlist.iloc[i][fill_val]
                else:
                    funcVals[rr, cc] = 1 
        return funcVals
    
    # CREATE A FUNCTION THAT REPLICATES THE EQ FILES BUT WITH VARIANCES IN K, A0, K/GAMMA
    def makeInhomogeneousTissue(self, A0, A0varratio, Kvarratio, LAMBDAvarratio):
        # read the file 
        filetext= open(self.file, 'r')
        content = filetext.read()
        filetext.close()
        lines = content.split('\n')
        endglobalindex = lines.index("endGlobal")
        
        A0index = lines.index('A0\t0')
        Kstdratioindex = lines.index('Kstdratio\t0')
        A0stdratioindex = lines.index('A0stdratio\t0')
        LABMDAstdratioindex = lines.index('LAMBDAstdratio\t0')


        A0split = lines[A0index].split('\t')
        A0split[1] = str(A0)
        lines[A0index] = A0split
        Kstdratiosplit = lines[Kstdratioindex].split('\t')
        Kstdratiosplit[1] = Kvarratio
        lines[Kstdratioindex] = Kstdratiosplit
        A0stdratiosplit = lines[A0stdratioindex].split('\t')
        A0stdratiosplit[1] = A0varratio

        splitus = self.file.split('_')
        
        A0findex = splitus.index('A0')
        splitus[A0findex+1] = str(A0)

        Kstdfindex = splitus.index('Kvarratio')
        splitus[Kstdfindex+1] = str(Kvarratio)

        A0stdfindex = splitus.index('A0varratio')
        splitus[A0stdfindex+1] = str(A0varratio)

        LAMBDAstdfindex = splitus.index('LAMBDAvarratio')
        splitus[LAMBDAstdfindex] = str(LAMBDAvarratio)
        new_name = ''
        for j in range(len(splitus[-1])):
            new_name = new_name + splitus[j] + '_'
        new_name = new_name + splitus[-1]

        newlines = ''
        for j in range(len(lines)):
            newlines = newlines + str(lines[j]) + '\n'
        
        new_file = open(new_name, "w")
        new_file.write(newlines)
        new_file.close()

    