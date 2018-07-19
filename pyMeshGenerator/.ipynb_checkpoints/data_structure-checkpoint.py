from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection




# Create points array
#------------------------------------------------
def global_points(N, rad, ang):
    from numpy import sin, cos, pi
    points = []
    ID = 0
    for j in range(1,N[0]):
        for i in range(N[1]):
            r     = j * rad / (N[1]-1)
            omega = i * ang / N[0] + j*pi
            points.append(Point(sin(2.0*omega*pi)*r, 
                                cos(2.0*omega*pi)*r,
                                ID))
            ID += 1
    return points



#------------------------------------------------
class Point:
    def __init__(self, x, y, ID):
        '''
        Class defines a single point in the mesh
        '''
        self.x = x
        self.y = y
        self.ID = ID
    

        
# Create quadtree for point dataset
#------------------------------------------------
class PointQuadtree:
    def __init__(self, points, bbox=None, **kwargs):
        
        # Determine bounding box 
        if (bbox==None):
            self.getBoundingBox(points)
        else:
            x_bb, y_bb = bbox[0], bbox[1]
            self.x_bbox = x_bb
            self.y_bbox = y_bb
            self.center = (0.5*sum(x_bb), 0.5*sum(y_bb))
            self.width  = x_bb[1]-x_bb[0]
            self.height = y_bb[1]-y_bb[0]
            self.bbox_IDs = [None, None, None, None]
        
        self.NorthEast = None
        self.NorthWest = None
        self.SouthWest = None
        self.SouthEast = None
        self.pointLimit = kwargs.get('limit', 4)
        self.points = []
        self.splitted = False
        
        # Distribute points
        pnt_ind = 0
        while not self.splitted and (pnt_ind < len(points)):
            
            if (len(self.points) >= self.pointLimit):
                self.splitTree(points)
            else:
                p = points[pnt_ind]
                if self.isInside(p):
                    self.points.append(p)
                pnt_ind += 1
    
    def splitTree(self, points):
        self.splitted = True
        self.points = []
        self.NorthEast = PointQuadtree(points, 
                                       limit=self.pointLimit, 
                                       bbox=[[self.center[0], self.center[0]+0.5*self.width],
                                             [self.center[1], self.center[1]+0.5*self.height]])
        self.NorthWest = PointQuadtree(points, 
                                       limit=self.pointLimit, 
                                       bbox=[[self.center[0]-0.5*self.width, self.center[0]],
                                             [self.center[1], self.center[1]+0.5*self.height]])
        self.SouthWest = PointQuadtree(points, 
                                       limit=self.pointLimit, 
                                       bbox=[[self.center[0]-0.5*self.width,  self.center[0]],
                                             [self.center[1]-0.5*self.height, self.center[1]]])
        self.SouthEast = PointQuadtree(points, 
                                       limit=self.pointLimit, 
                                       bbox=[[self.center[0], self.center[0]+0.5*self.width],
                                             [self.center[1]-0.5*self.height, self.center[1]]])
        
    def isInside(self, p):
        lowLimX = p.x >= self.center[0]-0.5*self.width
        hiLimX  = p.x <= self.center[0]+0.5*self.width
        lowLimY = p.y >= self.center[1]-0.5*self.height
        hiLimY  = p.y <= self.center[1]+0.5*self.height
        return (lowLimX and hiLimX and lowLimY and hiLimY)
    
    def doOverlap(self, l1, r1, l2, r2):
        if (l1[0] > r2[0] or l2[0] > r1[0]):
            return False
        if (l1[1] < r2[1] or l2[1] < r1[1]):
            return False
        return True
    
    def getPointsInRange(self, center, width, height):
        # Check if region overlaps with quadtree region
        l1 = (center[0]-0.5*width, center[1]+0.5*height)
        r1 = (center[0]+0.5*width, center[1]-0.5*height)
        l2 = (self.center[0]-0.5*self.width, self.center[1]+0.5*self.height)
        r2 = (self.center[0]+0.5*self.width, self.center[1]-0.5*self.height)
        
        pts = []
        
        if not self.doOverlap(l1,r1,l2,r2):
            return pts
        else:
            if self.splitted:
                pts.extend(self.NorthWest.getPointsInRange(center, width, height))
                pts.extend(self.NorthEast.getPointsInRange(center, width, height))
                pts.extend(self.SouthEast.getPointsInRange(center, width, height))
                pts.extend(self.SouthWest.getPointsInRange(center, width, height))
                
            else:
                if len(self.points) > 0:
                    for p in self.points:
                        pInRange = (p.x>=l1[0] and p.x<=r1[0] and p.y<=l1[1] and p.y>=r1[1])
                        if pInRange:
                            pts.append(p)

            return pts
                
    def getBoundingBox(self, points):
        x_bb, y_bb = [1.0e30,-1.0e30], [1.0e30,-1.0e30]
        x_bb_ID, y_bb_ID = [0,0], [0,0]
        for p in points:
            if p.x < x_bb[0]:
                x_bb[0] = p.x
            if p.x > x_bb[1]:
                x_bb[1] = p.x
            if p.y < y_bb[0]:
                y_bb[0] = p.y
            if p.y > y_bb[1]:
                y_bb[1] = p.y
                
        self.x_bbox = x_bb
        self.y_bbox = y_bb
        
        # Add bounding box points to global point array
        nPts = len(points)
        points.append(Point(x_bb[0], y_bb[0], nPts))
        points.append(Point(x_bb[1], y_bb[0], nPts+1))
        points.append(Point(x_bb[1], y_bb[1], nPts+2))
        points.append(Point(x_bb[0], y_bb[1], nPts+3))
        self.bbox_IDs = [nPts, nPts+1, nPts+2, nPts+3]
        
        self.center = (0.5*sum(x_bb), 0.5*sum(y_bb))
        self.width  = x_bb[1]-x_bb[0]
        self.height = y_bb[1]-y_bb[0]
        
        
    def plotQuadTree(self, patches):
        rect = mpatches.Rectangle((self.center[0]-0.5*self.width,
                                   self.center[1]-0.5*self.height), 
                                   self.width, self.height, 
                                   fc='none', ec='k', fill=0)
        patches.append(rect)
        if (self.NorthEast != None):
            self.NorthEast.plotQuadTree(patches)
        if (self.NorthWest != None):
            self.NorthWest.plotQuadTree(patches)
        if (self.SouthWest != None):
            self.SouthWest.plotQuadTree(patches)
        if (self.SouthEast != None):
            self.SouthEast.plotQuadTree(patches)