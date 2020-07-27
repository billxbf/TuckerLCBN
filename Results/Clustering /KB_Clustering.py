import numpy as np
import pandas as pd    
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial import distance_matrix


# return parallel minima of vectors
def pmin(dx, dy):
    n = len(dx)
    tmp = np.zeros(n)
    for i in range(n):
        if dx[i] < dy[i]:
            tmp[i] = dx[i]
        else:
            tmp[i] = dy[i]
    return tmp

# return points within the union of x and y
def which(boolarr):
    tmp = []
    for i in range(len(boolarr)):
        if boolarr[i]:
            tmp += [i]
    return tmp

# return matrix of Bx,z
def get_cont_mat(D):
    #D: Distance matrix, in csv format, transformed into numpy array
    #return B in numpy array 
    assert (np.shape(D)[0] == np.shape(D)[1])
    n = D.shape[0]
    B = np.zeros((n,n))
    for x in range(n):
        for y in range(n):
            dxy = D[x,y]
            dx = D[x]
            dy = D[y]
            Uxy = which(pmin(dx,dy) <= dxy) # find union of x,y
            wx =  1*(dx[Uxy] < dy[Uxy]) # label points closer to x as 1, otherwise 0
            B[x,Uxy] += 1/len(Uxy)*wx
    return (B/n)
                
   
class Graph: 
    # init function to declare class variables 
    def __init__(self,V): 
        self.V = V # number of vertices
        self.adj = [[] for i in range(V)] # adjacent list
  
    # do depth first search to find connected components
    def DFSUtil(self, temp, v, visited): 
        visited[v] = True # mark the current vertex as visited 
        temp.append(v)
        for i in self.adj[v]: 
            if visited[i] == False: 
                temp = self.DFSUtil(temp, i, visited) # update the list 
        return temp 
  
    # add an undirected edge 
    def addEdge(self, v, w): 
        self.adj[v].append(w) 
        self.adj[w].append(v) 
  
    # retrieve connected components in an undirected graph 
    def connectedComponents(self): 
        visited = [] 
        cc = [] 
        for i in range(self.V): 
            visited.append(False) 
        for v in range(self.V): 
            if visited[v] == False: 
                temp = [] 
                cc.append(self.DFSUtil(temp, v, visited)) 
        return cc 
            
def get_clusters(B):
    n = B.shape[0]
    Exp = 0 # Expected Bx,z
    for i in np.diag(B):
        Exp += i
    Exp = Exp/(2*n)

    g = Graph(n) 
    for x in range(n):
        for y in range(x+1,n):
            if min(B[x,y],B[y,x]) >= Exp:
                g.addEdge(x, y)
                cc = g.connectedComponents() 
    return cc, g.adj


def get_local_depth(B):
    LD = []
    n = B.shape[0]
    for i in range(n):
        LD += [round(np.sum(B[i,]),2)]
    
    return LD    
        

#%%
sns.color_palette(n_colors=1)
#Data = pd.read_csv('sample_data.csv')
#Data = Data.values
D = distance_matrix(df,Data)

B = get_cont_mat(D)
LD = get_local_depth(B) 
CC, E = get_clusters(B)    

for c in range(len(CC)):
    df = pd.DataFrame(Data[CC[c],], columns = ['x','y'])
    plt.scatter(x = "x", y = "y", data = df)

for i,ld in enumerate(LD):
    x = Data[i,0]
    y = Data[i,1]
    plt.text(x+0.3, y+0.3, ld, fontsize=9)
    e = E[i]
    for adj in e:
        if adj > i:
            plt.plot([Data[i,0], Data[adj,0]], [Data[i,1], Data[adj,1]], 
                     color = 'gray', linewidth=1)
plt.title('Parameter-Free Clustering via Partitioned Local Depth')
plt.savefig('Clustering.png')

for c in range(len(CC)):
    print('Cluster', c, ':',CC[c])
    
        
        