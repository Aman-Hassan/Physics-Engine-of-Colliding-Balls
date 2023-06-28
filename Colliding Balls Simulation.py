from cmath import inf

class Heap:
    class Node:
        def __init__(self,key,parent=None,r_child=None,l_child=None) -> None:
            self._k = key    #3 tuple of the form (i,t,T') where t is the time at which new collision can occur and T' is the time at which last collision for that block occured
            self._p = parent #parent of preset node
            self._rc = r_child #right child
            self._lc = l_child #left child
    
    def __init__(self):   #L contains the list of values
        self._heap = []  #This is list of pointers according to our reuired heap structure
        self._pointerlist= [1] #This is list of pointers according to index i
        self._n = 0  #no of nodes = 0 initially
        self._root = None #root of heap 
        self._tail = None #last element of heap
        
    
    def isempty(self):
        if (self.n == 0):
            return(True)
        else:
            return (False)
    
    def heapup(self,u):
        u1 = u
        
        if u1._p == None:
            return()

        while u1._p._k[1] >= u1._k[1]:
            if u1._p._k[1] == u1._k[1] and u1._p._k[0] > u1._k[0]:
                i1 = u1._k[0]
                i2 = u1._p._k[0]
                u1._p._k,u1._k= u1._k,u1._p._k
                self._pointerlist[i1+1],self._pointerlist[i2+1] = self._pointerlist[i2+1],self._pointerlist[i1+1]
                u1 = u1._p
            elif u1._p._k[1] == u1._k[1] and u1._p._k[0] < u1._k[0]:
                break
            else:
                i1 = u1._k[0]
                i2 = u1._p._k[0]
                u1._p._k,u1._k= u1._k,u1._p._k
                self._pointerlist[i1+1],self._pointerlist[i2+1] = self._pointerlist[i2+1],self._pointerlist[i1+1]
                u1 = u1._p
            if u1._p == None:
                break
        return()

    def heapdown(self,u):
        u1 = u
        v = None

        if u1._rc == None and u1._lc != None:
                v = u1._lc
                if v._k[1] < u1._k[1]:
                    i1 = u1._k[0]
                    i2 = v._k[0]
                    v._k,u1._k= u1._k,v._k
                    self._pointerlist[i1+1],self._pointerlist[i2+1] = self._pointerlist[i2+1],self._pointerlist[i1+1]
                    return()
                elif v._k[1] == u1._k[1] and v._k[0] < u1._k[0]:
                    i1 = u1._k[0]
                    i2 = v._k[0]
                    v._k,u1._k= u1._k,v._k
                    self._pointerlist[i1+1],self._pointerlist[i2+1] = self._pointerlist[i2+1],self._pointerlist[i1+1]
                    return()
                return()
        elif u1._rc == None and u1._lc == None:
            return()
       
        
        while u1._rc._k[1] <= u1._k[1] or u1._lc._k[1] <= u1._k[1]:
            if u1._rc._k[1] < u1._lc._k[1]:
                v = u1._rc
            elif u1._rc._k[1] > u1._lc._k[1]:
                v = u1._lc
            elif u1._rc._k[1] == u1._lc._k[1] and u1._rc._k[0] > u1._lc._k[0]:
                v = u1._lc
            else:
                v = u1._rc
            
            i1 = u1._k[0]
            i2 = v._k[0]

            if v._k[1] == u1._k[1] and v._k[0] < u1._k[0]:
                v._k,u1._k= u1._k,v._k
                u1 = v
                self._pointerlist[i1+1],self._pointerlist[i2+1] = self._pointerlist[i2+1],self._pointerlist[i1+1]
            elif v._k[1] == u1._k[1] and v._k[0] > u1._k[0]:
                break
            else:
                v._k,u1._k= u1._k,v._k
                u1 = v
                self._pointerlist[i1+1],self._pointerlist[i2+1] = self._pointerlist[i2+1],self._pointerlist[i1+1]
            if v._rc == None or v._lc == None:
                break
        
        if u1._rc == None and u1._lc != None:
                v = u1._lc
                i1 = u1._k[0]
                i2 = v._k[0]
                if v._k[1] < u1._k[1]:
                    v._k,u1._k= u1._k,v._k
                    self._pointerlist[i1+1],self._pointerlist[i2+1] = self._pointerlist[i2+1],self._pointerlist[i1+1]
                    return()
                elif v._k[1] == u1._k[1] and v._k[0] < u1._k[0]:
                    v._k,u1._k= u1._k,v._k
                    self._pointerlist[i1+1],self._pointerlist[i2+1] = self._pointerlist[i2+1],self._pointerlist[i1+1]
                    return()
                return()
        else:
            return()
        
        
    
    # def minofHeap(self):
    #     m = self._root._k
    #     index = m[0]
    #     t = self._heap.pop()
    #     tp = t._p
    #     if t == tp._rc:
    #         tp._rc = self._root
    #         tp._rc._k[1] = float(inf)
    #         tp._rc._k[2] = m[1]
    #         self._heap.append(tp._rc)
    #     else:
    #         tp._lc = self._root
    #         tp._lc._k[1] = float(inf)
    #         tp._lc._k[2] = m[1]
    #         self._heap.append(tp._lc)
    #     # if index + 1 <= self._n -1:
    #     #     self._pointerlist[index+1]._k[2] = m[1]
    #     self._root._k = t._k
    #     self.heapdown(self._root)
    #     return(m)
    
    def changekey(self,u,v): #u is node, v is value #Can be optimised by seperating change key to each of i t and t' (eg: in case of change of t' we dont need heap up or heapdown)
        x = u._k
        u._k = v
    
        if v[1]> x[1]:
            self.heapdown(u)
        else:
            self.heapup(u)
        return()

    def BuildHeap(self,L):
        l = len(L)
        self._n = l #no of nodes = 0 initially
        
        for i in range(l):
            u = self.Node(L[i])
            if i == 0:
                u._p= None
            else:
                u._p = self._heap[(i-1)//2]
            self._heap.append(u)
        
        for j in range(l):
            if 2*j + 2 <= l-1:
                self._heap[j]._lc = self._heap[2*j+1]
                self._heap[j]._rc = self._heap[2*j+2]
            elif 2*j + 1 <= l-1:
                self._heap[j]._lc = self._heap[2*j+1]
            else:
                break

        self._pointerlist = self._pointerlist + self._heap
        
        for p in range((l-2)//2,-1,-1):
            self.heapdown(self._heap[p])

        self._root = self._heap[0] #root of heap 
        self._tail = self._heap[l-1] #last element of heap
        
        return(self._pointerlist) 
    
    def print_tree(self):   #O(n)
        i=0
        import math
        # print()
        while i<self._n:
            print(self._heap[i]._k,end=" ")
            k=math.log(i+2,2)
            if k==float(int(k)):
                print()
            i+=1
      
    


# L = [(0,1),(1,4),(2,2),(3,4),(4,7),(5,9),(6,2),(9,3),(8,6),(7,5)]
# H1 = Heap()
# H = Heap.BuildHeap(H1,L)
# Heap.print_tree(H1)

def listCollisions(M,x,v,m,T): 
    L = []
    H1 = Heap()
    m1=0
    t=0
    # n=0
    l = len(M)

    d = []
    for i in range(l-1):
        if (v[i]-v[i+1]) > 0:
            t1 = (x[i+1]-x[i])/(v[i]-v[i+1])
            d.append([i,t1,float(0)]) #Initial list of 3-tuple required to create heap, NOTE: for collision to occur v[i]-v[i+1]>0
        else:
            d.append([i,float(inf),float(0)])
    d.append([i,float(inf),float(0)])
    L1 = Heap.BuildHeap(H1,d) #List of pointers based on index
    # Heap.print_tree(H1)
    
    i = 1
    while m1<m and t<=T:
        # print("iteration number:",i)
        # Heap.print_tree(H1)
        # print()


        #This part updates velocities and positions of the collision that will occur first in heap
        minimum = H1._root._k
        k = minimum[0] #index for collision
        t = minimum[1]  #time for collision
        t1 = minimum[2] #time of previous collision
        if t > T:
            break
        L.append((round(t,4),k,round(x[k]+v[k]*(t-t1),4)))
        msum = M[k] + M[k+1]
        mdif = M[k] - M[k+1]
        x_now=x[k]+v[k]*(t-t1)
        x[k] = x_now
        # x[k+1] = x_now
        v[k],v[k+1] = ( (mdif/msum)*v[k] + ((2*M[k+1])/msum)*v[k+1] ), ( ((2*M[k])/msum)*v[k] - (mdif/msum)*v[k+1] )
        x[k+1] = x_now #- v[k+1]*(t-H1._pointerlist[k+2]._k[2])
        
        # print("X-list:",x)
        # print("V-List:",v)
        Heap.changekey(H1,H1._root,[k,float(inf),t])
        m1+=1
        i+=1

        #This part updates values of possible time of collisions of k-1 and k+1
        p = 0

        # print("Working with index:",k)
        # print("time of collision:",t)
        # print("time of previous collision:",t1)
        # print()
        #for k-1
        if k-1 >= 0:
            point_k1 = H1._pointerlist[k]  #gives the pointer to the index correspoding to k-1th collision
            tk1_old = point_k1._k[2]         #gives the time of previous collision?
            if (v[k-1]-v[k]) > 0:
                tk1 = t + ((x[k]-(x[k-1] + v[k-1]*(t-tk1_old)))/(v[k-1]-v[k])) #updating the next possible collision time
            else:
                tk1 = float(inf)
            H1.changekey(point_k1,[k-1,tk1,tk1_old])  #changing the key , index and new time of collision is fine...what is old time of collision now?
            p = p+1

        #for k+1
        if k+2 <= l-1:
            point_k2 = H1._pointerlist[k+2]  #gives the pointer to index corresponding to k+1th collision
            point_k3 = H1._pointerlist[k+3]
            tk3_old = point_k3._k[2]
            if (v[k+1]-v[k+2]) > 0:
                # print("numerator:", ((x[k+2] + v[k+2]*(t-tk2_old))-x_now))
                tk2 = t + (((x[k+2] + v[k+2]*(t-tk3_old))-x_now)/(v[k+1]-v[k+2]))
            else:
                tk2 = float(inf)
            H1.changekey(point_k2,[k+1,tk2,t])
            p = p+1
        if k==l-2:
            # x[k+1] = x_now
            point_k2 = H1._pointerlist[k+2]  #gives the pointer to index corresponding to k+1th collision
            H1.changekey(point_k2,[k+1,float(inf),t])
        # if k == l-2:
        #     point_k2 = H1._pointerlist[k+2]  #gives the pointer to index corresponding to k+1th collision
        #     tk2_old = point_k2._k[2]
        #     if (v[k+1]-v[k+2]) > 0:
        #         print("numerator:", ((x[k+2] + v[k+2]*(t-t1))-x[k+1]))
        #         tk2 = t + (((x[k+2] + v[k+2]*(t-t1))-x[k+1])/(v[k+1]-v[k+2]))
        #     else:
        #         tk2 = float(inf)
        #     H1.changekey(point_k2,[k+1,tk2,tk2_old])


        if p ==0:
            break
        
    return(L)

    



# print(listCollisions([940.1440594570123, 342.32941684559046, 686.1000355388383, 520.8309066514597, 870.9632698994412, 727.2119773442081], [2.5912045650076445, 3.3979994719550377, 5.247957197003846, 5.383388625251065, 5.440818809376985, 6.415333653364417], [99.79672039800879, 94.19054127616612, 25.977729855078213, 25.5959601276192, 31.543951443609476, 25.267596192531126], 8, 4.531827813401554))


# r1=listCollisions([1.0, 5.0], [1.0, 2.0], [3.0, 5.0], 100, 100.0)
# # []
# r2=listCollisions([1.0, 1.0, 1.0, 1.0], [-2.0, -1.0, 1.0, 2.0], [0.0, -1.0, 1.0, 0.0], 5, 5.0)
# # [(1.0, 0, -2.0), (1.0, 2, 2.0)]
# r3=listCollisions([10000.0, 1.0, 100.0], [0.0, 1.0, 2.0], [0.0, 0.0, -1.0], 6, 10.0)
# # [(1.0, 1, 1.0), (1.505, 0, 0.0), (1.6756, 1, 0.3377), (1.7626, 0, -0.0001), (1.8163, 1, 0.2080), (1.8533, 0, -0.0002)]
# r4=listCollisions([10000.0, 1.0, 100.0], [0.0, 1.0, 2.0], [0.0, 0.0, -1.0], 100, 1.5)
# # [(1.0, 1, 1.0)]

# print(r1,r2,r3,r4,sep="\n")
# # print(r3)


# print("Starting")
# st = time.time()
# print(listCollisions([10000.0, 1.0, 100.0], [0.0, 1.0, 2.0], [0.0, 0.0, -1.0], 6, 10.0))
# [(1.0, 1, 1.0), (1.505, 0, 0.0), (1.6756, 1, 0.3377), (1.7626, 0, -0.0001), (1.8163, 1, 0.2080), (1.8533, 0, -0.0002)]
# print("Time Taken:",time.time()-st)
# print(listCollisions([940.1440594570123, 342.32941684559046, 686.1000355388383, 520.8309066514597, 870.9632698994412, 727.2119773442081], [2.5912045650076445, 3.3979994719550377, 5.247957197003846, 5.383388625251065, 5.440818809376985, 6.415333653364417], [99.79672039800879, 94.19054127616612, 25.977729855078213, 25.5959601276192, 31.543951443609476, 25.267596192531126], 8, 4.531827813401554))

# print(listCollisions([1, 1, 1, 1, 1, 1], [0, 1, 2, 3, 4, 5], [1, 0, 1, 0, 1, 0],8,100))