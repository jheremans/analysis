f=np.linspace(3.6,4.5,10)
deg=np.array([180.,60.,300.,420.])

matrix=np.zeros([5 , 11])

matrix[0,1:]=f[:]
matrix[1:,0]=deg[:]


for i in range(5)[1:5]:
    for j in range(11)[1:11]:
        matrix[i,j]=int(round(deg[i-1]/180*1/(2*f[j-1])*1e3,0)+11)
    
print matrix 