import numpy as np

u1=[0,2,2,2,2,1,1,2,2,2,3,2,2,0,0,0,2,2,2,2,2]
u2=[0,0,0,0,0,1,1,1,1,1,1,2,2,3,3,3,3,3,3,3,3]
v1=[0,1,1,1,3,3,3,3,3,0,0,0,1,1,1,2,2,2,2,3,3]
v2=[1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3]

def to_float(a, it_max=3):
	temp=[]
	for i in a:
		if type(i)!=list:
			temp.append(float(i))
		elif it_max>0:
			temp.append(i,it_max=itmax-1)
	return temp

u1=to_float(u1)
u2=to_float(u2)
v1=to_float(v1)
v2=to_float(v2)

a1=np.array([1,0,4,0,0,2,3,1,0,0,2,0,3,0,5,0])
a2=np.array([0,0,0,0,1,2,0,0,0,1,0,5,3,3,4,2])

a1=np.array(to_float(a1))/sum(a1)
a2=np.array(to_float(a2))/sum(a2)
m=np.array([[np.sqrt((i//4-j//4)**2+(i%4-j%4)**2)/4 for i in range(16)] for j in range(16)])

