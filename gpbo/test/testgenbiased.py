from gpbo.core import objectives

D=3
lb=[-1]*D
ub=[1]*D
q=[]
for i in range(10):
    q.append(objectives.genbiasedmat52ojf(D,lb,ub,0.5))

print [str([opt[1],opt[2]])+'\n' for opt in q]