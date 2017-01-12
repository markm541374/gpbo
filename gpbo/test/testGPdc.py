#!/usr/bin/env python2
#encoding: UTF-8

# demo the derivatives for the GP lib
def main():
    import scipy as sp
    from scipy import linalg as spl
    from matplotlib import pyplot as plt

    from gpbo.core import GPdc

    f,a = plt.subplots(3)
    ns = 200
    sup = sp.linspace(-1,1,ns)

    #points are a high noise obs, a low noise obs, a derivative obs and a second derivative obs
#    X = sp.array([[-0.8],[-0.25],[0.25],[0.8]])
#    Y = sp.array([[0.3],[-0.2],[2.5],[50.]])
#    S = sp.array([[1e-1],[1e-6],[1e-6],[1e-6]])
#    D = [[sp.NaN],[sp.NaN],[0],[0,0]]
    nn=30
    X = sp.array([sp.linspace(-1,1,nn)]).T
    Y = sp.array([map(lambda x:sp.sin(6*x), sp.linspace(-1,1,nn))]).T
    S = sp.array([[1e-1]*nn]).T
    D = [[sp.NaN]]*nn
    print X.shape
    print Y.shape
    print S.shape
    print len(D)
    #a[0].plot([-0.8,-0.8],[0.3-sp.sqrt(1e-1),0.3+sp.sqrt(1e-1)],'r')
    #a[0].plot([-0.8],[0.3],'ro')
    #a[0].plot([-0.25],[-0.2],'ro')
    #a[1].plot([0.25],[2.5],'ro')
    #a[2].plot([0.8],[50],'ro')
    for i in range(nn):
        a[0].plot(X[i,0],Y[i,0],'ro')
    k= GPdc.kernel(GPdc.SQUEXP, 1, sp.array([0.5, 0.2]))
    K = sp.empty([4,4])
    for i in xrange(4):
        for j in xrange(i,4):
            K[i,j] =K[j,i] = k(X[i,:],X[j,:],D[i],D[j])
        K[i,i]+=1e-6


    g = GPdc.GPcore(X, Y, S, D, GPdc.kernel(GPdc.SQUEXP, 1, sp.array([0.5, 0.2])))
    #g.printc()
    C= g.get_cho()
    #print C

    #print spl.cho_solve((C,True),sp.eye(4)).dot(K)
    for i,d in enumerate([[sp.NaN],[0],[0,0]]):
        m,v = g.infer_diag(sup,[d]*ns)

        null,V = g.infer_full_post(sup,[d]*ns)
        Vu = (m+sp.sqrt(sp.diag(V))).flatten()
        a[i].plot(sup, Vu.flatten(),'g')

        vl = (m-sp.sqrt(v)).flatten()
        vu = (m+sp.sqrt(v)).flatten()
        a[i].plot(sup,m.flatten())
        a[i].fill_between(sup,vl,vu,facecolor='LightBlue',edgecolor='LightBlue')
        a[i].set_ylabel(''.join(['d/dx ']*i)+'f')

        dr = g.draw_post(sup,[d]*ns,5)

        for j in range(5):
            a[i].plot(sup,dr[j,:],'lightpink')
    plt.show()
    return

if __name__=="__main__":
    main()
