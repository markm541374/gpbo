"""from robo.fmin import fabolas_fmin
import scipy as sp
import re

def runfab(Xl,Xu,Sworst,Sbest,f,n,D,ofpath='fabresults.csv'):


    sw=sp.log(100)
    sb=sp.log(10000)
    X_lower = sp.array([0.]*D+[sw])
    X_upper = sp.array([1.]*D+[sb])
    def fw(x,s):
        strue = Sworst + (Sbest-Sworst)*(sp.exp(s[0])-100)/9900
        xtrue = sp.empty(D)
        for i in xrange(D):
            xtrue[i] = Xl[i]+(Xu[i]-Xl[i])*x[0,i]
        print 'calling f with x:{} s:{} from fab => x:{} s:{} passed to ojf'.format(x,s,xtrue,strue)
        y,c = f(xtrue,strue)
        print 'ojf returned y:{} c:{}'.format(y,c)
        return sp.array([[y]]), sp.array([[c]])

    if True:
        O1 = fabolas_fmin(fw, X_lower, X_upper, num_iterations=max(0,(n-40)),rec='change2')
    augres(ofpath,Xl,Xu,Sbest,f)
    return

def augres(outfile,xl,xu,sbest,f):

    out = open(outfile,'w')
    text = open('results.csv').read()
    pattern = re.compile(r'".*?"', re.DOTALL)
    textnonl=pattern.sub(lambda x: x.group().replace('\n', ''), text)
    print textnonl
    lines = re.findall('(?:"[^"]*"|.)+', textnonl)

    for i,line in enumerate(lines):
        print line
        if i==0:
            out.write(line.strip('\r\n')+',xtrue,strue,c,incumbenttrue,trueyatxrecc\n')
        else:
            vals = line.split(',')
            con = [float(j) for j in vals[1].strip('[]').split()]
            xtrue = [xl[i]+(xu[i]-xl[i])*con[i] for i in xrange(len(con[:-1]))]
            strue = con[-1]
            ctrue = sp.exp(float(vals[10].strip('\n\r[]')))
            inc = [float(j) for j in vals[3].strip('[]').split()]
            incumbenttrue = [xl[i]+(xu[i]-xl[i])*inc[i] for i in xrange(len(con[:-1]))]
            yatxi = f(incumbenttrue,sbest)[0]
            out.write(line.strip('\r\n')+','+str(xtrue).replace(',',' ')+','+str(strue)+','+str(ctrue)+','+str(incumbenttrue).replace(',',' ')+','+str(yatxi)+'\n')
    out.close()
"""