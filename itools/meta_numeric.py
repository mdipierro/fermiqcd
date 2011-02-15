from math import *
from copy import copy, deepcopy
from random import *
import time

def Matrix(rows=0,cols=0):
    """Constuctor for a Matrix [[1,0],[0,1]]"""
    return [[0]*cols for k in range(rows)]

def pprint(A):
    print ' ['
    for line in A:
        print '  [',
        for item in line:
            print '%.3e,' % float(item),
            pass
        print '],'
        pass
    print ' ]'
    return
    
def rows(A):
    """Rows of a Matrix"""
    return len(A)

def cols(A):
    """Cols of a Matrix"""
    return len(A[0])

def add(A,B):
    """multiplies a number of Matrix A by a Matrix B"""
    if type(A)==type(1) or type(A)==type(1.0):
        C=deepcopy(B)
        for i in range(rows(B)):            
            C[i][i]+=A
            pass
        return C
    else:
        C=deepcopy(B)
        for i in range(rows(B)):            
            for j in range(cols(B)):            
                C[i][j]+=A[i][j]
                pass
            pass
        return C
    pass

def sub(A,B):
    """multiplies a number of Matrix A by a Matrix B"""
    C=deepcopy(A)
    for i in range(rows(A)):            
        for j in range(cols(A)):            
            C[i][j]-=B[i][j]
            pass
        pass
    return C


def multiply(A,B):
    """multiplies a number of Matrix A by a Matrix B"""
    if type(A)==type(1) or type(A)==type(1.0):
        C=Matrix(rows(B),cols(B))        
        for i in range(rows(B)):
            for j in range(cols(B)):
                C[i][j]=A*B[i][j]
                pass
            pass
        return C
    else:
        C=Matrix(rows(A),cols(B))        
        for i in range(rows(A)):
            for j in range(cols(B)):
                for k in range(cols(A)):
                    C[i][j]+=A[i][k]*B[k][j]
                    pass
                pass
            pass        
        return C
    pass

def inverse(A,checkpoint=None):
    """Computes the inverse of A using Gauss-Jordan emilimination"""
    A=deepcopy(A)
    n=rows(A)
    B=Matrix(n,n)
    for i in range(n): B[i][i]=1
    for c in range(n):
        if checkpoint: checkpoint('pivoting (%i) ...' % c)
        for r in range(c+1,n):
            if abs(A[r][c])>abs(A[c][c]):
                A[r],A[c],B[c],B[r]=A[c],A[r],B[r],B[c]
                pass
            pass
        p=float(A[c][c])
        for k in range(n):
            A[c][k],B[c][k]=float(A[c][k])/p,float(B[c][k])/p
            pass        
        for r in range(0,c)+range(c+1,n):
            p=float(A[r][c])
            for k in range(n):
                A[r][k]-=p*A[c][k]
                B[r][k]-=p*B[c][k]
                pass
            pass
        pass
    return B

def test_inverse():
    """Test for the inverse(A) function"""
    print "\n\nTesting inverse(A)........."
    A=[[1,2,3],[2,4,8],[1,3,7]]
    print "A=",A
    B=inverse(A)
    print "B=",B
    C=multiply(A,B)
    print "A*B=",C
    return

def transpose(A):
    """Transposed of A"""
    B=Matrix(cols(A),rows(A))
    for i in range(rows(B)):
        for j in range(cols(B)):
            B[i][j]=A[j][i]
            pass
        pass
    return B

def Cholesky(A):
    if A!=transpose(A): raise SyntaxError
    L=deepcopy(A)
    for k in range(cols(L)):
        if L[k][k]<=0: raise SyntaxError
        p=L[k][k]=sqrt(L[k][k])
        for i in range(k+1,rows(L)):        
            L[i][k]/=p
            pass
        for j in range(k+1,rows(L)):
            p=float(L[j][k])
            for i in range(k+1,rows(L)):
                L[i][j]-=p*L[i][k]
                pass
            pass
        pass
    for  i in range(rows(L)):
        for j in range(i+1,cols(L)):
            L[i][j]=0
            pass
        pass
    return L

def test_Cholesky():
    """Test for the inverse(A) function"""
    print "\n\nTesting Cholesky(A)........."
    A=[[4,2,1],[2,9,3],[1,3,16]]
    print "A=",A
    L=Cholesky(A)
    print "L=",L
    C=sub(multiply(L,transpose(L)),A)
    print "L*L^T-A=",C
    return

def identity(n):
    A=Matrix(n,n)
    for i in range(n): A[i][i]=1
    return A

def diagonal(v):
    n=len(v)
    A=Matrix(n,n)
    for i in range(n): A[i][i]=v[i]
    return A

def maxind(S,k):
    j=k+1
    for i in range(k+2,len(S[k])):
        if abs(S[k][i])>abs(S[k][j]): j=i
        pass
    return j

def Jacobi(A,checkpoint):
    """Returns U end e so that A=U*diagonal(e)*transposed(U)
       where i-column of U contains the eigenvector corresponding to
       the eigenvalue e[i] of A.

       from http://en.wikipedia.org/wiki/Jacobi_eigenvalue_algorithm
    """
    t0=time.time()
    n=rows(A)
    if n!=cols(A): raise SyntaxError
    S=Matrix(n,n)
    for i in range(n):
        for j in range(n):
            S[i][j]=float(A[i][j])
            pass
        pass
    E=identity(n)
    state=n
    ind=[maxind(S,k) for k in range(n)]
    e=[S[k][k] for k in range(n)]
    changed=[True for k in range(n)]
    iteration=0
    while state:
        if checkpoint: checkpoint('rotating vectors (%i) ...' % iteration)
        m=0
        for k in range(1,n-1):
            if abs(S[k][ind[k]])>abs(S[m][ind[m]]): m=k
            pass
        k,h=m,ind[m]
        p=S[k][h]
        y=(e[h]-e[k])/2
        t=abs(y)+sqrt(p*p+y*y)
        s=sqrt(p*p+t*t)
        c=t/s
        s=p/s
        t=p*p/t
        if y<0: s,t=-s,-t
        S[k][h]=0
        #update(k,-t)
        y=e[k]
        e[k]=y-t
        if changed[k] and y==e[k]: changed[k],state=False,state-1
        elif (not changed[k]) and y!=e[k]: changed[k],state=True,state+1
        #update(h,t)
        y=e[h]
        e[h]=y+t
        if changed[h] and y==e[h]: changed[h],state=False,state-1
        elif (not changed[h]) and y!=e[h]: changed[h],state=True,state+1

        for i in range(k):
            S[i][k],S[i][h]=c*S[i][k]-s*S[i][h],s*S[i][k]+c*S[i][h]
            pass
        for i in range(k+1,h):
            S[k][i],S[i][h]=c*S[k][i]-s*S[i][h],s*S[k][i]+c*S[i][h]
            pass
        for i in range(h+1,n):
            S[k][i],S[h][i]=c*S[k][i]-s*S[h][i],s*S[k][i]+c*S[h][i]
            pass
        for i in range(n):
            E[k][i],E[h][i]=c*E[k][i]-s*E[h][i],s*E[k][i]+c*E[h][i]
            pass
        ind[k],ind[h]=maxind(S,k),maxind(S,h)
        iteration+=1
        pass
    # SORT VECTORS
    for i in range(1,n):
        j=i
        while j>0 and e[j-1]>e[j]:
            e[j],e[j-1]=e[j-1],e[j]
            E[j],E[j-1]=E[j-1],E[j]
            j-=1
            pass
        pass
    # NORMALIZE VECTORS
    U=Matrix(n,n)
    for i in range(n):
        sum=0.0
        for j in range(n): sum+=E[i][j]**2;            
        sum=sqrt(sum)
        for j in range(n): U[j][i]=E[i][j]/sum
        pass
    return U,e    

def test_Jacobi():
    """Test the Jacobi algorithm"""
    print "Testing Jacobi on random matrices..."
    n=4
    A=Matrix(n,n)
    for k in range(3):
        for i in range(n):
            for j in range(i,n):
                A[i][j]=A[j][i]=gauss(10,10)
                pass
            pass
        print "A=",
        pprint(A)
        U,e=Jacobi(A)
        print "U*e*U^T-A=",
        pprint(sub(multiply(U,multiply(diagonal(e),transpose(U))),A))
        pass
    return

def fitting_function(f,C,x):
    sum=0.0
    i=0
    for func in f:
        sum+=func(x)*C[i][0]
        i+=1
        pass
    return sum
    

def fit(f,x,y,dy=None):
    """Linear fit of y[i]+/idy[i] using sum_j f[j](x[i])"""
    A=Matrix(len(x),len(f))
    B=Matrix(len(y),1)
    for i in range(rows(A)):
        w=1.0
        if dy: w=1.0/sqrt(dy[i])
        B[i][0]=w*y[i]
        for j in range(cols(A)):
            A[i][j]=w*f[j](x[i])
            pass
        pass
    C=multiply(inverse(multiply(transpose(A),A)),multiply(transpose(A),B))
    chi=sub(multiply(A,C),B)
    chi2=multiply(transpose(chi),chi)
    ff=lambda x,C=C,f=f,q=fitting_function: q(f,C,x)
    return C,chi2[0][0],ff

"""Sets of fitting functions for fit"""
CONSTANT=[1]
LINEAR=[lambda x: 1.0, lambda x: x]
QUADRATIC=[lambda x: 1.0, lambda x: x, lambda x: x*x]
CUBIC=[lambda x: 1.0, lambda x: x, lambda x: x*x, lambda x: x*x*x]
QUARTIC=[lambda x: 1.0, lambda x: x, lambda x: x*x, lambda x: x*x*x, lambda x: x*x*x*x]
def POLYNOMIAL(n):
    """Generic polynmial fitting function"""
    return [eval('lambda x: pow(x,%i)' % i) for i in range(n+1)]        
def EXPONENTIAL(n):
    """Generic exponential fitting function"""
    return [eval('lambda x: exp(x*%i)' % i) for i in range(n+1)]        

def test_fit():
    """Test for the fit function"""
    print "\n\nTesting fit(QUADRATIC,...)........."
    print "data generated using 5+0.8*k+0.3*k*k+gauss(0,1)"
    x=[]
    y=[]
    for k in range(100):
        x.append(k)
        y.append(5+0.8*k+0.3*k*k+gauss(0,1))
        pass
    a,chi2,ff=fit(QUADRATIC,x,y)
    print "f(x)=(",a[0][0],")+(",a[1][0],")*x+(",a[2][0],")*x*x"
    print "chi2=",chi2
    
def AR1filter(r):
    r.sort()    
    """Performs AR(1) filtering and eliminates auto-correlation"""
    x=[]
    y=[]
    for i in range(1,len(r)):
        x.append(r[i-1][1])
        y.append(r[i][1])
        pass
    try:
        a,chi2,ff=fit(LINEAR,x,y)
        a1=a[1][0]
    except:
        a=[]
        a1=0
    r2=[r[0]]
    for i in range(1,len(r)):
        r2.append((r[i][0],(r[i][1]-a1*r[i-1][1])/(1.0-a1)))
        pass
    r2.reverse()
    return r2,a

def test_AR1filter():
    """Test for AR1filter"""
    print "\n\nTesting AR1filter()........."
    y=[(0,0.05)]
    for k in range(1,30):
        y.append((k,0.05+y[k-1][1]*0.3+gauss(0,0.01)))
        pass
    yf,a=AR1filter(y)
    for k in range(30):
        print y[k],',',yf[k]
        pass
    return

def cov2cor(cov):
    n=rows(cov)
    sigma=[0]*n
    cor=Matrix(n,n)
    for i in range(n):
        sigma[i]=sqrt(cov[i][i])
        for j in range(0,i+1):
            cor[i][j]=cor[j][i]=cov[i][j]/sigma[i]/sigma[j]            
            pass
        pass
    return cor, sigma

def cor2cov(cor,sigma):
    n=rows(cor)
    cov=Matrix(n,n)    
    for i in range(n):
        for j in range(0,i+1):
            cov[i][j]=cov[j][i]=cor[i][j]*sigma[i]*sigma[j]
            pass
        pass
    return cov

def cor2cor(cor):
    """deprecated!"""
    n=rows(cor)
    for i in range(n):
        for j in range(0,i):
            a=cor[i][j]
            b=abs(a)
            if b>1: cor[i][j]=cor[j][i]=a/b
            pass
        pass
    return cor

def truncate_eigenvalues(A,delta=0.01,checkpoint=None):
    """Takes a symmetric matrix and relaces all eigenvalues with < delta with delta"""
    U,e1=Jacobi(A,checkpoint)
    e2=deepcopy(e1)
    for i in range(len(e2)):
        if e2[i]<delta: e2[i]=delta
        pass
    return multiply(U,multiply(diagonal(e2),transpose(U))),e1,e2

def truncate_eigenvalues_cor(cor,delta=0.01,checkpoint=None):
    """like truncate_eigenvalues but restores the diagonal elements to 1"""
    if len(cor)<2: return cor,[1],[1]
    cov2,e1,e2=truncate_eigenvalues(cor,delta,checkpoint)
    cor,sigma2=cov2cor(cov2) # restore 1 on diagonal
    return cor,e1,e2

def truncate_eigenvalues_cov(cov,delta=0.01,checkpoint=None):
    """projects into a valid correlation matrix"""
    if len(cov)<2: return cov,[],[]
    cor,sigma=cov2cor(cov)
    cor,e1,e2=truncate_eigenvalues_cor(cor,delta,checkpoint)
    cov=cor2cov(cor,sigma)
    return cov,e1,e2

def test_truncate_eigenvalues_cov():
    print "Testing truncate_eigenvalues_cov"
    n=5
    A=Matrix(n,n)    
    for i in range(n):
        A[i][i]=abs(gauss(10,10))
        for j in range(i+1,n):
            A[i][j]=A[j][i]=gauss(10,10)
            pass
        pass
    print 'A=',
    pprint(A)
    B,e1,e2=truncate_eigenvalues_cov(A)
    print 'B=',
    pprint(B)
    for i in range(n):
        print e1[i],e2[i]
    return

def CorrelationMatrix(table,do_AR1fiter=True):    
    """Takes a table of the form [[(date,return),...],[(date,return),...],...]
       and filters each time series (AR1fiter). computes covariance matrix,
       and removes negative eigenvalues"""    
    n=len(table)
    mean=[0]*n
    for series in table:
        r=[r for d,r in series]
        if do_AR1filter: a, r=AR1filter(r)
        for i in range(len(series)):
            series[i]=(series[i][d],r[i])
            mean[i]+=r[i]
            pass
        mean[i]/=len(series)
        pass

def normalize(v):
    v=deepcopy(v)
    norm=0.0
    for x in v: norm+=x[0]
    for x in v: x[0]=float(x[0])/abs(norm)
    return v

def risk_return(cov,returns,x_mar,r_free):
    n=len(returns)
    r_freev=[[r_free]]*n
    x_free=1.0-sum([x[0] for x in x_mar])
    r=multiply(transpose(returns),x_mar)[0][0]+x_free*r_free
    sigma=sqrt(multiply(transpose(x_mar),multiply(cov,x_mar))[0][0])
    return x_free, r, sigma

def Markowitz(returns, cov, r_free, checkpoint=None):    
    n=len(returns)    
    inv_cov=inverse(cov,checkpoint)
    r_freev=[[r_free]]*n
    ones=[[1.0]]*n
    returns=[[r] for r in returns]
    invSI=multiply(inv_cov,ones)
    IinvSI=multiply(transpose(ones),invSI)[0][0]
    invSr=multiply(inv_cov,returns)   
    IinvSr=multiply(transpose(ones),invSr)[0][0]
    rinvSr=multiply(transpose(returns),invSr)[0][0]
    # portfolio min variance
    x_min=normalize(invSI)
    r_min=multiply(transpose(returns),x_min)[0][0]
    sigma_min=sqrt(multiply(transpose(x_min),multiply(cov,x_min))[0][0])
     # portfolio Markowitz
    invSr=multiply(inv_cov,sub(returns,r_freev))
    x_mar=normalize(invSr)
    x_free,r_mar,sigma_mar=risk_return(cov,returns,x_mar,r_free)
    x_mar=[x[0] for x in x_mar]
    # parameters of hyperbola
    a2=1.0/IinvSI
    b2=(rinvSr*IinvSI-IinvSr**2)*a2*a2
    # hyperbola sigma^2/a2 - (r-r_min)^2/b2=1
    return x_mar, r_mar, sigma_mar, (r_min, a2, b2)
    

def mean(series):
    """series=[(datetime, rate),..]
    retuns mean and points (second col of series)"""
    sum=0.0
    rs=[]
    for d,r in series:
        rs.append(r)
        sum+=r
    return sum/len(series), rs
        
def covariance(series1,m1,series2,m2):
    """series1=[(datetime, rate1),..], m1 is mean of rates1
    series2=[(datetime, rate2),..], m1 is mean of rates2
    retuns covariance and points"""
    sum=0.0
    rs=[]
    i=j=k=0
    while i<len(series1) and j<len(series2):
        d1,r1=series1[i]
        d2,r2=series2[j]
        if d1<d2: j+=1
        elif d1>d2: i+=1
        else:
            x=(r1-m1)*(r2-m2)
            sum+=x
            rs.append(x)
            i,j,k=i+1,j+1,k+1
            pass
        pass
    if k>0:
        return sum/k, rs
    return 0.0,[]

def test_all():
    test_inverse()
    test_Cholesky()
    test_fit()
    test_AR1filter()
    test_Jacobi()
    test_truncate_eigenvalues_cov()

if __name__=='__main__': test_all()
