import sys

def printT(m, n):
    for i in range(m,n+1):
        for j in range(m,n+1):
            for k in range(m,n+1):
                sys.stdout.write("%d %d %d\n" % (i, j, k))


def printP(m, n):
    for i in range(m,n+1):
        for j in range(m,n+1):
            sys.stdout.write("%d %d %f\n" % (i, j, (1.0/(n+1-m))))

printT(0,3000)
