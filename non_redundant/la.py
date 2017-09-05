import numpy as np


def do_rotate(a, i, j, k, l, tau, s):
    g = a[i][j]
    h = a[k][l]
    a[i][j] = g - s * (h + g * tau)
    a[k][l] = h + s * (g - h * tau)


def jacobi(a, n, d, v):
    for ip in xrange(n):
        for iq in xrange(n):
            v[ip][iq] = 0.0
        v[ip][ip] = 1.0
    b = np.zeros(n)
    z = np.zeros(n)
    for ip in xrange(n):
        b[ip] = d[ip] = a[ip][ip]
        z[ip] = 0.0
    nrot = 0
    for i in xrange(50):
        sm = 0.0
        for ip in xrange(n - 1):
            for iq in xrange(ip + 1, n):
                sm += abs(a[ip][iq])
        if sm == 0.0:
            return nrot
        if i < 3:
            tresh = 0.2 * sm / (n * n)
        else:
            tresh = 0.0
        for ip in xrange(n - 1):
            for iq in xrange(ip + 1, n):
                g = 100.0 * abs(a[ip][iq])
                if i > 3 and abs(d[ip]) + g == abs(d[ip]) and abs(d[iq]) + g == abs(d[iq]):
                    a[ip][iq] = 0.0
                elif abs(a[ip][iq]) > tresh:
                    h = d[iq] - d[ip]
                    if abs(h) + g == abs(h):
                        t = a[ip][iq] / h
                    else:
                        theta = 0.5 * h / a[ip][iq]
                        t = 1.0 / (abs(theta) + np.sqrt(1.0 + theta * theta))
                        if theta < 0.0:
                            t = -t
                    c = 1.0 / np.sqrt(1 + t * t)
                    s = t * c
                    tau = s / (1.0 + c)
                    h = t * a[ip][iq]
                    z[ip] -= h
                    z[iq] += h
                    d[ip] -= h
                    d[iq] += h
                    a[ip][iq] = 0.0
                    for j in xrange(ip):
                        do_rotate(a, j, ip, j, iq, tau, s)
                    for j in xrange(ip + 1, iq):
                        do_rotate(a, ip, j, j, iq, tau, s)
                    for j in xrange(iq + 1, n):
                        do_rotate(a, ip, j, iq, j, tau, s)
                    for j in xrange(n):
                        do_rotate(v, j, ip, j, iq, tau, s)
                    nrot += 1
        for ip in xrange(n):
            b[ip] += z[ip]
            d[ip] = b[ip]
            z[ip] = 0.0
    assert False
