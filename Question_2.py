import pandas as pd
import numpy as np
from numpy.linalg import inv
import latex
from latex import start_latex,save_latex,print_objective_function,print_bnstar_matrix


A = np.asarray(pd.read_csv('csv_A.csv',header=None)).astype('float')
b = np.asarray(pd.read_csv('csv_b.csv',header=None)).astype('float')
c = np.asarray(pd.read_csv('csv_c.csv',header=None)).astype('float')
latex=start_latex()

A = np.hstack((A,np.identity(len(b))))

A = np.matrix(A)
b = np.matrix(b)
c = np.matrix(c)
latex = print_objective_function(A,b,c,latex)



print 'b =', b
print 'c =', c
print 'A =', A


def dual(A, b, c,latex):
    print "Inside Dual"


#neta and beta specification
    N = np.asarray(pd.read_csv('csv_A.csv',header=None)).astype('float')
    print 'N = ', N
    neta = [x+1 for x in xrange(N.shape[1])]
#neta = np.matrix(neta)
    print 'neta = ', neta
    latex=print_bnstar_matrix(N,latex,"N")

    B = np.identity(len(b))
    print 'B =', B
    beta = [x+neta[-1]+1 for x in xrange(B.shape[1])]
    print 'beta =', beta
    latex=print_bnstar_matrix(B,latex,"B")

#Initial values of basic variables
    Xbetastar = b
    print 'Xbetastar =', Xbetastar

#initial Non basic dual variables
    Cnetastar = c
    Znetastar = -(Cnetastar)
    print 'Znetastar =', Znetastar

#Starting Iterations
    while True:
        latex += (r"\newline \section{Iteration}")
        latex = print_bnstar_matrix(N, latex, "N")
        latex = print_bnstar_matrix(B, latex, "B")
    #Step 1:
        print 'Step 1:'

        if (Xbetastar < 0).any():
            print 'Current solution is not optimal, proceed further'
        #break
        else:
            print Xbetastar
            break
        print "The current solution is " + str(Xbetastar)
    #Step 2:
        print 'Step 2:'
        Xb_EnteringIndex = np.argmin(Xbetastar)

        i = beta[Xb_EnteringIndex]
        print i

    #Step 3:
        print 'Step 3:'
        eiSize = b.shape[0]
        print eiSize
        ei = np.zeros(eiSize)
        ei[Xb_EnteringIndex] = 1
        ei = np.matrix(ei).T
        B_inv = inv(B)
        tempdelta_zn = (np.matmul(B_inv, N))
        print tempdelta_zn
        tempdelta_Tzn = tempdelta_zn.T
        print tempdelta_Tzn
        delta_zn = np.matmul(tempdelta_Tzn, ei)
        delta_zn = -(delta_zn)
        print delta_zn
        print Znetastar
    # Step 4:
        print 'Step 4:'
        s_temp = [(delta_zn.item(x)/Znetastar.item(x)) for x in xrange(len(neta))]
        print Znetastar
        print delta_zn
        print s_temp
        s_index = s_temp.index(max(s_temp))
        smax = max(s_temp)**-1
        print smax

    # Step 5:
        print 'Step 5:'
        j = neta[s_temp.index(max(s_temp))]
        print j

    #Step 6:
        print 'Step 6:'
        ejSize = N.shape[1]
        ej = np.zeros(ejSize)
        ej[s_temp.index(max(s_temp))] = 1
        ej = np.matrix(ej)
        ej = ej.T
        B_inv = inv(B)
        tempdelta_xb = np.matmul(N, ej)
        delta_xb = np.matmul(B_inv, tempdelta_xb)
        print delta_xb

    # Step 7:
        print 'Step 7:'
        t = Xbetastar.item(Xb_EnteringIndex)/delta_xb.item(Xb_EnteringIndex)
        print t

    # Step 8:
        print 'Step 8:'
        Xbetastar = Xbetastar - (t*delta_xb)
        Znetastar = Znetastar - (smax*delta_zn)
        print Xbetastar
        print Znetastar

    # Step 9:
        print 'Step 9:'
        neta[s_index] = i
        beta[Xb_EnteringIndex] = j
        temp_matA = np.matrix(np.zeros(len(A)))
        temp_matA = temp_matA.T
        for i in range(len(beta)):
            if i == 0:
                temp_matA = temp_matA+A[:, beta[i]-1]
            else:
                temp_matA = np.append(temp_matA, A[:, int(beta[i])-1], axis=1)
        B = temp_matA

        temp_matN = np.matrix(np.zeros(len(A)))
        temp_matN = temp_matN.T
        print temp_matN

        for j in range(len(neta)):
            if j == 0:
                temp_matN = temp_matN + A[:, neta[j]-1]
            else:
                temp_matN = np.append(temp_matN, A[:, int(neta[j]) - 1], axis=1)
        N = temp_matN

        print B
        print N
        print "beta"+str(beta)
        print neta

    #new basic primal variables and nonbasic dual variables
        print Znetastar
        for p in range(len(Xbetastar)):
            if Xbetastar[p] == 0:
                Xbetastar[p] = t
        for q in range(len(Znetastar)):
            if Znetastar[q] == 0:
                Znetastar[q] = smax
                print Znetastar
                print Xbetastar
        latex = print_bnstar_matrix(np.matrix(neta), latex, "N indices")
        latex = print_bnstar_matrix(np.matrix(beta), latex, "B indices")
        latex = print_bnstar_matrix(Znetastar, latex, "z_n")
        latex = print_bnstar_matrix(Xbetastar, latex, "x_b")
    #Znetastar_min = np.amin(Znetastar)
        if (Xbetastar > 0).all():
            print 'Since Znetastar has all non-negative components, the current solution is optimal'
            print Xbetastar
            break
    return latex



latex=dual(A,b,c,latex)
latex+=r"\end{document}"

save_latex("file.tex",latex)



