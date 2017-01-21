import pandas as pd
import numpy as np
import Question_1
import Question_2
import latex
from latex import start_latex,save_latex,print_objective_function
from Question_1 import primal
from Question_2 import dual
from numpy.linalg import inv

def solver(A,b,c):
    bflag=True
    for i in range(len(b)):
        if b[i] < 0:
            bflag = False
    if bflag == True:
        primal(A, b, c)
    else:
        cflag = False
        for i in range(len(c)):
            if c[i] > 0:
                cflag = True
        if cflag == True:
            #phaseOne(A, b, c)
             crissCross(A,b,c)
        else:
            dual(A, b, c)

def phaseOne(A, b, c):
    print 'dual is called'
    # neta and beta specification
    N = np.asarray(pd.read_csv('csv_A.csv', header=None)).astype('float')
    print 'N = ', N
    neta = [x + 1 for x in xrange(N.shape[1])]
    print 'neta = ', neta

    B = np.identity(len(b))
    print 'B =', B
    beta = [x + neta[-1] + 1 for x in xrange(B.shape[1])]
    # beta = np.matrix(beta)
    print 'beta =', beta

    # Initial values of basic variables
    Xbetastar = b
    print 'Xbetastar =', Xbetastar

    # initial Non basic dual variables
    Cnetastar = c
    Znetastar = -(Cnetastar)
    print 'Znetastar =', Znetastar

    temp_cn = np.array(Cnetastar)
    temp_NB = np.array(neta)

    #making Cnetastar to an array of -1
    for i in range(len(Cnetastar)):
        Cnetastar[i] = -1
    Znetastar = -(Cnetastar)
    while True:
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
        print 'i =', i

    #Step 3:
        print 'Step 3:'
        eiSize = b.shape[0]
        print eiSize
        ei = np.zeros(eiSize)
        ei[Xb_EnteringIndex] = 1
        ei = np.matrix(ei).T
        B_inv = inv(B)
        tempdelta_zn = (np.matmul(B_inv, N))
        print 'tempdelta_zn =', tempdelta_zn
        tempdelta_Tzn = tempdelta_zn.T
        print 'tempdelta_Tzn =', tempdelta_Tzn
        delta_zn = np.matmul(tempdelta_Tzn, ei)
        delta_zn = -(delta_zn)
        print 'delta_zn =', delta_zn
        print 'Znetastar =', Znetastar
    # Step 4:
        print 'Step 4:'
        s_temp = [(delta_zn.item(x)/Znetastar.item(x)) for x in xrange(len(neta))]
        print 'Znetastar =', Znetastar
        print 'delta_zn =', delta_zn
        print 's_temp =', s_temp
        s_index = s_temp.index(max(s_temp))
        smax = max(s_temp)**-1
        print 'smax =', smax

    # Step 5:
        print 'Step 5:'
        j = neta[s_temp.index(max(s_temp))]
        print 'j =', j

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
        print 'delta_xb =', delta_xb

    # Step 7:
        print 'Step 7:'
        t = Xbetastar.item(Xb_EnteringIndex)/delta_xb.item(Xb_EnteringIndex)
        print 't =', t

    # Step 8:
        print 'Step 8:'
        Xbetastar = Xbetastar - (t*delta_xb)
        Znetastar = Znetastar - (smax*delta_zn)
        print 'Xbetastar =', Xbetastar
        print 'Znetastar =', Znetastar

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
        print 'temp_matN =', temp_matN

        for j in range(len(neta)):
            if j == 0:
                temp_matN = temp_matN + A[:, neta[j]-1]
            else:
                temp_matN = np.append(temp_matN, A[:, int(neta[j]) - 1], axis=1)
        N = temp_matN

        print B
        print N
        print "beta"+str(beta)
        print 'neta =', neta

    #new basic primal variables and nonbasic dual variables
        print Znetastar
        for p in range(len(Xbetastar)):
            if Xbetastar[p] == 0:
                Xbetastar[p] = t
        for q in range(len(Znetastar)):
            if Znetastar[q] == 0:
                Znetastar[q] = smax
                print 'Znetastar =', Znetastar
                print 'Xbetastar =', Xbetastar

        if (Xbetastar > 0).all():
            print 'Since Znetastar has all non-negative components, the current solution is optimal'
            print 'Xbetastar =', Xbetastar
            break

#Calling primal
    print 'primal is called'
    cb = np.matrix(np.zeros(len(beta)))
    cb = cb.T


    Cnetastar = np.matrix(np.zeros(len(neta)))
    Cnetastar = Cnetastar.T

    # Setting the CB Matrix

    for i in range(len(beta)):

        var = beta[i]
        print 'var =', var

        if temp_NB.__contains__(var):
            index = int(np.where(temp_NB == var)[0])
            cb[i] = temp_cn[index]
    print 'cb =', cb

    for j in range(len(neta)):

        var = neta[j]
        print 'var =', var

        if temp_NB.__contains__(var):
            index = int(np.where(temp_NB == var)[0])
            Cnetastar[j] = temp_cn[index]
    print 'Cnetastar =', Cnetastar

    # coefficients of objective function
    B_inv = inv(B)
    B_invN = np.matmul(B_inv, N)
    B_invNT = B_invN.T
    B_invN_cb = np.matmul(B_invNT, cb)
    temp_cb = -(B_invN_cb - Cnetastar)
    print 'temp_cb =', temp_cb
    Znetastar = -temp_cb

    while True:
    #Step 1:
        print 'Step 1:'
        if (Znetastar > 0).all():
            print 'Current solution is Optimal'
            break
        else:
            print 'Current solution is not optimal, proceed further'

    #Step 2:
        print 'Step 2:'
        Zns_EnteringIndex = np.argmin(Znetastar)
        j = neta[Zns_EnteringIndex]
        print 'j =', j

    #Step 3:
        print 'Step 3:'
        ejSize = N.shape[1]
        ej = np.zeros(ejSize)
        ej[Zns_EnteringIndex] = 1
        ej = np.matrix(ej)
        ej = ej.T
        B_inv = inv(B)
        tempdelta_xb = np.matmul(N, ej)
        delta_xb = np.matmul(B_inv, tempdelta_xb)
        print 'delta_xb =', delta_xb

    # Step 4:
        print 'Step 4:'
        SizeOfb = b.shape[0]
        t_temp = [(delta_xb.item(x)/Xbetastar.item(x)) for x in xrange(SizeOfb)]
        t_index = t_temp.index(max(t_temp))
        tmax = max(t_temp)**-1
        print 'tmax =', tmax

    # Step 5:
        print 'Step 5:'
        i = beta[t_temp.index(max(t_temp))]
    #i = beta[i]
        print 'i =', i

    #Step 6:
        print 'Step 6:'
        eiSize = b.shape[0]
        ei = np.zeros(eiSize)
        ei[t_index] = 1
        ei = np.matrix(ei).T
        tempdelta_zn = (np.matmul(B_inv, N))
        tempdelta_Tzn = tempdelta_zn.T
        delta_zn = np.matmul(tempdelta_Tzn, ei)
        delta_zn = -(delta_zn)
        print 'delta_zn =', delta_zn

    # Step 7:
        print 'Step 7:'
        s = Znetastar.item(Zns_EnteringIndex)/delta_zn.item(Zns_EnteringIndex)
        print 's =', s
        print 'Znetastar =', Znetastar

    # Step 8:
        print 'Step 8:'
        Xbetastar = Xbetastar - (tmax*delta_xb)
        Znetastar = Znetastar - (s*delta_zn)
        print 'Xbetastar =', Xbetastar
        print 'Znetastar =', Znetastar

    # Step 9:
        print 'Step 9:'
        beta[t_index] = j
        neta[Zns_EnteringIndex] = i
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
        print 'temp_matN =', temp_matN

        for j in range(len(neta)):
            if j == 0:
                temp_matN = temp_matN + A[:, neta[j]-1]
            else:
                temp_matN = np.append(temp_matN, A[:, int(neta[j]) - 1], axis=1)
        N = temp_matN

        print 'B =', B
        print 'N =', N
        print 'beta =', beta
        print 'neta =', neta

    #new basic primal variables and nonbasic dual variables
        print 'Znetastar =', Znetastar
        for p in range(len(Xbetastar)):
            if Xbetastar[p] == 0:
                Xbetastar[p] = tmax
        for q in range(len(Znetastar)):
            if Znetastar[q] == 0:
                Znetastar[q] = s
                print 'Znetastar =', Znetastar
                print 'Xbetastar =', Xbetastar

    #Znetastar_min = np.amin(Znetastar)
        if (Znetastar > 0).all():
            print 'Since Znetastar has all non-negative components, the current solution is optimal'
            break

def primal(A, b, c):


#neta and beta specification
    N = np.asarray(pd.read_csv('csv_A.csv',header=None)).astype('float')
    print 'N = ', N
    neta = [x+1 for x in xrange(N.shape[1])]
    print 'neta = ', neta

    B = np.identity(len(b))
    print 'B =', B
    beta = [x+neta[-1]+1 for x in xrange(B.shape[1])]
    print 'beta =', beta

#Initial values of basic variables
    Xbetastar = b
    print 'Xbetastar =', Xbetastar

#initial Non basic dual variables
    Cnetastar = c
    Znetastar = -(Cnetastar)
    print 'Znetastar =', Znetastar

#Starting Iterations
    while True:
    #Step 1:
        print 'Step 1:'
        if (Znetastar > 0).all():
            print 'Current solution is Optimal'
            break
        else:
            print 'Current solution is not optimal, proceed further'

    #Step 2:
        print 'Step 2:'
        Zns_EnteringIndex = np.argmin(Znetastar)
        j = neta[Zns_EnteringIndex]
        print j

    #Step 3:
        print 'Step 3:'
        ejSize = N.shape[1]
        ej = np.zeros(ejSize)
        ej[Zns_EnteringIndex] = 1
        ej = np.matrix(ej)
        ej = ej.T
        B_inv = inv(B)
        tempdelta_xb = np.matmul(N, ej)
        delta_xb = np.matmul(B_inv, tempdelta_xb)
        print delta_xb

    # Step 4:
        print 'Step 4:'
        SizeOfb = b.shape[0]
        t_temp = [(delta_xb.item(x)/Xbetastar.item(x)) for x in xrange(SizeOfb)]
        t_index = t_temp.index(max(t_temp))
        tmax = max(t_temp)**-1
        print tmax

    # Step 5:
        print 'Step 5:'
        i = beta[t_temp.index(max(t_temp))]
    #i = beta[i]
        print i

    #Step 6:
        print 'Step 6:'
        eiSize = b.shape[0]
        ei = np.zeros(eiSize)
        ei[t_index] = 1
        ei = np.matrix(ei).T
        tempdelta_zn = (np.matmul(B_inv, N))
        tempdelta_Tzn = tempdelta_zn.T
        delta_zn = np.matmul(tempdelta_Tzn, ei)
        delta_zn = -(delta_zn)
        print delta_zn

    # Step 7:
        print 'Step 7:'
        s = Znetastar.item(Zns_EnteringIndex)/delta_zn.item(Zns_EnteringIndex)
        print s
        print Znetastar

    # Step 8:
        print 'Step 8:'
        Xbetastar = Xbetastar - (tmax*delta_xb)
        Znetastar = Znetastar - (s*delta_zn)
        print Xbetastar
        print Znetastar

    # Step 9:
        print 'Step 9:'
        beta[t_index] = j
        neta[Zns_EnteringIndex] = i
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

        print 'B =', B
        print 'N =', N
        print 'beta =', beta
        print 'neta =', neta

    #new basic primal variables and nonbasic dual variables
        print Znetastar
        for p in range(len(Xbetastar)):
            if Xbetastar[p] == 0:
                Xbetastar[p] = tmax
        for q in range(len(Znetastar)):
            if Znetastar[q] == 0:
                Znetastar[q] = s
                print Znetastar
                print Xbetastar

    #Znetastar_min = np.amin(Znetastar)
        if (Znetastar > 0).all():
            print 'Since Znetastar has all non-negative components, the current solution is optimal'
            break

def dual(A, b, c):

#neta and beta specification
    N = np.asarray(pd.read_csv('csv_A.csv',header=None)).astype('float')
    print 'N = ', N
    neta = [x+1 for x in xrange(N.shape[1])]
#neta = np.matrix(neta)
    print 'neta = ', neta

    B = np.identity(len(b))
    print 'B =', B
    beta = [x+neta[-1]+1 for x in xrange(B.shape[1])]
    print 'beta =', beta

#Initial values of basic variables
    Xbetastar = b
    print 'Xbetastar =', Xbetastar

#initial Non basic dual variables
    Cnetastar = c
    Znetastar = -(Cnetastar)
    print 'Znetastar =', Znetastar

#Starting Iterations
    while True:
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

    #Znetastar_min = np.amin(Znetastar)
        if (Xbetastar > 0).all():
            print 'Since Znetastar has all non-negative components, the current solution is optimal'
            print Xbetastar
            break

def crissCross(A,b,c):

    print 'criss cross is called'
    inf = True
    # neta and beta specification
    N = np.asarray(pd.read_csv('csv_A.csv', header=None)).astype('float')
    print 'N = ', N
    neta = [x + 1 for x in xrange(N.shape[1])]
    # neta = np.matrix(neta)
    print 'neta = ', neta

    B = np.identity(len(b))
    print 'B =', B
    beta = [x + neta[-1] + 1 for x in xrange(B.shape[1])]
    print 'beta =', beta

    # Basic variables
    Xbetastar = b
    print 'Xbetastar =', Xbetastar

    #Non basic dual variables
    Cnetastar = c
    Znetastar = -(Cnetastar)
    print 'Znetastar =', Znetastar
    print neta
    while (inf):
        dual_inf = False
        primal_inf = False
        index = 0
        inf = False
        for i in range(len(A_T)):
            if i in beta:
                print "Checking for Basic Infeasiblity"
                for pos in range(len(beta)):
                    if beta[pos] == i+1:
                        index = pos

                if Xbetastar[index] < 0:
                    primal_inf = True
                    inf = True
                    break
            else:
                print "Checking for Dual Infeasibility"
                print neta
                for pos in range(len(neta)):
                    if neta[pos] == i+1:
                        index = pos

                if Znetastar[index] < 0:
                    print Znetastar
                    dual_inf = True
                    inf = True
                    break


        # Starting Iterations
        if dual_inf:
            #step 1
            j_pos=index
            j = neta[j_pos]
            print j

            # Step 2:
            print 'Step 2:'
            ejSize = N.shape[1]
            ej = np.zeros(ejSize)
            ej[j_pos] = 1
            ej = np.matrix(ej)
            ej = ej.T
            B_inv = inv(B)
            tempdelta_xb = np.matmul(N, ej)
            delta_xb = np.matmul(B_inv, tempdelta_xb)
            print delta_xb

            #step 3
            temp_pos = len(delta_xb)
            i_pos = 0
            for pos in range(len(delta_xb)):
                if delta_xb[pos] < 0:
                    if pos < temp_pos:
                        temp_pos = pos
                        i_pos = pos
            i = beta[i_pos]

            # step 4
            ratio = (delta_xb.item(i_pos)) / (Xbetastar.item(i_pos))
            if ratio != 0:
                t = 1 / ratio
            else:
                t = 0

            # Step 5:
            print 'Step 5:'
            eiSize = b.shape[0]
            ei = np.zeros(eiSize)
            ei[i_pos] = 1
            ei = np.matrix(ei).T
            tempdelta_zn = (np.matmul(B_inv, N))
            tempdelta_Tzn = tempdelta_zn.T
            delta_zn = np.matmul(tempdelta_Tzn, ei)
            delta_zn = -(delta_zn)
            print delta_zn

            # Step 6:
            print 'Step 7:'
            s = Znetastar.item(j_pos) / delta_zn.item(j_pos)
            print s
            print Znetastar

            # Step 7:
            print 'Step 7:'
            Xbetastar = Xbetastar - (t * delta_xb)
            Znetastar = Znetastar - (s * delta_zn)
            print Xbetastar
            print Znetastar

            # Step 8:
            print 'Step 8:'
            beta[i_pos] = j
            neta[j_pos] = i
            temp_matA = np.matrix(np.zeros(len(A)))
            temp_matA = temp_matA.T
            for i in range(len(beta)):
                if i == 0:
                    temp_matA = temp_matA + A[:, beta[i] - 1]
                else:
                    temp_matA = np.append(temp_matA, A[:, int(beta[i]) - 1], axis=1)
            B = temp_matA

            temp_matN = np.matrix(np.zeros(len(A)))
            temp_matN = temp_matN.T
            print temp_matN

            for j in range(len(neta)):
                if j == 0:
                    temp_matN = temp_matN + A[:, neta[j] - 1]
                else:
                    temp_matN = np.append(temp_matN, A[:, int(neta[j]) - 1], axis=1)
            N = temp_matN

            print 'B =', B
            print 'N =', N
            print 'beta =', beta
            print 'neta =', neta

            # new basic primal variables and nonbasic dual variables
            print Znetastar
            for p in range(len(Xbetastar)):
                if Xbetastar[p] == 0:
                    Xbetastar[p] = t
            for q in range(len(Znetastar)):
                if Znetastar[q] == 0:
                    Znetastar[q] = s
                    print Znetastar
                    print Xbetastar

                    # Znetastar_min = np.amin(Znetastar)
            if (Znetastar > 0).all():
                print 'Since Znetastar has all non-negative components, the current solution is optimal'
                break

        if primal_inf:
            # Step 1:
            print 'Step 1:'
            Xb_EnteringIndex=index
            i = beta[Xb_EnteringIndex]
            print i

            # Step 3:
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
            # Solving for s
            temp_pos = len(Znetastar)
            j_pos = 0
            for pos in range(len(delta_zn)):
                if delta_zn[pos] > 0:
                    if pos < temp_pos:
                        temp_pos = pos
                        j_pos = pos
            print "The position of exiting variable is " + str(j_pos)
            j = beta[j_pos]

            ratio = (delta_zn.item(j_pos)) / (Znetastar.item(j_pos))
            if ratio == 0:
                print "The dual is unbounded"
                break
            else:
                s = 1 / ratio
            print "The incoming variable is index " + str(i_value) + " and the outgoing variable is " + str(j)

            # Step 6:
            print 'Step 6:'
            ejSize = N.shape[1]
            ej = np.zeros(ejSize)
            ej[j_pos] = 1
            ej = np.matrix(ej)
            ej = ej.T
            B_inv = inv(B)
            tempdelta_xb = np.matmul(N, ej)
            delta_xb = np.matmul(B_inv, tempdelta_xb)
            print delta_xb

            # Step 7:
            print 'Step 7:'
            t = Xbetastar.item(Xb_EnteringIndex) / delta_xb.item(Xb_EnteringIndex)
            print t

            # Step 8:
            print 'Step 8:'
            Xbetastar = Xbetastar - (t * delta_xb)
            Znetastar = Znetastar - (s * delta_zn)
            print Xbetastar
            print Znetastar

            # Step 9:
            print 'Step 9:'
            neta[j_pos] = i
            beta[Xb_EnteringIndex] = j
            temp_matA = np.matrix(np.zeros(len(A)))
            temp_matA = temp_matA.T
            for i in range(len(beta)):
                if i == 0:
                    temp_matA = temp_matA + A[:, beta[i] - 1]
                else:
                    temp_matA = np.append(temp_matA, A[:, int(beta[i]) - 1], axis=1)
            B = temp_matA

            temp_matN = np.matrix(np.zeros(len(A)))
            temp_matN = temp_matN.T
            print temp_matN

            for j in range(len(neta)):
                if j == 0:
                    temp_matN = temp_matN + A[:, neta[j] - 1]
                else:
                    temp_matN = np.append(temp_matN, A[:, int(neta[j]) - 1], axis=1)
            N = temp_matN

            print B
            print N
            print "beta" + str(beta)
            print neta

            # new basic primal variables and nonbasic dual variables
            print Znetastar
            for p in range(len(Xbetastar)):
                if Xbetastar[p] == 0:
                    Xbetastar[p] = t
            for q in range(len(Znetastar)):
                if Znetastar[q] == 0:
                    Znetastar[q] = s
                    print Znetastar
                    print Xbetastar

                    # Znetastar_min = np.amin(Znetastar)
            if (Xbetastar > 0).all():
                print 'Since Znetastar has all non-negative components, the current solution is optimal'
                print Xbetastar
                break


A = np.asarray(pd.read_csv('csv_A.csv',header=None)).astype('float')
b = np.asarray(pd.read_csv('csv_b.csv',header=None)).astype('float')
c = np.asarray(pd.read_csv('csv_c.csv',header=None)).astype('float')
latex=start_latex()

A = np.hstack((A,np.identity(len(b))))

A = np.matrix(A)
b = np.matrix(b)
c = np.matrix(c)
latex=print_objective_function(A,b,c,latex)
A_T = A.T
save_latex("file.tex",latex)
solver(A, b, c)

