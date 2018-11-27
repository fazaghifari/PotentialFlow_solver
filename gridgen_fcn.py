import numpy as np
import matplotlib.pyplot as plt
from dataimport import importer
from copy import deepcopy

def pointsgen(file,domsize,npoints):

    _,_,_,airfoil = importer(file)
    halflen = int((len(airfoil)+1)/2)
    airfoil_up = airfoil[0:halflen,:]
    airfoil_up = airfoil_up[airfoil_up[:,0].argsort()[::1]]
    airfoil_low = airfoil[halflen-1:len(airfoil)+1,:]

    #Split upper and lower airfoil
    Xairfoil_up = (airfoil_up[:,0]);Yairfoil_up = (airfoil_up[:,1])
    Xairfoil_low = (airfoil_low[:,0]);Yairfoil_low = (airfoil_low[:,1])

    #Define domain
    chord = np.max(airfoil[:,0])
    dom_length = np.array([-domsize*chord, domsize*chord+chord])
    dom_width = np.array([-domsize*chord, domsize*chord])

    #Define number of nodes at each block
    xnode_out = npoints
    xnode_foil = npoints
    ynode = npoints

    #Array initialization
    #Block 1 (upleft block)
    mid1 = np.zeros(shape=[xnode_out,2]); far1 = deepcopy(mid1)
    side1 = np.zeros(shape=[ynode,2])
    #Block 2 (upright block)
    mid2 = deepcopy(mid1); far2 = deepcopy(far1);side2 = deepcopy(side1)
    #Block 3 (lowleft block)
    mid3 = deepcopy(mid1); far3 = deepcopy(far1);side3 = deepcopy(side1)
    #Block 4 (lowright block)
    mid4 = deepcopy(mid1); far4 = deepcopy(far1);side4 = deepcopy(side1)
    #Block 5 (upper airfoil)
    upfoil = np.zeros(shape=[xnode_foil,2]); upfar = deepcopy(upfoil)
    upleft = np.zeros(shape=[ynode,2]); upright = deepcopy(upleft)
    #Block 6 (lower airfoil)
    lowfoil = deepcopy(upfoil); lowfar = deepcopy(upfar)
    lowleft = deepcopy(upleft); lowright = deepcopy(upright)

    #Fill the inner blocks
    #Define coefficients
    z1 = -1;z2 = 1
    e1 = -1;e2 = 1
    d5 = 2 ;d6 = 2

    #Block 5
    block5 = [[[]for ii in range(xnode_foil)]for ii in range(ynode)]

    for ii in range(0,ynode): #side points
        eta5 =(np.exp(e2*(ii)/(ynode-1))-1)/(np.exp(e2)-1)
        upleft[ii,:] = np.array([0,(0*(1-eta5)+eta5*dom_width[1])])
        upright[ii, :] = np.array([1, (0 * (1 - eta5) + eta5 * dom_width[1])])
        block5[ii][0] = upleft[ii,:];block5[ii][xnode_foil-1] = upright[ii,:]

    for ii in range(0,xnode_foil): #far points
        zeta5 = (1-np.cos(np.pi*(ii)/(xnode_foil-1)))/d5
        upfar[ii,:] = np.array([(0*(1-zeta5)+zeta5*1),dom_width[1]])
        block5[0][ii] = upfar[ii,:]

    for ii in range(0,xnode_foil): #airfoil points
        bool5 = 0; k = 0
        upfoil[ii,0] = deepcopy(upfar[ii,0])
        while bool5 == 0:
            p = upfoil[ii,0]-Xairfoil_up[k]
            q = upfoil[ii,0]-Xairfoil_up[k+1]
            if p>0 and q<0:
                bool5 = 1
                x_in1 = Xairfoil_up [k]; y_in1 = Yairfoil_up[k]
                x_in2 = Xairfoil_up[k+1];y_in2 = Yairfoil_up[k+1]
                upfoil[ii,1] = y_in1 +((y_in2-y_in1)*(upfoil[ii,0]-x_in1)/(x_in2-x_in1))
            elif q==0:
                bool5 = 1
                upfoil[ii,1] = Yairfoil_up[k+1]
            elif p==0:
                bool5 = 1
                upfoil[ii,1] = Yairfoil_up[k]
            k = k+1
        block5[ynode-1][ii] = upfoil[ii,:]
    #Inner block using TFI
    for ii in range(1,xnode_foil-1):
        for jj in range(1,ynode-1):
            eta5 = (np.exp(e2*(jj)/(ynode-1))-1)/(np.exp(e2)-1)
            k = (1-eta5)*upfoil[ii,1] + eta5*upfar[ii,1]
            block5[jj][ii] = np.array([upfoil[ii,0],k])

    #Block 6
    block6 = [[[]for ii in range(xnode_foil)]for ii in range(ynode)]

    for ii in range(0,ynode): #side points
        eta6 =(np.exp(e1*(ii)/(ynode-1))-1)/(np.exp(e1)-1)
        lowleft[ii,:] = np.array([0,(dom_width[0]*(1-eta6)+eta6*0)])
        lowright[ii, :] = np.array([1, (dom_width[0]* (1 - eta6) + eta6 *0 )])
        block6[ii][0] = lowleft[ii,:];block6[ii][xnode_foil-1] = lowright[ii,:]

    for ii in range(0,xnode_foil): #far points
        zeta6 = (1-np.cos(np.pi*(ii)/(xnode_foil-1)))/d6
        lowfar[ii,:] = np.array([(0*(1-zeta6)+zeta6*1),dom_width[0]])
        block6[ynode-1][ii] = lowfar[ii,:]

    for ii in range(0,xnode_foil): #airfoil points
        bool6 = 0; k = 0
        lowfoil[ii,0] = deepcopy(lowfar[ii,0])
        while bool6 == 0:
            p = lowfoil[ii,0]-Xairfoil_low[k]
            q = lowfoil[ii,0]-Xairfoil_low[k+1]
            if p>0 and q<0:
                bool6 = 1
                x_in1 = Xairfoil_low[k]; y_in1 = Yairfoil_low[k]
                x_in2 = Xairfoil_low[k+1];y_in2 = Yairfoil_low[k+1]
                lowfoil[ii,1] = y_in1 +((y_in2-y_in1)*(lowfoil[ii,0]-x_in1)/(x_in2-x_in1))
            elif q==0:
                bool6 = 1
                lowfoil[ii,1] = Yairfoil_low[k+1]
            elif p==0:
                bool6 = 1
                lowfoil[ii,1] = Yairfoil_low[k]
            k = k+1
        block6[0][ii] = lowfoil[ii,:]
    #Inner block using TFI
    for ii in range(1,xnode_foil-1):
        for jj in range(1,ynode-1):
            eta6 = (np.exp(e1*(jj)/(ynode-1))-1)/(np.exp(e1)-1)
            k = (1-eta6)*lowfar[ii,1] + eta6*lowfoil[ii,1]
            block6[jj][ii] = np.array([lowfoil[ii,0],k])

    #Fill the outer box
    #Block1
    block1 = [[[]for ii in range(xnode_foil)]for ii in range(ynode)]
    for ii in range(0,xnode_out):
        zeta1 = (np.exp(z1*(ii)/(xnode_out-1))-1)/(np.exp(z1)-1)
        far1[ii,:] = np.array([(dom_length[0]*(1-zeta1)+zeta1*0), dom_width[1]])
        mid1[ii, :] = np.array([(dom_length[0] * (1 - zeta1) + zeta1 * 0), 0])
        block1[0][ii] = far1[ii,:]; block1[ynode-1][ii] = mid1[ii,:]

    for ii in range(0,ynode):
        eta1 = (np.exp(e1*(ii)/(ynode-1))-1)/(np.exp(e1)-1)
        side1[ii,:] = np.array([dom_length[0],(dom_width[1]*(1-eta1)+eta1*0)])
        block1[ii][0] = side1[ii,:];block1[ii][xnode_out-1] = upleft[ii,:]

    #Inner Block
    for jj in range(1,ynode-1):
        for ii in range (1,xnode_out-1):
            block1[jj][ii] = np.array([far1[ii,0],side1[jj,1]])

    # Block2
    block2 = [[[] for ii in range(xnode_foil)] for ii in range(ynode)]
    for ii in range(0, xnode_out):
        zeta2 = (np.exp(z2 * (ii) / (xnode_out - 1)) - 1) / (np.exp(z2) - 1)
        far2[ii, :] = np.array([( 1* (1 - zeta2) + zeta2 * dom_length[1]), dom_width[1]])
        mid2[ii, :] = np.array([( 1* (1 - zeta2) + zeta2 * dom_length[1]), 0])
        block2[0][ii] = far2[ii, :];block2[ynode - 1][ii] = mid2[ii, :]

    for ii in range(0, ynode):
        eta2 = (np.exp(e1 * (ii) / (ynode - 1)) - 1) / (np.exp(e1) - 1)
        side2[ii, :] = np.array([dom_length[1], (dom_width[1] * (1 - eta2) + eta2 * 0)])
        block2[ii][0] = upright[ii, :]; block2[ii][xnode_out - 1] = side2[ii, :]

    # Inner Block
    for jj in range(1, ynode - 1):
        for ii in range(1, xnode_out - 1):
            block2[jj][ii] = np.array([far2[ii, 0], side2[jj, 1]])

    #Block3
    block3 = [[[]for ii in range(xnode_foil)]for ii in range(ynode)]
    for ii in range(0,xnode_out):
        zeta3 = (np.exp(z1*(ii)/(xnode_out-1))-1)/(np.exp(z1)-1)
        far3[ii,:] = np.array([(dom_length[0]*(1-zeta3)+zeta3*0), dom_width[0]])
        mid3[ii, :] = np.array([(dom_length[0] * (1 - zeta3) + zeta3 * 0), 0])
        block3[0][ii] = mid3[ii,:]; block3[ynode-1][ii] = far3[ii,:]

    for ii in range(0,ynode):
        eta3 = (np.exp(e2*(ii)/(ynode-1))-1)/(np.exp(e2)-1)
        side3[ii,:] = np.array([dom_length[0],(0*(1-eta3)+eta3*dom_width[0])])
        block3[ii][0] = side3[ii,:];block3[ii][xnode_out-1] = lowleft[ii,:]

    #Inner Block
    for jj in range(1,ynode-1):
        for ii in range (1,xnode_out-1):
            block3[jj][ii] = np.array([far3[ii,0],side3[jj,1]])

    #Block4
    block4 = [[[]for ii in range(xnode_foil)]for ii in range(ynode)]
    for ii in range(0,xnode_out):
        zeta4 = (np.exp(z2*(ii)/(xnode_out-1))-1)/(np.exp(z2)-1)
        far4[ii,:] = np.array([(1*(1-zeta4)+zeta4*dom_length[1]), dom_width[0]])
        mid4[ii, :] = np.array([(1* (1 - zeta4) + zeta4 * dom_length[1]), 0])
        block4[0][ii] = mid4[ii,:]; block4[ynode-1][ii] = far4[ii,:]

    for ii in range(0,ynode):
        eta4 = (np.exp(e2*(ii)/(ynode-1))-1)/(np.exp(e2)-1)
        side4[ii,:] = np.array([dom_length[1],(0*(1-eta4)+eta4*dom_width[0])])
        block4[ii][0] = lowright[ii,:];block4[ii][xnode_out-1] = side4[ii,:]

    #Inner Block
    for jj in range(1,ynode-1):
        for ii in range (1,xnode_out-1):
            block4[jj][ii] = np.array([far4[ii,0],side4[jj,1]])


    #Block merging
    blockup = block1[:][:]+block5+block2[:][:]
    blocklow = block3[:][:]+block6+block4[:][:]
    block = blockup+blocklow

    flat_list = []
    for ii in block:
        for jj in ii:
            flat_list.append(jj)

    coord = np.zeros(shape=[len(flat_list),2])
    for ii in range(0,len(flat_list)):
        coord[ii,0] = flat_list[ii][0]
        coord[ii, 1] = flat_list[ii][1]
    coord = coord[coord[:, 1].argsort()[::1]]
    # print(flat_list[1799])
    plt.scatter(coord[:,0],coord[:,1],s=1)
    # plt.grid(True)
    plt.show()
    return coord
