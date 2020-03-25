from scipy import stats
import math
import numpy as np
import pandas as pd
import os
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")

# input directories
data_path = ''
anno_path = ''
# output directories
output_kernels = ''
output_distances = ''

try:
    os.mkdir(output_kernels)
except:
    pass

try:
    os.mkdir(output_distances)
except:
    pass


# this is the 
def calculateDistance(point,ppoint,immuneCell):
    
    if point[0]-ppoint[0]==0.0 and point[1]-ppoint[1]==0.0:
        distance = ((immuneCell[0]-point[0])**2+(immuneCell[1]-point[1])**2)**0.5
    else:
        if point[0]-ppoint[0]==0.0:
            normaly=0.0
            normalx=1.0
            grad=0.0
        elif point[1]-ppoint[1]==0.0:
            normaly=1.0
            normalx=0.0
            grad=0.0
        else:
            grad = (point[1]-ppoint[1])/(point[0]-ppoint[0])
            normalx = 1.0/(1.0**2+1/grad**2)**0.5
            normaly = -(1.0/grad)/(1.0**2+1.0/grad**2)**0.5

        dbpoints = ((ppoint[0]-point[0])**2+(ppoint[1]-point[1])**2)**0.5
        dimmuneToPoint = ((immuneCell[0]-point[0])**2+(immuneCell[1]-point[1])**2)**0.5
        distance = abs(normalx*(point[0]-immuneCell[0]) + normaly*(point[1]-immuneCell[1]))
        if dbpoints**2+distance**2 < dimmuneToPoint**2 :
            distance = dimmuneToPoint
        if distance>dimmuneToPoint:
            difference = distance-dimmuneToPoint
            if difference>1e-7:
                print("Differnce:      "+ str(difference))
                print("point:          "+str(point))
                print("previous point: "+str(ppoint))
                print("immunceCell:    "+str(immuneCell))
                print("dbpoints:       "+str(dbpoints))
                print("dimmuneToPoint: "+str(dimmuneToPoint))
                print("distance:       "+str(distance))
                print("normalx:        "+str(normalx))
                print("normaly:        "+str(normaly))
                print("grad:           "+str(grad))
        
    return distance


# this is the more naive but approximate way to calculate immune cell to tumor distances
def calculateSDistance(point,immuneCell):
    
    distance = ((immuneCell[0]-point[0])**2+(immuneCell[1]-point[1])**2)**0.5
    return distance


data = []
    
for f in os.listdir(data_path):
    f =data_path+f
    if f.endswith('.data'):
        anno_name = f.split('/')[-1].split('.')[0] + '.data'
        if anno_name in os.listdir(anno_path):
            anno_name = anno_path+anno_name
            data.append([f,anno_name])
            
icount=0
reso = 1000000.0

for d in data:
    icount=icount+1
    print('Run Id: '+str(icount))
    print('Working on '+d[0] )
    try:
        immuneCells = []
        idata=d[0]
        adata=d[1]
        out_name = idata.split('/')[-1].split('.')[0]
        
        df = pd.read_table(idata, sep = ',', comment='#', header=None)
        df = df.dropna()
        m1 = pd.to_numeric(df[0].values,errors='coerce')
        m2 = pd.to_numeric(df[1].values,errors='coerce')

        for ind,n1 in enumerate(m1):
                immuneCells.append((n1,m2[ind]))

        xmin = m1.min()
        xmax = m1.max()
        ymin = m2.min()
        ymax = m2.max()

        if ymin > ymax:
            print('strange things are going on')
            ymin,ymax=ymax,ymin

        ttannos = []
        annos=[]
        anno=[]
        # read in borders
        anno_file = pd.read_table(adata, header=None,sep=',')
        for index,row in anno_file.iterrows():
            if '#' in str(row) and 'TUMOR' in str(row).upper():
                annos.append(anno)
                anno=[]
            elif '#' in str(row):
                anno=[]
            else:
                x = int(float(row[0])*reso)
                y = int(float(row[1])*reso)
                if math.isnan(x) or math.isnan(y):
                    continue
                anno.append((x,y))
            
        print('Done reading csv file')

        #resolution of Kernel
        res = 100j 
        X, Y = np.mgrid[xmin:xmax:res, ymin:ymax:res]
        positions = np.vstack([X.ravel(), Y.ravel()])
        values = np.vstack([m1, m2])
        #initialise Kernel
        kernel = stats.gaussian_kde(values)
        results = kernel(positions)
        Z = np.reshape(results.T, X.shape)

        fig = plt.figure(figsize=(18,12), dpi=300)
        ax = fig.add_subplot(111)
        im = np.rot90(Z)
        plot=ax.imshow(im,
                cmap=plt.cm.gist_earth_r,
                extent=[xmin, xmax, ymin, ymax])
        fig.colorbar(plot)
        if container:
            pass
        else:
            ax.plot(m1, m2, 'k.', markersize=2, alpha=0.1)

        for an in annos:
            # check if an is disconnected
            remove_list=[]
            idx_list=[]
            included_points=[]
            size = len(annos)
            if an[-1] == an[0]:
                subannos = [an]
            else:
                # remove duplicated points that appear in a row
                for ida,point in enumerate(an):
                    if ida==0:
                        continue
                    ppoint = an[ida-1]
                    if ppoint == point:
                        remove_list.append(ida)
                an = [ele for i,ele in enumerate(an) if i not in remove_list]
                
                for point in an:
                    if an.count(point)==2 and point not in included_points:

                        ids = [i for i, e in enumerate(an) if e == point]
                        included_points.append(point)
                        idx_list.append(ids)

                if idx_list:
                    flattened_idx = []
                    for ent in idx_list:
                        flattened_idx.append(ent)
                    exclude=[]
                    for ii in idx_list:
                        for i in range(ii[0],ii[1]+1,1):
                            exclude.append(i)

                    not_exclude=[]
                    not_idx_list = []
                    for ii in exclude:
                        multi = exclude.count(ii)
                        if multi > 1:
                            for jj in idx_list:
                                for ij in range(jj[0],jj[1]+1,1):
                                    if ij == ii:
                                        if ii not in not_exclude:
                                            not_exclude.append(ii)
                                        if jj not in not_idx_list:
                                            not_idx_list.append(jj)


                    for ii in not_exclude:
                        exclude.remove(ii)
                    for ii in not_idx_list:
                        idx_list.remove(ii)
                    
                    subannos = [an[i[0]:i[1]+1] for i in idx_list]

                    if 0 not in flattened_idx or (len(an)-1) not in flattened_idx:
                        subann = [an[i] for i in range(len(an)) if i not in exclude]
                        subannos.append(subann)
                        
                else:
                    subannos = [an]
            # # remove empty subannos
            subannos = [subanno for subanno in subannos if subanno != []]

            for subanno in subannos:
                subanno = [ (p[0] / reso, p[1] / reso) for p in subanno]
                # append remaining tumor annos
                ttannos.append(subanno)
                
                poly = Polygon(subanno, fill=None, edgecolor = 'r', alpha=1.0)
                ax.add_patch(poly)
            
        ax.set_title(out_name)
        ax.set_xlim([xmin, xmax])
        ax.set_ylim([ymin, ymax])
        plt.savefig(outputs_kernel+out_name+'.pdf', format='pdf')
        plt.close(fig)
        plt.clf()
        
        print("Done plotting")
            
        # get min distance distribution
        print()
        print("Generate Distances")
        fout = open(output_distances+out_name+".data","w")
        ncells = len(immuneCells)
        plain_points=[]
        for tanno in ttannos:
            for point in tanno:
                plain_points.append(point)
        nannos = len(plain_points)
        print("Found " + str(ncells) + " immune cells.")
        print("Found " + str(nannos) + " tumor anno points.")

        useSimpleD = False
        # if ncells*nannos>(150000*150000) and nannos>150000:
        #     useSimpleD = True
        #     print("Using simple distances.")
        fout.write('# Used simple distance: '+str(useSimpleD))
        
        for idxI,immuneCell in enumerate(immuneCells):
            distance_min = np.inf
            for tanno in ttannos:
                for iidx,point in enumerate(tanno):
                    # get previous point (ppoint)
                    if iidx==0:
                        ppoint = tanno[-1]
                    else:
                        ppoint = tanno[iidx-1]
                    if useSimpleD:
                        distance = calculateSDistance(point,immuneCell)
                    else:
                        distance = calculateDistance(point,ppoint,immuneCell)

                    if distance<distance_min:
                        distance_min=distance

            fout.write(str(distance_min)+'\n')
            if idxI%1000==0 :
                frac = (idxI+1.0)/(ncells*1.0)
                print("Done with fraction "+str(frac)+" of cells")
        fout.close()
    except:
        print('Failed on '+d[0]+'.')
        pass
