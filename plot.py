from scipy import stats
import cython_distance
import math
import numpy as np
import pandas as pd
import os
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")


################
# start inputs
################

# outputdirectory for plots
out_dir_plots = './outputs_true/' 
try:
   os.mkdir(out_dir_plots)
except:
   pass

# output directory for distance data
out_dir ='/media/jan/LARA-02/PGP-pancreas-JAN/Projekt/distanceData/' 
try:
   os.mkdir(out_dir)
except:
   pass

#path to folder where immune cell data is stored
data_path ='/media/jan/LARA-02/PGP-pancreas-JAN/Projekt/immuneCellData/'
imm_ext='.data'
#path to folder where anno data is stored
anno_path ='/media/jan/LARA-02/PGP-pancreas-JAN/Projekt/AnnoData/'
anno_ext='.csv'

# selection by number
rang=[]
container = []
# selection by name
part=['PDAC-PGP20']
#only plots no distances
plotOnly=False
# if plot at all
plotit=True
useSimpleD = False

#######################
# end of inputs
#######################

if plotOnly:
   plotit=True
if plotOnly:
   print("Only plotting")

data = []

for f in os.listdir(data_path):
   f =data_path+f
   if f.endswith(imm_ext):
      anno_name = f.split('/')[-1].split('.')[0] + anno_ext
      if anno_name in os.listdir(anno_path):
         anno_name = anno_path+anno_name
         data.append([f,anno_name])
            
icount=0
reso = 1000000.0

for d in data:
    
   if container:
      not_pass=True
      for con in container:
         if str(con) in d[0].split('/')[-1]:
            not_pass=False

         if not_pass:
            continue
   elif rang:
      icount=icount+1
      if icount < rang[0] or icount >= rang[1]:
         continue
   elif part:
      not_pass = True
      for p in part:
         if str(p) in d[0].split('/')[-1]:
            not_pass=False
      if not_pass:
         continue
    
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

      if plotit:
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
                  # print(point)
            if idx_list:
               
               exclude=[]
               for ii in idx_list:
                  for i in range(ii[0],ii[1]+1,1):
                     exclude.append(i)

               not_exclude=[]
               not_idx_list = []
               for ii in exclude:
                  multi = exclude.count(ii)
                  # focus on overlapping subannotations
                  if multi > 1:
                     for jj in idx_list:
                        for ij in range(jj[0],jj[1]+1,1):
                           if ij == ii:
                              # if ii not in not_exclude:
                              #    not_exclude.append(ii)
                              if jj not in not_idx_list:
                                 not_idx_list.append(jj)


               # for ii in not_exclude:
               #    exclude.remove(ii)
               
               for ii in not_idx_list:
                  for i in range(ii[0],ii[1]+1,1):
                     try:
                        exclude.remove(i)
                     except ValueError:
                        pass                  
                  idx_list.remove(ii)
                                   
               subannos = [an[i[0]:i[1]+1] for i in idx_list]

               # if 0 not in flattened_idx or (len(an)-1) not in flattened_idx:

               # add all the remaining points as one annotation
               subann = [an[i] for i in range(len(an)) if i not in exclude]
               subannos.append(subann)
                        
               # subannos = [an[i:j] for i,j in
               #             zip([0] + idx_list, idx_list +
               #                 ([size] if idx_list[-1] != size else []))]
            else:
               subannos = [an]
         # # remove empty subannos
         subannos = [subanno for subanno in subannos if subanno != []]

         for subanno in subannos:
            subanno = [ (p[0] / reso, p[1] / reso) for p in subanno]
            # append remaining tumor annos
            ttannos.append(subanno)
            if plotit:
               poly = Polygon(subanno, fill=None, edgecolor = 'r', alpha=1.0)
               ax.add_patch(poly)
      if plotit:
         ax.set_title(out_name)
         ax.set_xlim([xmin, xmax])
         ax.set_ylim([ymin, ymax])
         plt.savefig(out_dir_plots+out_name+'.pdf', format='pdf')
         plt.close(fig)
         plt.clf()
      print("Done plotting")
      if plotOnly:
         continue
            
      # get min distance distribution
      print()
      print("Generate Distances")
      fout = open(out_dir+out_name+".data","w")
      ncells = len(immuneCells)
      plain_points=[]
      for tanno in ttannos:
         for point in tanno:
            plain_points.append(point)
      nannos = len(plain_points)
      print("Found " + str(ncells) + " immune cells.")
      print("Found " + str(nannos) + " tumor anno points.")

      outstr = []
      # if ncells*nannos>(150000*150000) and nannos>150000:
      #     useSimpleD = True
      #     print("Using simple distances.")
      fout.write('# Used simple distance: '+str(useSimpleD)+'\n')
        
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
                  distance = cython_distance.calculateSDistance(point,immuneCell)
               else:
                  distance = cython_distance.calculateDistance(point,ppoint,immuneCell)

               if distance<distance_min:
                  distance_min=distance

         os = "%.2f" % distance_min+'\n'
         outstr.append(os)
         if idxI%1000==0 and idxI>0:
            frac = (idxI+1)/(ncells*1.0)
            print("Done with fraction "+str(frac)+" of cells")
            fout.write(''.join(outstr))
            outstr=[]
                
      fout.write(''.join(outstr))
      outstr=[]        
      fout.close()
   except ValueError:
      print('Failed on '+d[0]+'.')
      pass
