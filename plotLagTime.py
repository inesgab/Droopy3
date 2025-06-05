# In[]:
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import os, glob, sys

module_path = os.path.abspath(os.path.join('./'))
if module_path not in sys.path:
    sys.path.append(module_path)
from functionsCleanPipeline import *

# In[]:
#get the halftime of the logistic growth to measure the lag
plt.close('all')
slope=[]
inc=[]

rootPathList = ['/Users/inesgabert/Documents/LBE/experiences/RFP_inoculum/']
source = '/Users/inesgabert/Documents/LBE/experiences/'
        
for rootPath in rootPathList:
    folder=sorted([join(join(rootPath,o),'analysis/') for o in listdir(rootPath) if (o[:3]=='201' or o[:3]=='EXP') and os.path.isdir(os.path.join(rootPath, o))])
    channelList =['RFP']
    #incCarte=parameters.incCarte.values #incertitude de la mesure par la carte hardware
    #timeThresh=parameters.timeThresh.values #heures max d'analyse
    #threshOutlayer=parameters.threshOutlayer.values  #heures pour mesure de contamination. Pour desactiver mettre une val negative
    
    #startTime=parameters.startTime.values #time to start analysing the curves to get ride of the weird shit at the beginning
    #nbReps=parameters.nbReps.values
    ctclr =-1
    calc = True    


            
    for channel in channelList:
        
        fig,ax = plt.subplots()
        fig2,ax2 = plt.subplots()
        fig3,ax3 = plt.subplots()
        
        for path in folder:
            [dropMap, df, label]=loadData(path)
            parameters= pd.read_csv(path+'parametersAnalysis.csv')
            parameters = parameters.set_index('channelList')
            parameters['display'] = False
            
            if calc==True:
                print('calculate the interpolation of data in a single time vector for all drop')
                poolDataInterpolate(dropMap, df, label, channel,path, parameters, savefig=False)#used after to measure halftime
                checkDataInterpolated(path, label,channel)   
            else:
                print('no calculation of interpolation') 
                
            print('plotting the results')
            thresh = getThresHalfTime(rootPath,source,label)
          
            for l in label:

                file = path+l+channel
                print(file)
                df = pd.read_csv(file+'Interp.csv')
                
                t = df['time'].to_list()
                df = df.drop('time', axis=1)
                ih = df[df>thresh].apply(pd.Series.first_valid_index )
                ih = ih.dropna()
                tl = [t[int(i)] for i in ih.values]
                dft = pddf({l:tl},index=ih.index) 
                dft.to_csv(path+'resultIndiv/'+l+'halftime'+channel+'.csv')
  

#In[]:

#plot halftime as a function of N0 to extract lag time t_lag = t_theta - (ln theta  – ln N_0)/Lambda
channelList =['RFP']
Vdrop = 4e-4 #ml


def calcLag(h,g,C0,threshN):
    lag = h -  np.log(threshN/C0)/g
    return lag


for rootPath in rootPathList:
    print('rootPath:', rootPath)
    print('listdir(rootPath):', listdir(rootPath))
    folder=sorted([join(join(rootPath,o),'analysis/') for o in listdir(rootPath) if (o[:3]=='201' or o[:3]=='EXP') and os.path.isdir(os.path.join(rootPath, o))])
    print('folder:', folder)
    fig,ax = plt.subplots()
    clr = matplotlib.cm.Dark2(np.linspace(0, 1, len(folder)))
    
    for channel in channelList:
        jitter=-0.1
        idxColor=0
        for path in folder:
            tlag = pddf()
            N0 = list()
            [dropMap, df, label]=loadData(path)
            [stdgRate, gRate, lag, stdlag, yld, stdyield, halfTime] = getDataHistogram(label,path,channel,True)
            label = sorted(label)
            for l in label:
                h = halfTime[[l+'_ht',l]].dropna()
                print(path,l)
                inoc = getValueInLabel(l,path)
                if inoc<1:
                    inoc = 1
                    
                C0 = inoc/Vdrop
                #print('{:.5E}'.format(C0))
                N0.append(getValueInLabel(l,path))
                threshN = calib(getThresHalfTime(rootPath,source,label),list(map(float,ast.literal_eval(parameters.calib[channel]))))
                t=pddf()
                t[l] =  calcLag(h[l+'_ht'],h[l],C0,threshN)
                print('tlag:', t[l], 'ht:', h[l+'_ht'],  'h:', h[l],  'C0:', C0, 'threshN:', threshN)
                
                tlag = pd.concat([tlag,t[l]], axis=1)
                
            pts = list()
            for (columnName, columnData) in tlag.items():
                p = columnData.values
                pts.append(p[~(pd.isnull(p))])
                        
            flierprops = dict(marker='o', markerfacecolor=clr[idxColor], markersize=2, markeredgecolor=clr[idxColor], alpha=0.5)
                        
            box1 = plt.boxplot(pts, positions=np.log(N0)+jitter,showfliers=True, widths=.2,whis=0.2, flierprops=flierprops)
            for item in ['boxes', 'whiskers', 'fliers', 'medians', 'caps']:
                plt.setp(box1[item], color=clr[idxColor])
                plt.setp(box1[item], alpha=0.8)
                plt.setp(box1[item], lw=2)
            jitter+=0.1
            idxColor+=1
     
            pathPlot = path+'plotIndiv/'
            if not os.path.exists(pathPlot):
                os.makedirs(pathPlot)
            tlag.to_csv(pathPlot+'result_tlagGood_'+ channel+'.csv')
    
    axis_font = {'fontname':'Courrier', 'size':'14'}
    plt.ylabel('lag (h)',**axis_font)
    yMin = 0 
    yMax = 9
    plt.ylim([yMin,yMax])
    grid_y_ticks_minor = np.arange(yMin, yMax, 1 )
    grid_y_ticks_major = np.arange(yMin, yMax, 2 )
    plt.xticks(np.log(N0))
    ax.set_xticklabels(list(N0), rotation=45, ha='right')
    plt.xlabel('inoculum on logscale') 
    ax.set_yticks(grid_y_ticks_major)
    ax.set_yticks(grid_y_ticks_minor, minor=True)
    ax.grid(which='both')
    date=rootPath.split('/')
    plt.title(date[-2]+' lag '+channel)
    fig.savefig(rootPath+'t_lagVSN0_'+channel+'.pdf', format='pdf', dpi=500, bbox_inches='tight')


# In[]:
#plot le boxplot du lag, gRate, yield en allant chercher les bons label directement et en les ordonnant
#sauve le plot qui montre tout les replicat dans rootpath
#remove the outlayers


for rootPath in rootPathList:
    folder=sorted([join(join(rootPath,o),'analysis/') for o in listdir(rootPath) if (o[:3]=='201' or o[:3]=='EXP') and os.path.isdir(os.path.join(rootPath, o)) ])


    #listParam=['stdgRate','gRate', 'lag', 'yld']
    #listParam=['gRate', 'lag', 'yld']
    listParam=['lag']
    channel='RFP'  #GFP RFP PVD
    logX = True
    xlabel = 'inoculum '+channel+ ' (cell/drop)'

    df_data_alldate = pddf()

    for kindOfData in listParam:

        for path in folder:
            print(path)
            
            df_data = pddf()

            a=path.split('/')
            date=a[-3]
            data=[]
            cleanedData=[]
            listData=[]
            listLabel=[]

            [nn, nnn, label]=loadData(path)
            [stdgRate, gRate, lag, stdlag, yld, stdyield, halfTime]=getDataHistogram(label,path,channel,fromLag=True)
            for l in label:

                #defini le yrange
                if kindOfData == 'gRate':
                    data=gRate
                    if '/CAA/' in path or '/CAAnoGly/' in path  or '/CAAPvdS/' in path or '/CAA1/' in path:
                        yMin=0.5
                        yMax=2
                    if 'M9gly' in path:
                        yMin=0
                        yMax=0.5
                    #yMax=1 #GFP
                    if 'BipyCAA' in path:
                        yMin=0.2
                        yMax=1
                    grid_y_ticks_minor = np.arange(yMin, yMax, 0.01 )
                    grid_y_ticks_major = np.arange(yMin, yMax, 0.1 )
                        
                elif kindOfData =='lag':
                    data=lag
                    yMin=0
                    yMax=23
                    ylabel = 'lag time (h)'
                    if 'BipyCAA' in path:
                        yMin=0
                        yMax=35
                    grid_y_ticks_minor = np.arange(yMin, yMax, 0.5 )
                    grid_y_ticks_major = np.arange(yMin, yMax, 1 )
                    
                elif kindOfData == 'yld':
                    data = yld
                    ylabel = 'maximal concentration (cell/ml)'
                    if '/CAA/' in path or '/CAAnoGly/' in path or '/CAAPvdS/' in path or '/CAA1/' in path or '/Maelle/' in path:
                        scale = 1e9
                        yMin= -.1*scale#5.7
                        yMax= 6*scale#6.8
                        data=calib(yld,list(map(float,ast.literal_eval(parameters.calib[channel])))) #convert fluo_area to cell/ml
                        grid_y_ticks_minor = np.arange(yMin, yMax, scale/2 )
                        grid_y_ticks_major = np.arange(yMin, yMax, scale )
                    if 'M9gly' in path:
                        yMin=5
                        yMax=7   
                        data=yld 
                    if 'BipyCAA' in path:
                        yMin=2
                        yMax=7
                        grid_y_ticks_minor = np.arange(yMin, yMax, 0.1 )
                        grid_y_ticks_major = np.arange(yMin, yMax, 0.2 )
                        
                elif kindOfData == 'stdgRate':
                    data=stdgRate
                    if '/CAA/' in path or '/CAAnoGly/' in path  or '/CAAPvdS/' in path or '/CAA1/' in path:
                        yMin=0
                        yMax=0.06
                    if 'M9gly' in path:
                        yMin=0
                        yMax=0.02
                    
                    if 'BipyCAA' in path:
                        yMin=0
                        yMax=0.06
                        
                    grid_y_ticks_minor = np.arange(yMin, yMax, 0.001 )
                    grid_y_ticks_major = np.arange(yMin, yMax, 0.01 )
                else :
                    print('do not get kind of data')

                #get the non nan data for every label
                cleanedData = [x for x in data[l] if str(x) != 'nan']
                listData.append(cleanedData)
                #listLabel.append(l+' N='+str(len(cleanedData))) not used get the label from df instead
                df_data[ getValueInLabel(l,path) ] =  pd.Series(cleanedData)
                

            df_data['date'] = pd.Series( df_data.shape[0] * [path.split('/')[-3]]) 
            df_data = df_data.reindex(natsorted(df_data.columns), axis=1)
            
            df_data_alldate = pd.concat([df_data_alldate,df_data])
            df_plot = df_data_alldate.melt(id_vars='date', value_vars=df_data.columns)

            


        df_data_alldate.to_csv(rootPath+'_df_'+ kindOfData +'_'+channel+'.csv')
       
        
       
        fig, ax = plt.subplots()
        sns.boxplot(x='variable', y='value', hue='date', data=df_plot, showfliers = False)
        
     
        
   
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
        ax.set_xlabel(xlabel) 
        ax.set_ylim([yMin,yMax])
        ax.set_ylabel(ylabel)
        ax.set_yticks(grid_y_ticks_major)
        ax.set_yticks(grid_y_ticks_minor, minor=True)
        ax.grid(which='both')


        fig.savefig(rootPath+'_boxplot_'+ kindOfData +'_'+channel+'.pdf',bbox_inches='tight', format='pdf')
            
          
        df_plot2 = pddf()
        df_data_mean = pddf()
        df_data_std = pddf()
        
        df_data_mean = df_data_alldate.groupby('date').mean()
        df_data_mean['date'] = df_data_mean.index
        df_plot2 = df_data_mean.melt(id_vars='date', value_vars=df_data_mean.columns, value_name='mean')

        df_data_std = df_data_alldate.groupby('date').std()
        df_data_std['date'] = df_data_std.index
        df_data_std = df_data_std.melt(id_vars='date', value_vars=df_data_std.columns, value_name='std')
        df_plot2['std'] = df_data_std['std']
        
  
        
        fig,ax = plt.subplots()
        
        lvls = df_plot2.date.unique()
        for i in lvls:
            ax.errorbar(x = df_plot2[ df_plot2['date']==i]["variable"],
                        y = df_plot2[ df_plot2['date']==i]["mean"], 
                        yerr = df_plot2[ df_plot2['date']==i]["std"],
                        label=i, 
                        linestyle = '',
                        fmt='o')
            
        ax.legend()
        
        ax.set_xlabel(xlabel) 
        ax.set_ylim([yMin,yMax])
        ax.set_ylabel(ylabel)

        
        fig.savefig(rootPath+'_errorbar_'+ kindOfData +'_'+channel+'.pdf',bbox_inches='tight', format='pdf')
