#!/usr/bin/python3.6

# Copyrigth 2019 Saverio Vicario 
# This file is part of PPresil-BM.

    # PPresil-BM is free software: you can redistribute it and/or modify
    # it under the terms of the GNU General Public License as published by
    # the Free Software Foundation, either version 3 of the License, or
    # (at your option) any later version.

    # PPresil-BM is distributed in the hope that it will be useful,
    # but WITHOUT ANY WARRANTY; without even the implied warranty of
    # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    # GNU General Public License for more details.

    # You should have received a copy of the GNU General Public License
    # along with PPresil-BM.  If not, see <http://www.gnu.org/licenses/>.



from sklearn.preprocessing import FunctionTransformer, PolynomialFeatures
from linear_model import BayesianLinearModel, _model_evidence, _negative_log_marginal_likelihood
from patsy.contrasts import Sum, Treatment
import copy
import numpy as np
import pandas as pd
import xarray
import xarray as xr
from pyproj import CRS, transform
from scipy.stats import circmean, circvar


	

def LoadXarray(path, maxvalue=10000,navalue=-32768, bandNameFormatter=pd.to_datetime):
    TDG=xarray.open_rasterio(path)
    TDG.coords["band"]=bandNameFormatter([x.strip() for x in TDG.attrs['band_names'].split(",")])
    if TDG.dtype.name.find("int")>-1:
        attrs=TDG.attrs
        TDG=TDG.where(TDG!=navalue)
        TDG=TDG/maxvalue
        TDG.attrs=attrs
    if TDG.attrs["coordinate_system_string"].find("UTM")!=-1:
        TDG.coords["x"].attrs["units"]=TDG.coords["y"].attrs["units"]="meters"
    TDG=TDG.rename({"band":"Time","x":"X","y":"Y"})
    TDG.coords["X"].attrs["long_name"]="Longitude"
    TDG.coords["Y"].attrs["long_name"]="Latitude"
    return TDG

def MakeLinearSinus(t, k=[1,2,3], trend=False, YAM=False, freq=365, offset=0):
    """
    Actual Harmonic model expressed as  transformation of time variable necessary 
    to fit the data within a linear framework
    """
    N=len(t)
    K=np.repeat(np.array([k]),N, axis=0)
    fix=2*np.pi/freq
    Fix=t.reshape(N,1)*fix*K
    #print(Fix.shape)
    names=[]
    if trend:
        Xm=np.concatenate([np.array([1]*N).reshape(N,1),
                           t.reshape(N,1),np.sin(Fix), np.cos(Fix)],axis=1)
        names=["intercept","trend"]+["sin"+str(k) for k in K ]+["cos"+str(k) for k in K ]
    elif YAM:
        year=(((t+offset)/365).astype("int")).flatten()
        #year = Sum().code_without_intercept(list(set(year))).matrix[year,:]
        #print(Treatment().code_without_intercept(list(set(year))))
        #To deal with year that have no observation make an index that guide observation to correct contrast
        ex=Treatment().code_without_intercept(list(set(year)))
        #print(year[0])
        #print(ex.matrix.shape)
        cast=pd.Series(np.arange(ex.matrix.shape[0]), index=[0]+[int(x.split(".")[1][:-1]) for x in ex.column_suffixes])
        #print(cast)
        #inc=cast.loc[year]
        year = ex.matrix[cast.loc[year],:]
        #print(year.shape[1])
        #year=sm.tools.categorical(year, drop=True)
        Xm=np.concatenate([np.array([1]*N).reshape(N,1),
                           year,np.sin(Fix), np.cos(Fix)],axis=1)
        names=["intercept"]+["year"+str(k) for k in 2+np.arange(year.shape[1]) ]+["sin"+str(kk) for kk in k ]+["cos"+str(kk) for kk in k ]
    else:
        Xm=np.concatenate([np.array([1]*N).reshape(N,1),
                           np.sin(Fix), np.cos(Fix)],axis=1)
        names=["intercept"]+["sin"+str(kk) for kk in k ]+["cos"+str(kk) for kk in k ]
    return Xm,names

def RunnerBayesian(values,times, change=True,PP=False, name=False, k=[1,2,3],
                   fulloutput=False,reps=100, Event=False, BFT=1.6,freqsubsample=10,reformat_dict={}, Expected=False, Loutput=None, dense=False, offset=0, minmaxFeat=False, **kwargs):
    """
    Run a pixel columns and get output
    """
    freq=365
    Values=values.astype("float").squeeze()
    if reformat_dict:
        Values=values.astype("float").squeeze()/reformat_dict["maxvalues"]
        Values[reformat_dict["nan"]==values.squeeze()]=np.nan
    #print(Values.shape)
    fullyear=np.ceil(times[-1]/freq).astype("int")
    times=times[np.isfinite(Values)]
    VValues=Values[np.isfinite(Values)]
    #delta=(times[0]-times[-1])
    #if int(delta/freq)==int(delta/freq):
    #    times=times[:-1]
    #    values=values[:-1]
    if dense:
        Model=SingleModel
    elif change:
        Model=MultiModel
    else:
        Model=MonoModel
    try:
        models, evidences, check=Model(times,VValues,freq=freq, k=k,BFT=BFT, offset=offset)
    except:
        #to handle nan-only columns
        return np.array([np.nan]*Loutput)
    modelA=models
    if Expected:
        Years=np.ceil(times.max()/freq).astype("int")
        nobs=np.arange(0,Years*freq,freqsubsample)
        tempE, tempS= models.predict(nobs,variance=True)
        temp=np.concatenate([tempE,tempS]).flatten()
        if name:
            NameExp=np.concatenate(["Expected_"+pd.Series(nobs).astype("str"),"ExpectedSD_"+pd.Series(nobs).astype("str")])
            temp=[temp,NameExp]
        if fulloutput:
            return temp, models
        return temp
    Res, Rest=FeatureBayes3(modelA, obsdays=times, obsvalue=values, freq=freq,PP=PP, Event=Event,reps=reps,freqsubsample=freqsubsample, YearsL=fullyear, offset=offset, minmaxFeat=minmaxFeat)
    #R=Res.stack()
    #R[R==-1]=np.nan
    R=Res
    if PP:
        RT=Rest.stack()
        R=pd.concat([R,RT])
    if fulloutput:
        return models, Res, Rest,evidences, check
    if name:
        name=["_".join(x) for x in R.index]
        if PP:
            name+=["_".join(x) for x in RT.index]
        return R, name
    else:
        return R

def PartialModelEvidence2(model, T,y):
    """
    Calculate the increase in evidence adding new data, without modifying original model
    """
    Imodel=copy.deepcopy(model)
    Imodel.update(T[:,np.newaxis],y[:,np.newaxis])
    base=0
    if model._BayesianLinearModel__initialised:
        base=model.evidence()
    return Imodel.evidence()-base
def SingleModel(t,YY,freq,k=[1,2,3], BFT=1.6, poly=False, offset=0):
    """
    Use a different YSM (year seasonal model) for each year
    To be used on temporally dense data
    """
    HarmonicA = FunctionTransformer(lambda x: MakeLinearSinus(x, k=k,YAM=False)[0], validate=False)
    if poly:
        from sklearn.preprocessing import PolynomialFeatures
        #HarmonicA = FunctionTransformer(lambda x: PolynomialFeatures(x.shape[0]-3))
        Polymodel = FunctionTransformer(lambda x: PolynomialFeatures(3).fit_transform(x), validate=False)
        Polymodel=PolynomialFeatures(3)
        #PolyModel = BayesianLinearModel(basis=lambda x: polybasis(x, t))
    yearStart=0
    models=[]
    overlap=0
    for year in np.arange(yearStart,np.ceil(t.max()/freq)):
        timefilter=((t+offset)>(-overlap+freq*year))&((t+offset)<(+overlap+freq*(year+1)))
        T=t[timefilter]
        if len(T)==0:
            print(year)
            continue
        #print(T+np.array(["2009-01-01"]).astype("datetime64[D]"))
        y=YY[timefilter]
        modelA=BayesianLinearModel(basis=HarmonicA.fit_transform)
        try:
            modelA.update(T[:,np.newaxis],y[:,np.newaxis])
        except :
            modelA=BayesianLinearModel(basis=Polymodel.fit_transform)
            modelA.update(T[:,np.newaxis],y[:,np.newaxis])
        models.append(modelA)
    
    return ModelList(models, offset=offset), [],[]

def MonoModel(t,YY,freq,k=[1,2,3], BFT=1.6, YAM=True, offset=0):
    """
    Use a YAM ( Year Anomaly model ) across all years
    """
    HarmonicA = FunctionTransformer(lambda x: MakeLinearSinus(x, k=k,YAM=YAM)[0], validate=False)
    modelA=BayesianLinearModel(basis=HarmonicA.fit_transform)
    modelA.update(t[:,np.newaxis],YY[:,np.newaxis])
    return modelA, [],[]
def MultiModel(t,YY,freq=365, k=[1,2,3],BFT=1.6, offset=0):
    """
    Fit a model year by year switching between the previous year  model, a new Seasonal model and a  new flat trend model
    """
    HarmonicA = FunctionTransformer(lambda x: MakeLinearSinus(x, k=k,YAM=True)[0], validate=False)
    modelA=BayesianLinearModel(basis=HarmonicA.fit_transform)
    nameparam=MakeLinearSinus(t[:,np.newaxis], k=k,YAM=True, freq=365)[1]
    noyearparam=[x for x,y in enumerate(nameparam) if y.find("year")!=0]
    modelA.update(t[:,np.newaxis],YY[:,np.newaxis])
    Harmonic = FunctionTransformer(lambda x: MakeLinearSinus(x, k=k,YAM=False)[0], validate=False)
    g=len(t)
    startingPoint={"location":modelA.location[noyearparam],
                    "dispersion":modelA.dispersion[noyearparam,:][:,noyearparam]*(g**0.5),
                    "scale":modelA.scale,
                   "shape":modelA.shape}
    
    #     startingPoint={"location":modelA.location,
    #                 "dispersion":modelA.dispersion*g,
    #                 "scale":modelA.scale,
    #                "shape":modelA.shape}
    model=BayesianLinearModel(basis=Harmonic.fit_transform, **startingPoint)
    yearStart=0
    #therdshold for bayesian factor, assuming natural log of the ratio
    #BFT=4
    evidences=[np.array([-np.inf])]
    evits=[]
    eviss=[]
    evifs=[]
    maxvalue=YY.max()
    cont=[]
    models=[]
    import copy
    for year in np.arange(yearStart,np.ceil(t.max()/freq)):
    #for year in np.arange(3):
        overlap=0
        timefilter=((t+offset)>(-overlap+freq*year))&((t+offset)<(+overlap+freq*(year+1)))
        T=t[timefilter]
        y=YY[timefilter]
        evit=PartialModelEvidence2(model, T,y)
        evits.append(evit)
        #New flat model
        evif=np.nan
        if len(T)>3:
            modelF=BayesianLinearModel(basis=lambda x:np.concatenate([np.array([1]*len(x))[:,np.newaxis],x], axis=1))
            modelF.update(T[:,np.newaxis],y[:,np.newaxis])
            if modelF.location[0]<0.2:
                evif=modelF.evidence()
        evifs.append(evif)
        evis=np.nan
        #NewHarmonic model
        if len(T)>2:
            #modelS=BayesianLinearModel(basis=Harmonic.fit_transform)
            #startingPoint={}
            modelS=BayesianLinearModel(basis=Harmonic.fit_transform, **startingPoint)
        
            modelS.update(T[:,np.newaxis],y[:,np.newaxis])
            evis=modelS.evidence()
        eviss.append(evis)
        #Decision
        if evis-evit>BFT:
            model=modelS
            cont.append("new")
        #added this rule
        elif evit-evis >BFT:
            model=model
            model.update(T[:,np.newaxis],y[:,np.newaxis])
            cont.append("old")

        elif evif-evit>BFT:
            model=modelF
            cont.append("flat")
        #modified the uncertain outcome
        else:
            model=modelS
            #print(model.location)
            cont.append("new")

        evi=model.evidence()
        models.append(model)    
        evidences.append(np.array([evi]))
    evidences=[np.concatenate(evidences, axis=0)]
    evidences=evidences[0][1:]
    return ModelList(models, offset=offset), evidences, (eviss,evits, evifs, cont)
def FeatureBayes3(modelA, obsdays, obsvalue, freq=365,PP=True, Event=False,reps=100,freqsubsample=5, YearsL=None, offset=0, minmaxFeat=False):
    """
    Extract 4 descriptors of the time series from a model or multiyear model and perfrom a PP test of goodness of fit 
    """

    t=obsdays
    YY=obsvalue
    Years=np.ceil(t.max()/freq).astype("int")
    freqsubsample=freqsubsample
    
    daystoobs=np.arange(0,Years*freq, freqsubsample)
    X=modelA._BayesianLinearModel__basis(daystoobs)
    
    try:
        params=Randomizer(modelA,daystoobs,reps)
        #modelA._predict_mean
    except ValueError:
        print(obsvalue)
    
    #Dati attesi per la frequenza di osservazione ideale
    pred0=pd.DataFrame(np.nansum(X[:,np.newaxis,:]*params,axis=2).T)
    Res=[]
    Rest=[]
    varyear=0
    meansdyear=[]
    totdays=pred0.shape[1]
    year=0
    #Statistiche sulle date realmente osservati
    if PP or Event:
        X=modelA._BayesianLinearModel__basis(t)
        #pred0t=pd.DataFrame(np.array([np.nansum(X*Randomizer(modelA,t),axis=1) for x in np.arange(reps)]))
        paramst=Randomizer(modelA,t,reps)
        pred0t=pd.DataFrame(np.nansum(X[:,np.newaxis,:]*paramst,axis=2).T)
        #pred0t=pred0.iloc[:,t-1]
        yeart=((t+offset)/freq).astype("int")
    #giorni idealmente osservati
    yearobs=((daystoobs+offset)/freq).astype("int")
    for year in np.arange(Years):
        def Loadingres(name, index, value, year=year,Res=Res):
            Res.append([name, index,str(year), value])
        #valori idealmente osservato in un dato anno
        f=pred0.iloc[:,yearobs==year].values
        daystoobsY=daystoobs[yearobs==year]
        # or Event?
        if PP or Event:
            ft=pred0t.iloc[:,yeart==year].values
            if PP:
                ftt=YY[yeart==year]
                if (yeart==year).sum()< 3: 
                    Rest.append((np.nan,np.nan,np.nan))
                    continue
                M=ft.mean(axis=1)
                Mm=M.mean()
                STD=ft.std(axis=1)
                STDm=STD.mean()
                argMAX=ft.argmax(axis=1)
                argMAXm=argMAX.mean()
                rest=((np.abs(M-Mm)>=np.abs(ftt.mean()-Mm)).sum(),
                      (np.abs(STD-STDm)>=np.abs(ftt.std()-STDm)).sum(),(np.abs(argMAX-argMAXm)>=np.abs(ftt.argmax()-argMAXm)).sum())
                Rest.append(rest)
            if Event:
                #in caso di due sole osservazioni in un anno meglio saltare
                if (yeart==year).sum()< 3: 
                    Loadingres(name="Mevent",index="mean",value=np.na)
                    Loadingres(name="Mevent",index="sd",value=np.na)
                    continue
                tt=t[yeart==year]
                stat=np.diff(ft,axis=1)/np.diff(tt)
                try:
                    temp=(-stat).argmax(axis=1)
                except ValueError:
                    pass
                #average between the two date for which difference was estimated
                eventf=(tt[temp]+tt[temp+1])/2
                Loadingres(name="Mevent",index="mean",value=eventf.mean()-(year*freq))
                Loadingres(name="Mevent",index="sd",value=eventf.std())
            
        if Event:
            stat=np.diff(f,axis=1)
            temp=(-stat).argmax(axis=1)
            MaxTemp=(-stat).max(axis=1)
            eventf=(daystoobsY[temp]+daystoobsY[temp+1])/2
            Loadingres(name="event",index="mean",value=eventf.mean()-(year*freq))
            Loadingres(name="event",index="sd",value=eventf.std())
            Loadingres(name="eventVal",index="mean",value=MaxTemp.mean())
            Loadingres(name="eventVal",index="sd",value=MaxTemp.std())

        ff=f.mean(axis=1)
        meansdyear.append(ff)
        Loadingres(name="mean",index="mean",value=ff.mean())
        Loadingres(name="mean",index="sd",value=ff.std())
        
        ff=f.std(axis=1)
        Loadingres(name="stdintra",index="mean",value=ff.mean())
        Loadingres(name="stdintra",index="sd",value=ff.std())
        varyear+=(ff**2)*f.shape[1]/totdays
        
        ff=daystoobsY[f.argmax(axis=1)]
        ffmaxpos=ff-year*freq
        Loadingres(name="maxpos",index="mean",value=circmean(ffmaxpos, high=365,low=0))
        Loadingres(name="maxpos",index="meanmean",value=ffmaxpos.mean())
        Loadingres(name="maxpos",index="sd",value=circvar(ffmaxpos, high=365,low=0)**0.5)
        
        if minmaxFeat:
            fval=f.max(axis=1)
            Loadingres(name="maxposVal",index="mean",value=fval.mean())
            Loadingres(name="maxposVal",index="sd",value=fval.std())
            
            ffmin=daystoobsY[(-f).argmax(axis=1)]
            ffminpos=ffmin-year*freq
            Loadingres(name="minpos",index="mean",value=circmean(ffminpos, high=365,low=0))
            Loadingres(name="minpos",index="sd",value=circvar(ffminpos, high=365,low=0)**0.5)

            fval=f.min(axis=1)
            Loadingres(name="minposVal",index="mean",value=fval.mean())
            Loadingres(name="minposVal",index="sd",value=fval.std())
    
    varinter=(np.var(pred0.values,axis=1)-varyear)**0.5

    Loadingres(name="sdinter",index="mean",year="all",value=varinter.mean())
    Loadingres(name="sdinter",index="sd",year="all",value=varinter.std())

    meansdinter=(np.concatenate([x[np.newaxis] for x in meansdyear])).std(axis=1)
    Loadingres(name="sdinterm",index="mean",year="all",value=meansdinter.mean())
    Loadingres(name="sdinterm",index="sd",year="all",value=meansdinter.std())
    Res=pd.DataFrame(Res)
    Res.columns=["stat","index","year","value"]
    Res=Res.set_index(["stat","index","year"])["value"]
    if PP:
        Rest=pd.DataFrame(Rest,columns=["mean","stdintra","maxPos"], index=np.arange(Years).astype("U2"))
        Rest=Rest/reps
    return Res, Rest

def FeatureBayes2(modelA, obsdays, obsvalue, freq=365,PP=True, Event=False,reps=100,freqsubsample=5, YearsL=None, offset=0, minmaxFeat=False):
    """
    Extract 4 descriptors of the time series from a model or multiyear model and perfrom a PP test of goodness of fit 
    """
    t=obsdays
    YY=obsvalue
    Years=np.ceil(t.max()/freq).astype("int")
    freqsubsample=freqsubsample
    #daystoobs=np.arange(Years*freq)
    daystoobs=np.arange(0,Years*freq, freqsubsample)
    X=modelA._BayesianLinearModel__basis(daystoobs)
    #pred0=pd.DataFrame(np.array([np.nansum(X*Randomizer(modelA,daystoobs),axis=1) for x in np.arange(reps)]))
    #params=modelA.random(samples=1000)
    #print(X.shape,Years,freqsubsample)
    try:
        params=Randomizer(modelA,daystoobs,reps)
        #modelA._predict_mean
    except ValueError:
        print(obsvalue)
    pred0=pd.DataFrame(np.nansum(X[:,np.newaxis,:]*params,axis=2).T)
    Res=[]
    Rest=[]
    varyear=0
    meansdyear=[]
    totdays=pred0.shape[1]
    if PP or Event:
        X=modelA._BayesianLinearModel__basis(t)
        #pred0t=pd.DataFrame(np.array([np.nansum(X*Randomizer(modelA,t),axis=1) for x in np.arange(reps)]))
        paramst=Randomizer(modelA,t,reps)
        pred0t=pd.DataFrame(np.nansum(X[:,np.newaxis,:]*paramst,axis=2).T)
        #pred0t=pred0.iloc[:,t-1]
        yeart=((t+offset)/freq).astype("int")
    yearobs=((daystoobs+offset)/freq).astype("int")
    for year in np.arange(Years):
        res=[]
        #f=pred0.iloc[:,year*freq:(year+1)*freq].values
        f=pred0.iloc[:,yearobs==year].values
        daystoobsY=daystoobs[yearobs==year]
        # or Event?
        if PP or Event:
            ft=pred0t.iloc[:,yeart==year].values
            if PP:
                ftt=YY[yeart==year]
                if (yeart==year).sum()< 3: 
                    Rest.append((np.nan,np.nan,np.nan))
                    continue
                M=ft.mean(axis=1)
                Mm=M.mean()
                STD=ft.std(axis=1)
                STDm=STD.mean()
                argMAX=ft.argmax(axis=1)
                argMAXm=argMAX.mean()
                rest=((np.abs(M-Mm)>=np.abs(ftt.mean()-Mm)).sum(),
                      (np.abs(STD-STDm)>=np.abs(ftt.std()-STDm)).sum(),(np.abs(argMAX-argMAXm)>=np.abs(ftt.argmax()-argMAXm)).sum())
                Rest.append(rest)
            if Event:
                #in caso di due sole osservazioni in un anno meglio saltare
                if (yeart==year).sum()< 3: 
                    res.append(pd.Series([np.nan,np.nan], name="event", index=["mean","sd"]))
                    continue
                tt=t[yeart==year]
                stat=np.diff(ft,axis=1)/np.diff(tt)
                try:
                    temp=(-stat).argmax(axis=1)
                except ValueError:
                    #print(stat)
                    #print(ft)
                    #print(t,YY)
                    pass
                #print(daystoobsY.shape,temp.max())
                #average between the two date for which difference was estimated
                eventf=(tt[temp]+tt[temp+1])/2
                res.append(pd.Series([eventf.mean()-(year*freq),eventf.std()], name="Mevent", index=["mean","sd"]))
            
        stat=np.diff(f,axis=1)
        temp=(-stat).argmax(axis=1)
        MaxTemp=(-stat).max(axis=1)
        if Event:
            eventf=(daystoobsY[temp]+daystoobsY[temp+1])/2
            res.append(pd.Series([eventf.mean()-(year*freq),eventf.std()], name="event", index=["mean","sd"]))
            res.append(pd.Series([MaxTemp.mean(),MaxTemp.std()], name="eventVal", index=["mean","sd"]))
        ff=f.mean(axis=1)
        meansdyear.append(ff)
        res.append(pd.Series([ff.mean(),ff.std()], name="mean", index=["mean","sd"]))
        ff=f.std(axis=1)
        varyear+=(ff**2)*f.shape[1]/totdays
        res.append(pd.Series([ff.mean(),ff.std()], name="stdintra", index=["mean","sd"]))
        #res.append(pd.Series([f.idxmax(axis=1).mean(),f.idxmax(axis=1).std()], name="maxpos", index=["mean","sd"]))
        ff=daystoobsY[f.argmax(axis=1)]
        ffmaxpos=ff-year*freq
        #provo con media circolare
        res.append(pd.Series([circmean(ffmaxpos, high=365,low=0),circvar(ffmaxpos, high=365,low=0)**0.5], name="maxpos", index=["mean","sd"]))
        #correzione grezza per na che appaiano nel mucchio
        res[-1][pd.isnull]=-1
        ##forse sostituire mean con median per dare risultati sensati in casi di più di un picco.
        #res.app end(pd.Series([ffmaxpos.mean(),ffmaxpos.std()], name="maxpos", index=["mean","sd"]))
        if minmaxFeat:
            fval=f.max(axis=1)
            res.append(pd.Series([fval.mean(),fval.std()], name="maxposVal", index=["mean","sd"]))
            ffmin=daystoobsY[(-f).argmax(axis=1)]
            ffminpos=ffmin-year*freq
            #provo con media circolare
            res.append(pd.Series([circmean(ffminpos, high=365,low=0),circvar(ffminpos, high=365,low=0)**0.5], name="minpos", index=["mean","sd"]))
            ##forse sostituire mean con median per dare risultati sensati in casi di più di un picco.
            #res.append(pd.Series([ffmaxpos.mean(),ffmaxpos.std()], name="maxpos", index=["mean","sd"]))
            fval=f.min(axis=1)
            res.append(pd.Series([fval.mean(),fval.std()], name="minposVal", index=["mean","sd"]))
        Res.append(pd.DataFrame(res).stack())
        
    Res=pd.concat(Res, keys=np.arange(Years).astype("U2"), axis=1)
    varinter=(np.var(pred0.values,axis=1)-varyear)**0.5
    #Res.loc[("sdinter","mean"),str(Years)]=varinter.mean()
    #Res.loc[("sdinter","sd"),str(Years)]=varinter.std()
    Res.loc[("sdinter","mean"),"all"]=varinter.mean()
    Res.loc[("sdinter","sd"),"all"]=varinter.std()
    #meansdinter=np.concatenate(meansdinter)
    meansdinter=(np.concatenate([x[np.newaxis] for x in meansdyear])).std(axis=1)
    Res.loc[("sdinterm","mean"),"all"]=meansdinter.mean()
    Res.loc[("sdinterm","sd"),"all"]=meansdinter.std()
    for y in np.arange(Years+1,YearsL+1):
        Res.loc[:,str(y)]=-1
        Res.loc["sdinter",str(y)]=np.nan
        Res.loc["sdinterm",str(y)]=np.nan
    if PP:
        Rest=pd.DataFrame(Rest,columns=["mean","stdintra","maxPos"], index=np.arange(Years).astype("U2"))
        Rest=Rest/reps
    return Res, Rest

def weighted_quantile(values, sample_weight,quantiles, 
                      values_sorted=False, old_style=False):
    """ Very close to numpy.percentile, but supports weights.
    NOTE: quantiles should be in [0, 1]!
    :param values: numpy.array with data
    :param quantiles: array-like with many quantiles needed
    :param sample_weight: array-like of the same length as `array`
    :param values_sorted: bool, if True, then will avoid sorting of
        initial array
    :param old_style: if True, will correct output to be consistent
        with numpy.percentile.
    :return: numpy.array with computed quantiles.
    """
    #print(values.shape,sample_weight.shape)
    values=values.flatten()
    sample_weight=sample_weight.flatten()
    #print(values.shape,sample_weight.shape)
    finvalue=np.isfinite(values)
    finweigth=np.isfinite(sample_weight)
    values=values[finvalue&finweigth]
    sample_weight=sample_weight[finvalue&finweigth]
    values = np.array(values)
    quantiles = np.array(quantiles)
    if sample_weight is None:
        sample_weight = np.ones(len(values))
    sample_weight = np.array(sample_weight)
    assert np.all(quantiles >= 0) and np.all(quantiles <= 1), \
        'quantiles should be in [0, 1]'

    if not values_sorted:
        sorter = np.argsort(values)
        values = values[sorter]
        sample_weight = sample_weight[sorter]

    weighted_quantiles = np.cumsum(sample_weight) - 0.5 * sample_weight
    if old_style:
        # To be convenient with numpy.percentile
        weighted_quantiles -= weighted_quantiles[0]
        weighted_quantiles /= weighted_quantiles[-1]
    else:
        weighted_quantiles /= np.sum(sample_weight)
    return np.interp(quantiles, weighted_quantiles, values)

class HarmBayMod(object):
    """
    Class that handle runs of the model
    """
    def __init__(self, xarrayD=None , reformat_dict={"maxvalues":1,"nan":np.inf}, geotransform=True, bandNameFormatter=pd.to_datetime):
        """
        initialize: it assume that first dimension is time, then latitude, finally longitude
        """
        self.PP=False
        self.Event=False
        self.Expected=False
        self.change=False
        self.memoryR=False
        self.variableLarge=["maxpos","event"]
        self.RawData=xarrayD
        self.reformat_dict=reformat_dict
        if self.RawData is None:
            return 
        assert self.RawData.dims==('Time', 'Y', 'X'), "Actual Dimension "+str(self.RawData.dims) +" do not match expectation:('Time', 'Y', 'X')"
        
        if geotransform:
            #Transform, maybe not relevant, check only if there are rotation
            self.nt,self.ny, self.nx=self.RawData.shape
            if "transform" in self.RawData.attrs:
                self.transform=self.RawData.attrs["transform"]
                if len(self.transform)==9:
                    self.transform=self.transform[:6]
                self.dx,mx,self.xmin,my,self.dy,self.ymax=self.transform
            else:
                mx=0
                my=0
            temp=self.RawData.X.diff(dim="X").values
            if not all((temp[0]-temp)==0):
                print("X dimension is not regular")
            self.dx=temp[0]
            temp=self.RawData.Y.diff(dim="Y").values
            if not all((temp[0]-temp)==0):
                print("Y dimension is not regular")
            self.dy=temp[0]
            self.xmin=self.RawData.X.values.min()-(self.dx/2)
            self.ymax=self.RawData.Y.values.max()-(self.dy/2)
            self.transform=np.array([self.dx,mx,self.xmin,my,self.dy,self.ymax])
            assert (mx==my)&(mx==0), "no rotation expected"
            #CRS
            #double conversion to handle wkt that have missing point in the non geometric description
            crs=None
            try:
                for crs_attr in [k  for k in self.RawData.attrs.keys() if k in ["spatial_ref","coordinate_system_string","crs_wkt","crs"]]:
                    crs=CRS(CRS(self.RawData.attrs[crs_attr]).to_proj4())
                    break
                if (not crs):
                    if ("crs" in self.RawData.coords):
                        for crs_attr in [k  for k in self.RawData.crs.attrs.keys() if k in ["spatial_ref","coordinate_system_string","crs_wkt"]]:
                            crs=CRS(CRS(self.RawData.crs.attrs[crs_attr]).to_proj4()) # Fixed Saverio's typo here
                            break
            except Exception:
                crs = None
                
            # Guard against the NoneType crash
            if crs is not None:
                self.crs=crs.to_wkt()
                self.RawData.attrs["pyproj_srs"]="epsg:"+str(crs.to_epsg())
            else:
                self.crs="Local_Cartesian"
                self.RawData.attrs["pyproj_srs"]="epsg:27700" # Safe fallback dummy value
        #Calculating time difference
        self.start=self.RawData[:,0,0].Time[0].values.astype('datetime64[Y]')
        if "band_names" not in self.RawData.attrs:
            self.tff=(self.RawData[:,0,0].Time-self.start).values.astype('timedelta64[D]').astype("int")
        else:
            time=bandNameFormatter([x.strip() for x in self.RawData.attrs['band_names'].split(",")])
            self.start=time[0:1].values[0].astype("datetime64[Y]")
            self.tff=(time-self.start).values.astype('timedelta64[D]').astype("int")
        self.RawData.coords["tff"]=("Time",self.tff)
    
    def Load(self,filename):
        self.Data=xr.open_dataset(filename)
        
    def RunStackBayesian(self, **args):
        """
        Run over a cube
        """
        #RunStackBayesian(self, memoryR=False, PP=False, Event=False, change=True, reps=100, Expected=False, dense=False,freqsubsample=15):
        for key in args:
            self.__dict__[key]=args[key]
        #self.PP=PP
        #self.Event=Event
        #self.change=change
        #self.Expected=Expected
        #self.dense=dense
        #self.freqsubsample=freqsubsample
        #xarrayD=self.RawData
        
        #tff=(MSAVIOLD[:,0,0].band-np.array([startdate]).astype('datetime64[ns]')).values.astype('timedelta64[D]').astype("int")
        dummy=RunnerBayesian(values=np.random.random(len(self.tff)),times=self.tff, name=False, PP=self.PP, change=self.change,Event=self.Event, Expected=self.Expected, dense=self.dense)
        print("dummy")
        print(dummy.shape[0])
        
        self.Data=np.apply_along_axis(arr=self.RawData.values, axis=0, func1d=RunnerBayesian, times=self.tff, change=self.change,PP=self.PP, Event=self.Event,reformat_dict=self.reformat_dict, reps=self.reps, Expected=self.Expected, Loutput=dummy.shape[0], dense=self.dense)
        if self.memoryR:
            self.MemoryRelease()
        return self.Data
    def DaskBayesian3(self, values,Loutput):
        #DaskBayesian3(self, values, times, reformat_dict,fulloutput,name, Event,Expected, change, reps, Loutput,freqsubsample, dense):
        return np.apply_along_axis(arr=values, axis=-1, func1d=RunnerBayesian, times=self.tff,reformat_dict=self.reformat_dict, 
                                    Loutput=Loutput, Event=self.Event,Expected=self.Expected,change=self.change, reps=self.reps,freqsubsample=self.freqsubsample, dense=self.dense)
    
    def DaskRunstackBayesian(self,**args):
        #DaskRunstackBayesian(self,PP=False, Event=False, change=True, reps=100, freqsubsample=15, Expected=False, dense=True)
        for key in args:
            self.__dict__[key]=args[key]
        #self.PP=PP
        #self.Event=Event
        #self.change=change
        #self.Expected=Expected
        #self.freqsubsample=freqsubsample
        #self.dense=dense
        dummy=RunnerBayesian(values=np.random.random(len(self.tff)),times=self.tff, name=False, PP=self.PP, change=self.change, Event=self.Event, Expected=self.Expected, freqsubsample=self.freqsubsample, dense=self.dense)
        self.Data=xr.apply_ufunc(self.DaskBayesian3,
                  self.RawData,
                  kwargs={ "Loutput":dummy.shape[0]},
                  input_core_dims=[["Time"]],output_core_dims=[["stat"]],
                  dask='parallelized',output_dtypes=[self.RawData.dtype],
                  output_sizes={"stat":dummy.shape[0]}).transpose("stat","Y","X")
        #self.Data.to_netcdf("Temp.nc")
    def MemoryRelease(self):
        """
        erase data to release memory after run to alleviate 
        """
        del self.RawData
    def FormatBayes(self, addWGS84=False, empty=False):
        """
        Build a coherent netcdf Dataset  with multiple variables from the single array of output
        """
        Name=RunnerBayesian(values=np.random.random(len(self.tff)),times=self.tff, name=True, PP=self.PP, change=self.change,Event=self.Event, Expected=self.Expected)[1]
        tff=self.tff
        if not self.Expected:
            start=self.start.astype("int")+1970
            var,stat,Year=map(pd.Series,zip(*[x.split("_") for x in Name]))
            ClassName=var+"_"+stat
        else:
            start=self.start.astype("datetime64[D]")
            ClassName,Year=map(pd.Series,zip(*[x.split("_") for x in Name]))
        
        if not empty:
            BDFBayesx=self.Data
        else:
            #BDFBayesx=np.zeros((len(ClassName),*self.RawData.shape[1:]))
            BDFBayesx=np.memmap("tempo.mm",dtype='float32', mode='w+',shape=(len(ClassName),*self.RawData.shape[1:]))
            BDFBayesx[:]=0
            
        if BDFBayesx.__class__!=xr.DataArray([2]).__class__:
            BDFBayesx=xarray.DataArray(BDFBayesx, dims=["stat","Y","X"])
            #BDFBayesx=BDFBayesx.rename(dim_0="Year",dim_1="y",dim_2="x")
            #produce random data to extract name
            
            #Build grid of affine trasform
            X=np.arange(self.nx)*self.dx+self.xmin+self.dx/2
            Y=np.arange(self.ny)*self.dy+self.ymax+self.dy/2
            #Annotate raw output
            BDFBayesx.coords["Y"]=(("Y"),Y)
            BDFBayesx.coords["X"]=(("X"),X)
        BDFBayesx.coords["X"].attrs["long_name"]="Easting"
        BDFBayesx.coords["Y"].attrs["long_name"]="Northing"
        BDFBayesx.coords["Time"]=(("stat"),Year)
        #Time=list(set(Year.astype("int")))
        #Time.sort()
        #BDFBa=xr.Dataset(data_vars={"X":X,"Y":Y,"Time":Time})
        BDFBa=xr.Dataset()
        for Class in [ x for x in ClassName if x not in ["sdinter_mean","sdinter_sd","sdinterm_mean","sdinterm_sd"]]:
            #BDFBa[Class]=BDFBayesx.where(BDFBayesx.ClassName==Class, drop=True)
            #BDFBa[Class]=BDFBayesx[BDFBayesx.ClassName==Class,:,:]
            #BDFBa[Class]=(("Time","Y","X"),np.memmap("tempo_"+Class+".mm", dtype='float32', mode='w+', shape=(len(set(Year)),len(Y),len(X))))
            BDFBa[Class]=BDFBayesx.sel({"stat":ClassName==Class})
            BDFBa[Class].attrs["grid_mapping"]="crs"
        #del BDFBa.coords["ClassName"]
        #BDFBa.coords["stat"]=(("stat"),Year)
        #Add variables with no time dimension
        for Class in [ x for x in ClassName if x in ["sdinter_mean","sdinter_sd","sdinterm_mean","sdinterm_sd"]]:
            #BDFBa[Class]=(("Y","X"),np.memmap("tempo_"+Class+".mm", dtype='float32', mode='w+', shape=(len(Y),len(X))))
            #BDFBa[Class].values[:]=(("Y","X"),BDFBayesx.sel({"stat":ClassName==Class})[0,:,:].values)
            BDFBa[Class]=(("Y","X"),BDFBayesx.sel({"stat":ClassName==Class})[0,:,:])
            BDFBa[Class].attrs["grid_mapping"]="crs"
        BDFBa=BDFBa.rename({"stat":"Time"})
        BDFBa.coords["Time"]=BDFBa.coords["Time"].values.astype("int")+start
        #print(start)
        self.Data=BDFBa
        BDFBa.coords["crs"]=CRS(CRS(self.crs).to_epsg()).name
        BDFBa.coords["crs"].attrs["crs_wkt"]=self.crs
        BDFBa.coords["crs"].attrs["spatial_ref"]=self.crs
        BDFBa.coords["crs"].attrs["transform"]=self.transform
        BDFBa.attrs["spatial_ref"]=self.crs
        BDFBa.attrs["pyproj_srs"]="epsg:"+str(CRS(self.crs).to_epsg())
        BDFBa.attrs["transform"]=self.transform
        BDFBa.attrs["GeoTransform"]=self.transform
        BDFBa.attrs["Conventions"]="CF-1.6"
        if addWGS84:
            Xsq,Ysq=np.meshgrid(BDFBa.X,BDFBa.Y)
            LAT,LONG=transform(CRS(self.crs),CRS("epsg:4326"),x=Xsq,y=Ysq)
            BDFBa.coords["Lat"]=(("Y","X"),LAT)
            BDFBa.coords["Long"]=(("Y","X"),LONG)
            BDFBa.coords["Long"].attrs["long_name"]="Longitude"
            BDFBa.coords["Lat"].attrs["long_name"]="Latitude"
            BDFBa.coords["Lat"].attrs["units"]="degrees_north"
            BDFBa.coords["Long"].attrs["units"]="degrees_east"

    def Store(self, namefile, default=True):
        enc={}
        MMaxs={
        "Expected":1,
        "ExpectedSD":1,
        "mean_mean":0.8,
        "maxposVal_mean":0.8,
        "minposVal_mean":0.8,
        "maxposVal_sd":0.3,
        "minposVal_sd":0.3,
        "mean_sd":0.3,
        "stdintra_mean":0.8,
        "stdintra_sd":0.3,
        "maxpos_mean":365,
        "maxpos_sd":180,
        "minpos_mean":365,
        "minpos_sd":180,
        "sdinter_mean":0.3,
        "sdinter_sd":0.3,
        "sdinterm_mean":0.3,
        "sdinterm_sd":0.3
        }
        for x in self.Data.data_vars.keys():
            Max=self.Data[x].max().values
            Min=self.Data[x].min().values
            MMax=np.max([np.abs(Min), np.abs(Max)])
            enc[x]={'_FillValue': -9999,'dtype': 'int16'}
            if default:
                enc[x]["scale_factor"]=MMax/10000
            else:
                enc[x]["scale_factor"]=MMaxs[x]/10000
        self.Data.to_netcdf(namefile+".nc", encoding=enc)
        
    def OLDFormatBayes(self, namefile, addWGS84=False):
        """
        Build a coherent netcdf Dataset  with multiple variables from the single array of output
        """
        tff=self.tff
        if not self.Expected:
            start=self.start.astype("int")+1970
        else:
            start=self.start.astype("datetime64[D]")
        BDFBayesx=self.Data
        if self.Data.__class__!=xr.DataArray([2]).__class__:
            BDFBayesx=xarray.DataArray(self.Data, dims=["stat","Y","X"])
            #BDFBayesx=BDFBayesx.rename(dim_0="Year",dim_1="y",dim_2="x")
            #produce random data to extract name
            
            #Build grid of affine trasform
            X=np.arange(self.nx)*self.dx+self.xmin+self.dx/2
            Y=np.arange(self.ny)*self.dy+self.ymax+self.dy/2
            #Annotate raw output
            BDFBayesx.coords["Y"]=(("Y"),Y)
            BDFBayesx.coords["X"]=(("X"),X)
        Name=RunnerBayesian(values=np.random.random(len(self.tff)),times=self.tff, name=True, PP=self.PP, change=self.change,Event=self.Event, Expected=self.Expected)[1]
        BDFBayesx.coords["X"].attrs["long_name"]="Easting"
        BDFBayesx.coords["Y"].attrs["long_name"]="Northing"
        #Year=[ int(x.split("_")[-1])+start for x in Name]
        Year=[ x.split("_")[-1] for x in Name]
        BDFBayesx.coords["stat"]=(("stat"),Year)
        ClassName=[ "_".join(x.split("_")[:-1]) for x in Name]
        BDFBayesx.coords["ClassName"]=(("stat"),ClassName)
        #Split raw output in the variables to build composite netcdf
        BDFBa=xarray.Dataset()
        
        for Class in [ x for x in BDFBayesx.ClassName if x not in ["sdinter_mean","sdinter_sd","sdinterm_mean","sdinterm_sd"]]:
            Class=str(Class.values)    
            #BDFBa[Class]=BDFBayesx.where(BDFBayesx.ClassName==Class, drop=True)
            #BDFBa[Class]=BDFBayesx[BDFBayesx.ClassName==Class,:,:]
            BDFBa[Class]=BDFBayesx.sel({"stat":BDFBayesx.ClassName==Class})
            BDFBa[Class].attrs["grid_mapping"]="crs"
        del BDFBa.coords["ClassName"]
        #format year, once removed stat not related to time
        Year=[ int(x.split("_")[-1])+start for x in BDFBa.stat.values]
        BDFBa.coords["stat"]=(("stat"),Year)
        #Add variables with no time dimension
        for Class in [ x for x in ClassName if x in ["sdinter_mean","sdinter_sd","sdinterm_mean","sdinterm_sd"]]:
            BDFBa[Class]=(("Y","X"),BDFBayesx.sel({"stat":BDFBayesx.ClassName==Class})[0,:,:])
            BDFBa[Class].attrs["grid_mapping"]="crs"
        BDFBa=BDFBa.rename({"stat":"Time"})
        BDFBa.coords["Time"]=BDFBa.coords["Time"].values.astype("int")+start
        print(start)
        #for notime in ["sdinter","sdinterm"]:
        #    for stat in ["mean","sd"]:
        #        BDFBa[notime+"_"+stat]=(("Y","X"),BDFBayesx[BDFBayesx.ClassName==notime+"_"+stat,:,:][0,:,:])
        #        BDFBa[notime+"_"+stat].attrs["grid_mapping"]="crs"
        BDFBa.coords["crs"]=0
        BDFBa.coords["crs"].attrs["crs_wkt"]=self.crs
        BDFBa.coords["crs"].attrs["spatial_ref"]=self.crs
        BDFBa.coords["crs"].attrs["transform"]=self.transform
        BDFBa.attrs["spatial_ref"]=self.crs
        BDFBa.attrs["pyproj_srs"]="epsg:"+str(CRS(self.crs).to_epsg())
        BDFBa.attrs["transform"]=self.transform
        BDFBa.attrs["GeoTransform"]=self.transform
        BDFBa.attrs["Conventions"]="CF-1.6"
        if addWGS84:
            Xsq,Ysq=np.meshgrid(BDFBa.X,BDFBa.Y)
            LAT,LONG=transform(CRS(self.crs),CRS("epsg:4326"),x=Xsq,y=Ysq)
            BDFBa.coords["Lat"]=(("Y","X"),LAT)
            BDFBa.coords["Long"]=(("Y","X"),LONG)
            BDFBa.coords["Long"].attrs["long_name"]="Longitude"
            BDFBa.coords["Lat"].attrs["long_name"]="Latitude"
            BDFBa.coords["Lat"].attrs["units"]="degrees_north"
            BDFBa.coords["Long"].attrs["units"]="degrees_east"
        #return BDFBa
        #encode and save
        #enc={}
        #enc.update([ (x,{'dtype': 'int16', 'scale_factor': BDFBa[x].max().round().astype("int").values/10000, '_FillValue': -9999}) for x in BDFBa.data_vars.keys()])
        enc={}
        for x in BDFBa.data_vars.keys():
            Max=BDFBa[x].max().values
            Min=BDFBa[x].min().values
            MMax=np.max([np.abs(Min), np.abs(Max)])
            enc[x]={'_FillValue': -9999,'dtype': 'int16'}
            enc[x]["scale_factor"]=MMax/10000.
        #for var in self.variableLarge:
        #for var in ["maxpos","event"]:
        #    enc[var+'_mean']['scale_factor']=1
        #    enc[var+'_sd']['scale_factor']=1
        self.Data=BDFBa
        BDFBa.to_netcdf(namefile+".nc", encoding=enc)
        #return BDFBa
    def ExpectPlot(self,sel=None, isel=None, mask=None,quantiles=[0.25,0.75],fontsize=20):
        da=self.Data.isel(isel).sel(sel)
        if mask is not None :
            da=da.where(mask)
        from matplotlib import pyplot as plt
        SumSD=da["ExpectedSD"].reduce(np.nansum, dim=["Y","X"])
        W=da["ExpectedSD"]/SumSD
        res=[]
        for i in np.arange(da.Time.shape[0]):
            res.append(weighted_quantile(da["Expected"].values[i,:],W.values[i],quantiles=quantiles)[np.newaxis,:])
        Q=xr.DataArray(np.concatenate(res),dims=["Time","quantile"], coords=[da.Time,quantiles ])
        Q.plot(x="Time", hue="quantile", size=10, aspect=3)
        ax=plt.gca()
        (da["Expected"]*W).reduce(np.nansum, dim=["Y","X"]).plot(label="mean",ax=ax)
        seasons=["-03-21","-06-21","-09-21","-12-21"]
        year="-01-01"
        startyear=da.Time.min().values.astype("datetime64[Y]").astype("str").astype("int")
        endyear=da.Time.max().values.astype("datetime64[Y]").astype("str").astype("int")
        Seasons=np.array(sum([[str(x)+y for y in seasons] for x in np.arange(startyear,endyear+1)],[])).astype("datetime64[D]")
        years=np.array([str(x)+year for x in np.arange(startyear,endyear+1)]).astype("datetime64[D]")
        MIN=Q.min().item()
        MAX=Q.max().item()
        plt.vlines(Seasons, MIN,MAX, linestyles="dashed",color="pink")
        plt.vlines(years, MIN,MAX,color="black")
        ax=plt.gca()
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(fontsize)
        #ax.legend()
        return Q
    def ObservedPlot(self,sel=None, isel=None, mask=None,quantiles=[0.25,0.75],fontsize=20):
        da=self.RawData.isel(isel).sel(sel)
        if mask is not None :
            da=da.where(mask)
        from matplotlib import pyplot as plt
        Test=da.copy()
        temp=Test.values.astype("float")
        temp[temp==self.reformat_dict["nan"]]=np.nan
        temp=temp/self.reformat_dict["maxvalues"]
        Test.values=temp
        Test.coords["Time"]=("Time",self.start.astype("datetime64[D]")+self.tff)
        Test.mean( dim=["Y","X"]).plot(label="mean",size=10, aspect=3)
        ax=plt.gca()
        Test=Test.quantile( dim=["Y","X"], q=[0.25,0.75])
        Test.plot(hue="quantile",ax=ax)

        seasons=["-03-21","-06-21","-09-21","-12-21"]
        year="-01-01"
        startyear=Test.Time.min().values.astype("datetime64[Y]").astype("str").astype("int")
        endyear=Test.Time.max().values.astype("datetime64[Y]").astype("str").astype("int")
        Seasons=np.array(sum([[str(x)+y for y in seasons] for x in np.arange(startyear,endyear+1)],[])).astype("datetime64[D]")
        years=np.array([str(x)+year for x in np.arange(startyear,endyear+1)]).astype("datetime64[D]")
        MIN=float(Test.min())
        MAX=float(Test.max())
        plt.vlines(Seasons, MIN,MAX, linestyles="dashed",color="pink")
        plt.vlines(years, MIN,MAX,color="black")
        ax=plt.gca()
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(20)
        
    def TimeTrendPlot(self,sel=None, isel=None, mask=None,quantiles=[0.25,0.5,0.75],stat="mean", fontsize=20, size=8, aspect=3):
        from matplotlib import pyplot as plt
        import seaborn as sns
        da=self.Data.isel(isel).sel(sel)
        if mask is not None :
            da=da.where(mask)
        SumSD=da[stat+"_sd"].reduce(np.nansum, dim=["Y","X"])
        W=(1/(da[stat+"_sd"]**2))/SumSD
        
        res=[]
        for i in np.arange(da.Time.shape[0]):
            res.append(weighted_quantile(da[stat+"_mean"].values[i,:],W.values[i],quantiles=quantiles)[np.newaxis,:])
        Q=xr.DataArray(np.concatenate(res),dims=["Time","quantile"], coords=[da.Time,quantiles ])
        Q=Q.to_dataframe(stat).reset_index()
        xx=sns.color_palette("icefire",n_colors=len(quantiles))
        sns.relplot(x="Time",y=stat, hue="quantile", size=size, aspect=aspect,data=Q,kind="line", palette=xx)
        ax=plt.gca()
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(fontsize)
        return Q
        
    def SinglePixelPlot(self,isel=None, sel=None,reformat_dict={"nan":-32768, "maxvalues":10000}, change=True, sizeplot=10, aspect=4,fontsize=20):
        from matplotlib import pyplot as plt
        values=self.RawData.isel(isel).sel(sel).squeeze()
        times=values.tff
        Mod, RES, rest, evidences, chack=RunnerBayesian(values.values,times.values,reformat_dict=reformat_dict,fulloutput=True, Loutput=4,PP=False, change=change,Event=False, Expected=False)
        MM=Mod.predict(times.values,np.arange(0,100)/150)[2]
        MM=xr.DataArray(MM, dims=["VI","Time"], coords=[ np.arange(0,100)/150,self.start.astype("datetime64[D]")+times.values])
        temp=values.astype("float")
        temp[temp==reformat_dict["nan"]]=np.nan
        temp=temp/reformat_dict["maxvalues"]
        
        MM.plot(size=sizeplot, aspect=aspect)
        plt.plot(self.start.astype("datetime64[D]")+values.tff,Mod.predict(times.values), color="red")
        plt.plot(self.start.astype("datetime64[D]")+values.tff,temp, color="blue", marker="o")
        
        seasons=["-03-21","-06-21","-09-21","-12-21"]
        year="-01-01"
        startyear=(self.start.astype("datetime64[D]")+times[0]).values.astype("datetime64[Y]").astype("str").astype("int")
        endyear=(self.start.astype("datetime64[D]")+times[-1]).values.astype("datetime64[Y]").astype("str").astype("int")
        Seasons=np.array(sum([[str(x)+y for y in seasons] for x in np.arange(startyear,endyear)],[])).astype("datetime64[D]")
        years=np.array([str(x)+year for x in np.arange(startyear,endyear)]).astype("datetime64[D]")
        MIN=0
        MAX=100/150
        plt.vlines(Seasons, MIN,MAX, linestyles="dashed",color="green")
        plt.vlines(years, MIN,MAX,color="white")
        plt.title(str(isel)+str(sel))
        ax=plt.gca()
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(fontsize)
        return Mod
        



class ModelList(object):
    """
    Class of multiyear model it simulate basic feature of Bayesian linear model class
    """
    def __init__(self, modellist=[], freq=365, offset=0):
        self.offset=offset
        self.freq=freq
        self.modellist=modellist
        self.maxparam=max([x.location.shape[0] for x in modellist]+[0])
    def predict(self,X, y=None, variance=False):
        """
        Calculate posterior predictive values.
        """
        res={}
        yeart=((X+self.offset)/self.freq).astype("int").flatten()
        for year, model in enumerate(self.modellist):
            out=model.predict(X[yeart==year][:,np.newaxis], y, variance)
            if out.__class__==np.array([]).__class__:
                #print(out.shape)
                out=[out]
            for n,o in enumerate(out):
                res.setdefault(n,[]).append(o)
        #print(res)
        Res=[]
        axis=0
        for r in res:
            if (r==2): axis=1
            Res.append(np.concatenate(res[r],axis=axis))
        if len(Res)==1:
            Res=Res[0]
        return Res
            
    def random(self, obs, samples):
        """
        Sample of parameters from posterior distribution
        """
        res=[]
        yeart=((obs+self.offset)/self.freq).astype("int")
        for year, model  in enumerate(self.modellist):
            temp=model.random(samples)
            if temp.shape[1]<self.maxparam:
                resn=np.array([[np.nan]*self.maxparam]*samples)
                resn[:,:temp.shape[1]]=temp
                temp=resn
            res.append(np.resize(temp,(len(obs[yeart==year]),samples,self.maxparam)))
        #return res
        return np.concatenate(res, axis=0)
    def _BayesianLinearModel__basis(self,obs):
        """
        Expand predictor 
        """
        res=[]
        yeart=((obs+self.offset)/self.freq).astype("int")
        for year,model in enumerate(self.modellist):
            temp=model._BayesianLinearModel__basis(obs[yeart==year][:,np.newaxis])
            if temp.shape[1]<self.maxparam:
                resn=np.zeros((temp.shape[0],self.maxparam))
                resn.fill(np.nan)
                resn[:temp.shape[0],:temp.shape[1]]=temp
                temp=resn
            res.append(temp)
        return np.concatenate(res, axis=0)

def Randomizer(model, nobs, nsamples=1):
    """
    dispatch at the correct class's method
    """
    if model.__class__==ModelList().__class__:
        return model.random(nobs, samples=nsamples)
    else:
        temp=model.random(samples=nsamples)[np.newaxis,:,:]
        return temp

def storefm(BDGBM,ARG):
    Time=BDGBM.Data.Time.shape[0]
    Step=ARG.Step
    #print("Step",Step)
    for var,df in BDGBM.Data.data_vars.items():
        print("out",var,df.shape)
        if (len(df.shape))==2:
            out = np.memmap(ARG.suffix+"_"+var+".mm", dtype='float32', mode='r+', shape=(N,M))
            out[nn:min((nn+Step),N),mm:min((mm+Step),M)]=df
        else:
            out = np.memmap(ARG.suffix+"_"+var+".mm", dtype='float32', mode='r+', shape=(Time,N,M))
            print(nn,N,mm,M,Step)
            print("map", out[:,nn:min((nn+Step),N),mm:min((mm+Step),M)].shape)
            out[:,nn:min((nn+Step),N),mm:min((mm+Step),M)]=df
    # certify that job was completed in stardoutput
    T1=time.time()
    T=str(T1-T0)
    #print("jobitem "+str(n)+" "+str(ARG.n)+" "+T)
    out = np.memmap(ARG.suffix+"_count"+".mm", dtype='float32', mode='r+', shape=(N,M))
    out[nn:min((nn+Step),N),mm:min((mm+Step),M)]=T
    
def storenetcdf(BDGBM,ARG):
    #X=BDGBM.Data.X.shape[0]
    #Y=BDGBM.Data.Y.shape[0]
    #Time=BDGBM.Data.Time.shape[0]
    Step=ARG.Step
    OUT=xr.open_dataset(ARG.suffix+".nc")
    for var in BDGBM.Data.data_vars:
        out=OUT[var]
        if (len(out.shape))==2:
            out[nn:min((nn+Step),N),mm:min((mm+Step),M)]=BDGBM.Data[var]
        else:
            out[:,nn:min((nn+Step),N),mm:min((mm+Step),M)]=BDGBM.Data[var]

    OUT.to_netcdf(ARG.suffix+".nc",mode="a")
    # certify that job was completed in stardoutput
    T1=time.time()
    print("jobitem "+str(n)+" "+str(ARG.n)+" "+str(T1-T0))


if "__main__"==__name__:
    import argparse
    import time
    parser = argparse.ArgumentParser(description='PPresil options', prog="PPresil-BM")
    parser.add_argument("--inputfile", action="store", help="path to input file")
    parser.add_argument("--suffix", dest="suffix",action="store", help="suffix for output")
    parser.add_argument("--step", dest="Step", type=int, action="store", help="size of chunk to analyze")
    parser.add_argument("-n", action="store", type=int, help="n-element of size stepxstep to analyze")
    parser.add_argument("--noffset", action="store", default=0, type=int, help="offset to start counting n-element of size stepxstep to analyze")
    parser.add_argument("--change", action="store_true", help="allow change of seasonality")
    parser.add_argument("--Expected", dest="Expected", action="store_true", help="outout expected value and do not estimate annual phenological statistics")
    parser.add_argument("--Event", dest="Event",action="store_true", help="store the date and value of the maximum change within each year")
    parser.add_argument("--PP", action="store_true" ,help="perform Posterior predictive test [experimental]")
    parser.add_argument("--dense", action="store_true", help="assume that time series is dense and each year can be fitted with its own seasonality model" )
    parser.add_argument("--freqsubsample",type=int, action="store",default=15, help=" frequency of sampling from model in days used in order to estimate phemolgical statistics" )
    parser.add_argument("--reps",default=500, type=int,dest="reps",action="store", help="replicated of model sampling in order to estimate standard variation of prediction")
    parser.add_argument("--dask", action="store_true", help="perform calculation using dask functionality and not numpy")
    parser.add_argument("--netcdf", action="store_true", help="if option n is used, ask for final output as netcdf and not naked array in memory map file")
    parser.add_argument("--minmaxFeat", dest="minmaxFeat",action="store_true", help="add to summary statistics min and max value and min position")
    parser.add_argument("--small", dest="small",action="store_true", help="save on single file even of option n was used")
    parser.add_argument("--bbox", dest="bbox",action="store", help="bounding box in the form west,south,east,north with unit as defined in input file CRS")
    parser.add_argument("--graph", dest="graph",action="store", help="output also a graphic visualization of results")
    parser.add_argument("--mask", dest="mask",action="store", help="mask input with a shapefile")
    
    T0=time.time()
    ARG=parser.parse_args()
    
    #Define BBox
    west=south=east=north=None
    if ARG.bbox:
        west,south,east,north=map(int,ARG.bbox.split(","))
    #Define opener
    if ARG.inputfile.split(".")[-1]=="nc":
        opener=xr.open_dataarray
    else: opener=xr.open_rasterio
    
    # open file with bbox
    ds = xr.open_dataset(ARG.inputfile)
    data_vars = [v for v in ds.data_vars if v not in ['crs', 'spatial_ref']]
    var_name = data_vars[0] 
    
    # 1. Rename dimensions
    A = ds[var_name].rename({"Z1": "Time", "northing": "Y", "easting": "X"})
    
    # 2. Reconstruct the lost phenology dates 
    # This builds a perfect monthly timeline matching your 278 bands
    import pandas as pd
    A['Time'] = pd.date_range(start='2000-01-01', periods=len(A.Time), freq='MS') + pd.Timedelta(days=14)
    
    # 3. Steal the CRS directly from the Shapefile mask
    import geopandas as gpd
    V = gpd.read_file(ARG.mask)
    if hasattr(V.crs, 'to_proj4'):
        A.attrs['proj4'] = V.crs.to_proj4()
    else:
        from pyproj import Proj
        A.attrs['proj4'] = Proj(V.crs).srs
        
    # 4. Clean up the confusing modern tags from R
    if 'grid_mapping' in A.attrs:
        del A.attrs['grid_mapping']
    if 'pyproj_srs' in A.attrs:
        del A.attrs['pyproj_srs']
        
    # A = A.loc[:,north:south,west:east]
    
    #Correct dimension name
    Time, Y, X=A.dims
    DICT={}
    if Time=="band":
        DICT[Time]="Time"
    if X in ["x", "longitude", "Longitude"]:
        DICT[X]="X"
    if Y in ["y", "latitude", "Latitude"]:
        DICT[Y]="Y"
    A=A.rename(DICT)
        #Apply mask
    if ARG.mask:
        #I use HarmBayMod to format CRS and apply mask correctly using salem
        BDGBM=HarmBayMod(A, reformat_dict={"nan":-32768,"maxvalues":10000}, bandNameFormatter=lambda y: pd.Series(y))
        import geopandas as gpd
        import salem
        V=gpd.read_file(ARG.mask)
        VV=V.to_crs(BDGBM.RawData.pyproj_srs)
        BB=VV.bounds
        Z=BDGBM.RawData.salem.grid.region_of_interest(geometry=V.geometry[0])
        BDGBM.RawData.coords["BB"]=(("Y","X"),Z)
        Z=BDGBM.RawData.BB.loc[BB.maxy[0]:BB.miny[0],BB.minx[0]:BB.maxx[0]]
        A=A.loc[:,BB.maxy[0]:BB.miny[0],BB.minx[0]:BB.maxx[0]].where(Z)
    
    #chunking if needed
    if ARG.dask:
        DICT={Time:-1,Y:ARG.Step,X:ARG.Step}
        A=A.chunk(DICT)
    #print(A)
    
    if ARG.n is not None:
        from itertools import count, product
        n=ARG.n+ARG.noffset
        Step=ARG.Step
        N,M=A.shape[1:]
        print(ARG.n,n)
        
        G=product(range(0,N,Step),range(0,M,Step))
        for i in range(n+1):
                nn,mm=next(G)
        A=A[:,nn:min((nn+Step),N),mm:min((mm+Step),M)]
        print(nn,mm,Step, A.shape)
    

    BDGBM=HarmBayMod(A, reformat_dict={"nan":-32768,"maxvalues":10000}, bandNameFormatter=lambda y: pd.Series(y))
    #Get all options on model 
    BDGBM.__dict__.update(ARG.__dict__)
    if ARG.dask:
        BDGBM.DaskRunstackBayesian()
    else:
        BDGBM.RunStackBayesian()
    BDGBM.FormatBayes()
    if (ARG.n is not None)&(not ARG.small):
        if ARG.netcdf:
            storenetcdf(BDGBM,ARG)
        else:
            storefm(BDGBM,ARG)
    else:
        BDGBM.Store(ARG.suffix)
    if ARG.graph:
        from matplotlib import pyplot as plt
        BDGBM.Data["stdinter_mean"].plot()
        plt.savefig("stdinter_mean.png")


