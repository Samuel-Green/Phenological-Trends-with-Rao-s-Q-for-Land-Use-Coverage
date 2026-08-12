import subprocess
import os
import json
import pandas as pd
import xarray as xr
import numpy as np
import geopandas as gpd

#print(os.listdir("./DataInput"))
print(os.listdir("./"))
#subprocess.run(["unzip","./DataInput/"]+os.listdir("./DataInput")[0])
subprocess.run(["unzip","DataInput.zip", "-d","./DataInput"])
print(os.listdir("./"))

try:
    inputfile=[ x for x in os.listdir("./DataInput") if x.split(".")[-1].lower() in ["tif", "tiff", "nc","envi"]][0]
except IndexError:
    print("Input file should have extension name tif, nc or envi")

ARG=json.load(open("vlabparams.json","r"))
if "mask.zip" in os.listdir("./"):
    subprocess.run(["unzip","mask.zip","-d","./MASK"])
    try:
        L=os.listdir("./MASK")
    except FileNotFoundError:
        path=None
        ARG["mask"]=False
    else:
        if len(L)==1:
            path="./MASK/"+L[0]
        else:
            path="./MASK/"+[ x for x in L if x.split(".")[-1].lower() in ["envi","shp"]][0]
    if path:
        try:
            V=gpd.read_file(path)
        except Exception as inst:
            if str(type(inst))!="<class 'fiona.errors.DriverError'>":
                raise(KeyError)
            ARG["mask"]="false"
        else:
            ARG["mask"]=path
print("ARG")
print(ARG)
arg=""
for k,v in ARG.items():
    if (v is False)|(v in ["False","false","F"]):
        continue
    else:
        if (v is True)|(v=="true"):
            v=""
        arg+=" --"+" ".join([k,str(v)])

print("arguments")
print(arg)
arg=""
for k,v in ARG.items():
    if (v is False)|(v in ["False","false","F"]):
        continue
    else:
        if (v is True)|(v=="true"):
            v=""
        arg+=" --"+" ".join([k,str(v)])

print("arguments")
print(arg)

process = subprocess.Popen(["python3", "PPresilBayes.py","--input","./DataInput/"+inputfile]+arg.split(), stdout=subprocess.PIPE)
while True:
    output = process.stdout.readline()
    if output == b'' and process.poll() is not None:
        break
    if output:
        print( output.strip())
rc = process.poll()
#subprocess.run(["python3", "PPresilBayes.py","--input","./DataInput/"+inputfile]+arg.split(),stdout=subprocess.PIPE).stdout


#WMS
import xarray as xr
A=xr.open_dataset(ARG["suffix"]+".nc")
A["sdinter_mean"]*=10000
A["sdinter_mean"]=A["sdinter_mean"].astype("int16")
A["sdinter_mean"].to_netcdf(ARG["suffix"]+"sdinter_mean.nc")
Min,Max=A.sdinter_mean.quantile(q=[0.05,0.975]).values
bins=10
Z=pd.read_csv("paletteViridis.csv", header=None)
Z[0]=np.linspace(Min,Max,bins).astype("int")
Z.to_csv("paletteViridisMod.csv", index=False, header=False)


CMD="gdaldem color-relief -alpha -nearest_color_entry NETCDF:{suffix}sdinter_mean.nc:sdinter_mean paletteViridisMod.csv {suffix}sdinter_mean.png".format(suffix=ARG["suffix"])
subprocess.run(CMD.split())

#CMD="rgb2pct.py {suffix}sdinter_mean.png {suffix}sdinter_mean_pct.png".format(suffix=ARG["suffix"])
#subprocess.run(CMD.split())
#gdalwarp -s_srs_
#subprocess.run(["gdal_translate", 'NETCDF:'+ARG["suffix"]+".nc:sdinter_mean", ARG["suffix"]+"sdinter_mean.tif"])




subprocess.run(["gdal2tiles.py", "-z","2-14", "-c", "CNR-IIA/vicario 2020", "-s",  A.pyproj_srs, ARG["suffix"]+"sdinter_mean.png"])


ID=os.environ['BENGINERUNID']
print(ID)
CMD='/root/.local/bin/aws s3 sync {suffix}sdinter_mean s3://testtilebucket/{ID}{V1}/3857 --acl public-read'.format(ID=ID, V1="V1",suffix=ARG["suffix"])
subprocess.run(CMD.split())

D={"url":"https://eddipf6en8.execute-api.us-east-1.amazonaws.com/alpha/wms?", "protocol":"urn:ogc:serviceType:WebMapService:1.1.1:HTTP", "name":ID+"V1"} 
#Putting the legend on the map
Z[0]/=10000
Ra=pd.DataFrame()
Ra["U"]=Z.iloc[:,0].diff()[1:]/2+Z.iloc[:-1,0].values
Ra["L"]=Z.iloc[:,0].diff()[:-1]/2+Z.iloc[:-1,0].values
#Ra.index=Ra.index.astype("str")
Ra.loc[0,["U","L"]]=[np.nan,Ra.iloc[0,0]]
Ra=Ra.sort_index()
Ra=Ra.join(Z)
def Form(x):
    temp={'type': 'square', 'color': 'rgba(72,33,114,255)', 'text': '0 - 59' }
    #print(x)
    temp["text"]='{:.2e} - {:.2e}'.format(x.values[0],x.values[1])
    temp['color']='rgba({0},{1},{2},{3})'.format(*x[-4:])
    return temp
Legend=list(Ra.apply(lambda x: Form(x), axis=1))
D["legendList"]=Legend
with open('wms.json', 'w') as outfile:
    json.dump(D, outfile)

print(D)

subprocess.run(["zip","Output.zip",ARG["suffix"]+".nc"])
