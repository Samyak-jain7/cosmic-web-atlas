"""Rebuild the atlas data files from the raw catalogues.

Reads   tools/raw/   (filled by tools/fetch_catalogues.sh)
Writes  near.txt, gal.txt  (base64 packed binary, loaded by index.html)
        and refreshes the <script id="meta"> block inside index.html.

Usage:  python3 tools/build_data.py
"""
import numpy as np, pandas as pd, json, warnings, os, re, base64
from astropy.io import fits
warnings.filterwarnings('ignore')
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW=os.path.join(ROOT,'tools','raw')+os.sep
OUT=ROOT+os.sep
rng=np.random.default_rng(7)
MPC=3.2616; D2R=np.pi/180
def scene(ra,dec,d):
    a=np.asarray(ra)*D2R; e=np.asarray(dec)*D2R; d=np.asarray(d)
    x=d*np.cos(e)*np.cos(a); y=d*np.cos(e)*np.sin(a); z=d*np.sin(e)
    return np.stack([x,z,-y],axis=-1)
# comoving distance table (Planck-like: H0 67.8, Om 0.308), Mly
zz=np.linspace(0,0.7,7001); E=1/np.sqrt(0.308*(1+zz)**3+0.692)
DC=np.concatenate([[0],np.cumsum((E[1:]+E[:-1])/2*np.diff(zz))])*299792.458/67.8*MPC
def Dc(z): return np.interp(z,zz,DC)
# heliocentric -> CMB frame (dipole apex l=264.02, b=48.25, 369.8 km/s)
GM=np.array([[-0.0548755604,-0.8734370902,-0.4838350155],[0.4941094279,-0.4448296300,0.7469822445],[-0.8676661490,-0.1980763734,0.4559837762]])
l,b=264.02*D2R,48.25*D2R; apex=GM.T@np.array([np.cos(b)*np.cos(l),np.cos(b)*np.sin(l),np.sin(b)])
def zcmb(ra,dec,z):
    a=np.asarray(ra)*D2R; e=np.asarray(dec)*D2R
    u=np.stack([np.cos(e)*np.cos(a),np.cos(e)*np.sin(a),np.sin(e)],axis=-1)
    return z+369.8/299792.458*(u@apex)
meta={}
# ---------- stars ----------
h=pd.read_csv(RAW+'hyg.csv',low_memory=False)
s=h[(h.dist>0)&(h.dist<100000)].copy()
chk=scene(s.ra*15,s.dec,s.dist); own=np.stack([s.x,s.z,-s.y],axis=-1)
print('HYG xyz convention max diff (pc):',np.abs(chk-own).max())
P=(own*3.2616e-6).astype('<f4')
ci=np.clip((s.ci.fillna(0.65).values+0.4)/2.4*255,0,255).astype('u1')
am=np.clip((s.absmag.values+8)/24*255,0,255).astype('u1')
idx=rng.permutation(len(s)); P,ci,am=P[idx],ci[idx],am[idx]
meta['stars']=dict(n=int(len(s)))
star_blob=P.tobytes()+ci.tobytes()+am.tobytes()
want=['Proxima Centauri','Rigil Kentaurus',"Barnard's Star",'Wolf 359','Lalande 21185','Sirius','Ran','Procyon','Altair','Vega','Fomalhaut','Pollux','Arcturus','Capella','Castor','Aldebaran','Regulus','Achernar','Canopus','Spica','Acrux','Hadar','Betelgeuse','Antares','Rigel','Deneb','Polaris','Alcyone','Bellatrix','Alnilam','Mizar','Algol']
nice={'Rigil Kentaurus':'Alpha Centauri','Ran':'Epsilon Eridani','Alcyone':'Alcyone (Pleiades)'}
pins=[]
for nm in want:
    r=s[s.proper==nm]
    if len(r)==0: print('missing star',nm); continue
    r=r.iloc[0]; p=np.array([r.x,r.z,-r.y])*3.2616e-6
    pins.append(dict(n=nice.get(nm,nm),p=[float('%.4g'%v) for v in p],ly=round(float(r.dist)*3.2616,2 if r.dist<10 else 0),mag=float(r.mag),sp=str(r.spect) if isinstance(r.spect,str) else ''))
meta['starPins']=pins
# ---------- Local Group / Local Volume ----------
def mly(dm): return 10**((dm+5)/5)*3.2616e-6
d=pd.read_csv(RAW+'lv_dwarf_all.csv'); d=d[d.distance_modulus.notna()&d.ra.notna()].copy(); d['D']=mly(d.distance_modulus)
g=pd.concat([pd.read_csv(RAW+'lv_gc_harris.csv'),pd.read_csv(RAW+'lv_gc_mw_new.csv')]); g=g[g.distance_modulus.notna()].copy(); g['D']=mly(g.distance_modulus)
dp=scene(d.ra,d.dec,d.D); gp=scene(g.ra,g.dec,g.D)
mv=d.M_V.fillna(-8).clip(-19,0)
meta['dwarfs']=dict(n=int(len(d)),p=[float('%.5g'%v) for v in dp.ravel()],mv=[round(float(v),1) for v in mv],max=float(d.D.max()),nmw=int((d.host=='mw').sum()),nm31=int((d.host=='m_031').sum()))
meta['gcs']=dict(n=int(len(g)),p=[float('%.4g'%v) for v in gp.ravel()])
def lvpin(key,label,kind='dwarf'):
    src=d if kind=='dwarf' else g
    r=src[src.key==key]
    if len(r)==0: print('missing',key); return
    r=r.iloc[0]; p=scene(r.ra,r.dec,r.D)
    return dict(n=label,p=[float('%.5g'%v) for v in p],mly=float('%.3g'%r.D),ra=round(float(r.ra),2),dec=round(float(r.dec),2))
lgp=[]
for k,lab in [('lmc','Large Magellanic Cloud'),('smc','Small Magellanic Cloud'),('sagittarius_1','Sagittarius dwarf'),('fornax_1','Fornax dwarf'),('sculptor_1','Sculptor dwarf'),('leo_1','Leo I'),('draco_1','Draco dwarf'),('m_032','M32'),('ngc_0205','M110'),('ic_0010','IC 10'),('ngc_6822','NGC 6822'),('ic_1613','IC 1613'),('wlm','WLM'),('ngc_0055','NGC 55'),('ngc_0300','NGC 300'),('ngc_5194','M51 Whirlpool')]:
    q=lvpin(k,lab)
    if q: lgp.append(q)
meta['lgPins']=lgp
gcp=[]
for k,lab in [('ngc_5139','Omega Centauri'),('ngc_0104','47 Tucanae'),('ngc_6205','M13'),('ngc_6121','M4'),('ngc_7078','M15'),('ngc_2419','NGC 2419')]:
    q=lvpin(k,lab,'gc')
    if q: gcp.append(q)
meta['gcPins']=gcp
for k in ['lmc','smc']: print(k, d[d.key==k][['ra','dec','D']].values)
# ---------- OpenNGC ----------
n=pd.read_csv(RAW+'NGC.csv',sep=';',low_memory=False)
gal=n[n.Type.isin(['G','GPair','GTrpl','GGroup'])&n.Redshift.notna()].copy()
def hms(sx):
    a=[float(v) for v in sx.split(':')]; return (a[0]+a[1]/60+a[2]/3600)*15
def dms(sx):
    sg=-1 if sx.strip().startswith('-') else 1; a=[abs(float(v)) for v in sx.split(':')]; return sg*(a[0]+a[1]/60+a[2]/3600)
gal['ra']=gal.RA.map(hms); gal['dec']=gal.Dec.map(dms)
# hosts in the Local Volume DB: place at the median measured distance of their satellites
hostd=d.groupby('host').D.median()
u=np.stack([np.cos(gal.dec*D2R)*np.cos(gal.ra*D2R),np.cos(gal.dec*D2R)*np.sin(gal.ra*D2R),np.sin(gal.dec*D2R)],axis=-1)
m87=np.array([np.cos(12.39*D2R)*np.cos(187.706*D2R),np.cos(12.39*D2R)*np.sin(187.706*D2R),np.sin(12.39*D2R)])
sep=np.degrees(np.arccos(np.clip(u@m87,-1,1)))
zc=zcmb(gal.ra,gal.dec,gal.Redshift.values)
virgo=(sep<6)&(gal.Redshift.values*299792.458<2800)
D=Dc(np.clip(zc,0,0.69)); D[virgo]=55.62+rng.normal(0,3.0,virgo.sum())
keep=virgo|(zc>0.002)
gal=gal[keep]; D=D[keep]; print('NGC kept',len(gal),'virgo members',int(virgo.sum()))
np_=scene(gal.ra,gal.dec,D).astype('<f4')
early=gal.Hubble.fillna('').str.match(r'^(E|S0)').values.astype('u1')
meta['ngc']=dict(n=int(len(gal)),virgo=int(virgo.sum()))
ngc_blob=np_.tobytes()+early.tobytes()
gp2=[]
def ngcpin(name,label,dist=None):
    r=n[n.Name==name]
    if len(r)==0: print('missing',name); return
    r=r.iloc[0]; ra=hms(r.RA); dec=dms(r.Dec)
    how='z'
    if dist is None:
        zc1=float(zcmb(ra,dec,float(r.Redshift))); dist=float(Dc(zc1))
    else: how='sat'
    p=scene(ra,dec,dist)
    gp2.append(dict(n=label,p=[float('%.5g'%v) for v in p],mly=float('%.3g'%dist),ra=round(ra,2),dec=round(dec,2),how=how))
print('host medians Mly',{k:round(float(hostd[k]),2) for k in ['m_031','m_033','m_081','ngc_5128','ngc_0253','m_083','ngc_4594','m_101','ngc_3379'] if k in hostd.index})
for nm,lab,hk in [('NGC0598','Triangulum Galaxy','m_033'),('NGC3031','M81','m_081'),('NGC5128','Centaurus A','ngc_5128'),('NGC0253','Sculptor Galaxy','ngc_0253'),('NGC5236','M83','m_083'),('NGC4594','Sombrero Galaxy','ngc_4594'),('NGC5457','M101','m_101'),('NGC3379','M105 (Leo I group)','ngc_3379')]:
    if hk in hostd.index: ngcpin(nm,lab,float(hostd[hk]))
    else: print('no host',hk)
for nm,lab in [('NGC4486','M87'),('NGC1316','Fornax A'),('NGC4889','NGC 4889 (Coma)'),('NGC1275','Perseus A'),('NGC6166','NGC 6166'),('NGC3311','NGC 3311 (Hydra)'),('NGC4696','NGC 4696 (Centaurus)')]:
    ngcpin(nm,lab)
gp2=[q for q in gp2]
# M87 is a Virgo member: pin it at the cluster distance
for q in gp2:
    if q['n']=='M87': q['p']=[float('%.5g'%v) for v in scene(187.706,12.391,55.62)]; q['mly']=55.6; q['how']='virgo'
meta['galPins']=gp2
# ---------- NSA (z<0.05) + SDSS DR8 (z>=0.05) ----------
a=np.load(RAW+'nsa.npy')
m=(a['ZDIST']>0.002)&(a['ZDIST']<0.05)
nra,ndec,nz,nd4=a['RA'][m],a['DEC'][m],a['ZDIST'][m],a['D4000'][m]
parts=[]
for f in ['sdss1.fit.gz','sdss2.fit.gz']:
    with fits.open(RAW+f) as hd:
        t=hd[1].data
        parts.append((np.array(t['ra'],'f8'),np.array(t['dec'],'f8'),np.array(t['z'],'f8'),np.array(t['zErr'],'f8'),np.array(t['d4000'],'f8')))
sra,sdec,sz,serr,sd4=[np.concatenate(x) for x in zip(*parts)]
ms=(sz>=0.05)&(sz<0.5)&(serr>0)&(serr<0.002)
sra,sdec,sz,sd4=sra[ms],sdec[ms],sz[ms],sd4[ms]
szc=zcmb(sra,sdec,sz)
ra=np.concatenate([nra,sra]); dec=np.concatenate([ndec,sdec]); zc_=np.concatenate([nz,szc]); d4=np.concatenate([nd4,sd4])
print('NSA',len(nra),'SDSS',len(sra),'total',len(ra),'z max',zc_.max())
Dg=Dc(np.clip(zc_,0,0.69)); pos=scene(ra,dec,Dg)
SC=6500/32767.0
q=np.clip(np.round(pos/SC),-32767,32767).astype('<i2')
c=np.clip((np.nan_to_num(d4,nan=1.5)-1.1)/(2.1-1.1)*255,0,255).astype('u1')
idx=rng.permutation(len(q)); q,c=q[idx],c[idx]
meta['gal']=dict(n=int(len(q)),scale=SC,nsa=int(len(nra)),sdss=int(len(sra)),zmax=0.5,dmax=float(Dg.max()))
# sanity: overdensities at sourced positions
def count(ra0,dec0,d0,rad):
    c0=scene(ra0,dec0,d0); return int((np.linalg.norm(pos-c0,axis=1)<rad).sum())
vol=lambda r:4/3*np.pi*r**3
for nm,(r0,d0,D0) in {'Coma':(194.95,27.98,336),'Hercules A2151':(241.31,17.75,528),'Bootes void':(222.5,46,700),'random field':(170,40,500)}.items():
    for rad in (25,60): print(nm,rad,'Mly ->',count(r0,d0,D0,rad),'gal; density',round(count(r0,d0,D0,rad)/vol(rad)*1e6,1),'per 1e6 Mly^3')
sg=(dec>-4)&(dec<6)&(ra>130)&(ra<206)
hist,edges=np.histogram(Dg[sg],bins=np.arange(400,1700,100)); print('equatorial strip RA130-206 counts by distance (Mly):',dict(zip(edges[:-1].astype(int).tolist(),hist.tolist())))
near=star_blob+ngc_blob; galb=q.tobytes()+c.tobytes()
meta['off']=dict(starP=0,starC=len(P.tobytes()),starM=len(P.tobytes())+len(ci),ngcP=len(star_blob),ngcE=len(star_blob)+len(np_.tobytes()))
open(OUT+'near.txt','w').write(base64.b64encode(near).decode())
open(OUT+'gal.txt','w').write(base64.b64encode(galb).decode())
html=open(OUT+'index.html',encoding='utf-8').read()
pat=re.compile(r'(<script id="meta" type="application/json">)(.*?)(</script>)',re.S)
assert pat.search(html),'meta block not found in index.html'
html=pat.sub(lambda m:m.group(1)+json.dumps(meta,separators=(',',':'))+m.group(3),html,count=1)
open(OUT+'index.html','w',encoding='utf-8').write(html)
for f in ['near.txt','gal.txt','index.html']: print(f,os.path.getsize(OUT+f),'bytes')
