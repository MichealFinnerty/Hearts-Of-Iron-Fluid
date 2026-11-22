#!/usr/bin/env python3
"""Builds map-data.js (the Magma Europe province map) from Natural Earth public-domain data.
Usage: python3 build_map.py [density]   (density 1.0 = default, smaller = more provinces)
Needs: shapely>=2.1, numpy, pillow (preview only). Input files are downloaded to /tmp first (see README line in index.html)."""
import json, math, os, random, sys, time, urllib.request
import numpy as np, shapely
from shapely.geometry import shape, Point, Polygon, MultiPolygon, MultiPoint, box
from shapely.ops import unary_union, voronoi_diagram
from shapely.prepared import prep
from shapely import STRtree
random.seed(11); np.random.seed(11)
T0=time.time()
def log(*a): print(round(time.time()-T0),*a,flush=True)
EXT=(-11.5,34.5,60.0,71.0); BOX=box(*EXT); Q=0.005
K=float(sys.argv[1]) if len(sys.argv)>1 else 1.0
RAW='https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/'
def feats(n):
    p='/tmp/%s.geojson'%n
    if not os.path.exists(p): urllib.request.urlretrieve(RAW+n+'.geojson',p)
    return json.load(open(p,encoding='utf-8'))['features']
def polys(g):
    if g is None or g.is_empty: return []
    if g.geom_type=='Polygon': return [g]
    if g.geom_type in('MultiPolygon','GeometryCollection'): return [p for x in g.geoms for p in polys(x)]
    return []
def mp(g):
    ps=polys(g); return MultiPolygon(ps) if ps else MultiPolygon()
def valid(g): return g if g.is_valid else shapely.make_valid(g)
GROUPS=[('GER','German Reich','#4a5058','Germany|Austria|Czechia|Czech Republic|Poland|Netherlands|Belgium|Luxembourg|Denmark|Norway|France|Greece|Republic of Serbia|Serbia|Slovenia'),
('SOV','Soviet Union','#7b1531','Russia|Ukraine|Belarus|Lithuania|Latvia|Estonia|Moldova|Georgia|Armenia|Azerbaijan|Kazakhstan'),
('FIN','Finland','#a9bdd1','Finland'),('ROM','Romania','#c9ba62','Romania'),('HUN','Hungary','#e8a4a6','Hungary'),('SVK','Slovakia','#8da3c4','Slovakia'),
('ITA','Italy','#3f8d5a','Italy|Albania|Montenegro|Kosovo|San Marino|Vatican|Malta'),('BUL','Bulgaria','#6aa04c','Bulgaria|North Macedonia|Macedonia'),
('CRO','Croatia','#b8603f','Croatia|Bosnia and Herzegovina'),('VIC','Vichy France','#4d77b3',''),('ESP','Spain','#c7a54a','Spain|Andorra'),('POR','Portugal','#3f9296','Portugal'),
('SWI','Switzerland','#c35f5f','Switzerland|Liechtenstein'),('SWE','Sweden','#6f9fd0','Sweden'),('TUR','Turkey','#a3523d','Turkey'),('GBR','United Kingdom','#c95a82','United Kingdom'),('IRL','Ireland','#3da16b','Ireland')]
GK=[g[0] for g in GROUPS]
name2g={n:g[0] for g in GROUPS for n in g[3].split('|') if n}
# ---- countries (1941 owners are built from modern borders plus a few overrides)
items=[]
for f in feats('ne_10m_admin_0_countries'):
    g=mp(valid(shape(f['geometry'])).intersection(BOX))
    if not g.is_empty: items.append((f['properties']['ADMIN'],g))
geoms=[g for _,g in items]
try:
    simp=list(shapely.coverage_simplify(geoms,0.012)); log('coverage_simplify ok')
except Exception as e:
    log('coverage_simplify failed',e); simp=[g.simplify(0.006) for g in geoms]
pieces=[(name2g.get(n,'OTH'),n,mp(valid(g))) for (n,_),g in zip(items,simp)]
def ov(poly,to,src):
    p=Polygon(poly); new=[]
    for i,(k,n,g) in enumerate(pieces):
        if n in src:
            take=mp(g.intersection(p))
            if not take.is_empty: new.append((to,n,take)); pieces[i]=(k,n,mp(g.difference(p)))
    pieces.extend(new)
ov([(19.2,53.7),(22.88,53.7),(22.88,55.4),(19.2,55.4)],'GER',['Russia'])
ov([(20.9,55.05),(22.05,55.05),(22.05,55.9),(20.9,55.9)],'GER',['Lithuania'])
ov([(28.9,69.0),(31.3,69.0),(31.3,70.4),(28.9,70.4)],'FIN',['Russia'])
ov([(22.0,47.8),(24.75,47.8),(24.75,49.0),(22.0,49.0)],'HUN',['Ukraine'])
ov([(21.5,46.3),(24.0,46.1),(25.6,46.5),(26.5,47.0),(26.5,48.3),(21.5,48.3)],'HUN',['Romania'])
ov([(22.55,54.7),(22.55,53.95),(22.7,53.5),(23.0,53.1),(23.4,52.6),(24.4,52.6),(24.4,54.7)],'SOV',['Poland'])
ov([(6.3,46.15),(4.85,46.8),(3.3,46.6),(2.1,47.2),(0.6,46.9),(0.2,46.0),(-0.4,44.8),(-1.5,43.3),(-1.5,40.0),(10.0,40.0),(10.0,46.15)],'VIC',['France'])
G={}
for k,n,g in pieces: G.setdefault(k,[]).append(g)
for k in G: G[k]=mp(unary_union(G[k]))
LAND=mp(unary_union([G[k] for k in G]))
log('groups',{k:round(G[k].area) for k in G})
# ---- province seeds: dart throwing with a density that thins out towards the east and north
def rf(x,y): return K*min(.8,.30+.011*max(0,x-12)+.02*max(0,y-58))
def seeds_for(pg):
    pts=[];grid={};cs=.85
    for part in sorted([p for p in polys(pg) if p.area>=.04],key=lambda p:-p.area):
        pp=prep(part);minx,miny,maxx,maxy=part.bounds;fails=0;tries=0;mine=0
        while fails<2500 and tries<400000:
            tries+=1;x=random.uniform(minx,maxx);y=random.uniform(miny,maxy)
            if not pp.contains(Point(x,y)):continue
            r=rf(x,y);gx,gy=int(x//cs),int(y//cs);bad=False
            for i in range(gx-2,gx+3):
                for j in range(gy-2,gy+3):
                    for (px,py,pr) in grid.get((i,j),()):
                        d=(r+pr)/2
                        if (px-x)**2+(py-y)**2<d*d:bad=True;break
                    if bad:break
                if bad:break
            if bad:fails+=1;continue
            fails=0;mine+=1;pts.append((x,y));grid.setdefault((gx,gy),[]).append((x,y,r))
        if mine==0:
            rp=part.representative_point();pts.append((rp.x,rp.y));grid.setdefault((int(rp.x//cs),int(rp.y//cs)),[]).append((rp.x,rp.y,rf(rp.x,rp.y)))
    return pts
def cells(pg,pts):
    vd=voronoi_diagram(MultiPoint(pts),envelope=box(*pg.bounds).buffer(2));tree=STRtree([Point(p) for p in pts]);out={}
    for c in vd.geoms:
        idx=tree.query(c,predicate='contains')
        if len(idx)==0:continue
        try: piece=mp(c.intersection(valid(shapely.clip_by_rect(pg,*c.bounds))))
        except Exception: continue
        if not piece.is_empty: out[int(idx[0])]=piece
    return out
PR=[]
for gi,k in enumerate(GK):
    pg=G.get(k)
    if pg is None or pg.is_empty:continue
    pts=seeds_for(pg);c1=cells(pg,pts);new=[]
    for i,p in enumerate(pts):
        pc=c1.get(i)
        if pc is None:new.append(p);continue
        cen=pc.centroid;new.append((cen.x,cen.y) if pc.contains(cen) else (pts[i] if pc.contains(Point(*pts[i])) else (lambda r:(r.x,r.y))(pc.representative_point())))
    c2=cells(pg,new)
    n0=len(PR)
    for i,pc in c2.items():
        ps=[p for p in polys(pc) if p.area>=.0015]
        if not ps:continue
        pc=MultiPolygon(ps);rp=pc.representative_point();PR.append({'g':gi,'poly':pc,'c':(rp.x,rp.y)})
    log(k,len(PR)-n0,'provinces')
log('total provinces',len(PR))
# ---- terrain
reg=feats('ne_10m_geography_regions_polys')
mnt=unary_union([valid(shape(f['geometry'])) for f in reg if f['properties'].get('FEATURECLA')=='Range/mtn']).intersection(BOX)
plat=unary_union([valid(shape(f['geometry'])) for f in reg if f['properties'].get('FEATURECLA')=='Plateau']).intersection(BOX)
marsh=unary_union([box(26,51.2,31.8,52.9),box(3.5,51.0,7.3,53.6),box(28.5,44.3,30.0,45.5),box(21,53.3,24,54.2)])
mp_,pl_,ma_=prep(mnt),prep(plat),prep(marsh)
cities=[]
for f in feats('ne_10m_populated_places_simple'):
    p=f['properties'];x,y=f['geometry']['coordinates'][:2]
    if EXT[0]<x<EXT[2] and EXT[1]<y<EXT[3] and (p['pop_max']>=500000 or p['featurecla']=='Admin-0 capital'):cities.append((p['name'],x,y,p['pop_max'],p['featurecla']=='Admin-0 capital'))
ctree=STRtree([Point(c[1],c[2]) for c in cities])
def h01(i,s=0):return ((i*2654435761+s*40503)%1000003)/1000003
for i,pr in enumerate(PR):
    poly=pr['poly'];x,y=pr['c'];a=poly.area;t=0
    if poly.intersection(mnt).area>.5*a:t=3
    elif poly.intersection(plat).area>.55*a:t=2
    elif poly.intersection(marsh).area>.5*a:t=4
    else:
        pf=.85 if (y>=56 and x>=10) else .4 if y>=56 else .55 if (50<=y<56 and x>=22) else .3 if (y>=46 and x>=18) else .14
        if y>=65:pf=.5
        t=1 if h01(i)<pf else 0
    idx=ctree.query(poly,predicate='contains')
    if len(idx) and max(cities[j][3] for j in idx)>=1000000:t=5
    pr['t']=t
# ---- states: k-means of province centres inside every country
def kmeans(X,k,it=25):
    n=len(X)
    if k>=n:return list(range(n))
    C=[X[np.random.randint(n)]]
    for _ in range(k-1):
        d=np.min([((X-c)**2).sum(1) for c in C],axis=0);C.append(X[np.random.choice(n,p=d/d.sum())])
    C=np.array(C)
    for _ in range(it):
        lab=np.argmin(((X[:,None,:]-C[None,:,:])**2).sum(2),axis=1)
        for j in range(k):
            m=lab==j
            if m.any():C[j]=X[m].mean(0)
    return lab
ST=[]
for gi in range(len(GK)):
    ids=[i for i,p in enumerate(PR) if p['g']==gi]
    if not ids:continue
    X=np.array([PR[i]['c'] for i in ids]);lab=kmeans(X,max(1,round(len(ids)/6)))
    base=len(ST);m={}
    for i,l in zip(ids,lab):m.setdefault(int(l),[]).append(i)
    for l,mem in m.items():
        sid=len(ST);ST.append({'g':gi,'mem':mem})
        for i in mem:PR[i]['s']=sid
for s in ST:
    cx=np.mean([PR[i]['c'][0] for i in s['mem']]);cy=np.mean([PR[i]['c'][1] for i in s['mem']]);s['c']=(cx,cy)
    best=None
    for c in cities:
        d=(c[1]-cx)**2+(c[2]-cy)**2
        if d<1.2 and (best is None or c[3]>best[3]):best=c
    s['n']=best[0] if best else None
for i,s in enumerate(ST):
    if s['n'] is None:s['n']='%s %d'%(GROUPS[s['g']][1],i)
log('states',len(ST))
# ---- adjacency (provinces that share a border, also across country borders)
pt=STRtree([p['poly'] for p in PR]);adj=[]
for i,p in enumerate(PR):
    for j in pt.query(p['poly'].buffer(.012)):
        j=int(j)
        if j<=i:continue
        q=PR[j]['poly']
        if p['poly'].distance(q)<.006 and p['poly'].buffer(.012).intersection(q).area>4e-5:adj+=[i,j]
log('adjacency pairs',len(adj)//2)
# ---- outlines and decoration layers
def enc(coords,close=True):
    pts=[];last=None
    for x,y in coords:
        qx=round((x-EXT[0])/Q);qy=round((EXT[3]-y)/Q)
        if (qx,qy)!=last:pts.append((qx,qy));last=(qx,qy)
    if close and len(pts)>1 and pts[0]==pts[-1]:pts.pop()
    if len(pts)<(3 if close else 2):return None
    out=[pts[0][0],pts[0][1]]
    for a,b in zip(pts,pts[1:]):out+=[b[0]-a[0],b[1]-a[1]]
    return out
def rings(g,tol=0):
    out=[]
    for p in polys(g if tol==0 else g.simplify(tol)):
        for r in [p.exterior]+list(p.interiors):
            e=enc(r.coords)
            if e:out.append(e)
    return out
CB=LAND.boundary.buffer(.007)
def lns(g,tol=.004,sub=False):
    out=[];b=(g.simplify(tol) if tol else g).boundary
    if sub:b=b.difference(CB)
    for l in ([b] if b.geom_type=='LineString' else list(getattr(b,'geoms',[]))):
        if l.geom_type=='LineString' and len(l.coords)>=2:
            e=enc(list(l.coords),False)
            if e:out.append(e)
    return out
stateLines=[];
for s in ST:
    u=mp(unary_union([PR[i]['poly'] for i in s['mem']]));stateLines+=lns(u,.004,True)
countryLines=[]
for gi in range(len(GK)):
    pg=G.get(GK[gi])
    if pg is not None and not pg.is_empty:countryLines+=lns(pg,.004,True)
coast=lns(LAND,.006)
landR=rings(LAND,.006)
lakes=[]
for f in feats('ne_10m_lakes'):
    g=valid(shape(f['geometry'])).intersection(BOX)
    for p in polys(g):
        if p.area>.12:lakes+=rings(p,.01)
rivers=[]
for f in feats('ne_10m_rivers_lake_centerlines'):
    if (f['properties'].get('scalerank') or 99)>6:continue
    g=shape(f['geometry']).intersection(BOX)
    for l in ([g] if g.geom_type=='LineString' else list(getattr(g,'geoms',[]))):
        if l.geom_type=='LineString':
            e=enc(list(l.simplify(.008).coords),False)
            if e:rivers.append(e)
# country labels
labs={}
for gi,k in enumerate(GK):
    pg=G.get(k)
    if pg is None or pg.is_empty:continue
    big=max(polys(pg),key=lambda p:p.area);rp=big.representative_point();cen=big.centroid
    lx,ly=(cen.x,cen.y) if big.contains(cen) else (rp.x,rp.y)
    co=np.array(big.exterior.coords);co=co-co.mean(0);vals,vecs=np.linalg.eigh(np.cov(co.T));v=vecs[:,1]
    ang=-math.degrees(math.atan2(v[1],v[0]));ang=(ang+90)%180-90
    if abs(ang)>40:ang=0
    ang=max(-22,min(22,ang))
    b=big.bounds;nm=GROUPS[gi][1].upper()
    size=min((b[2]-b[0])*.8/(len(nm)*.78),(b[3]-b[1])*.5,5.5 if k=='SOV' else 3.2);size=max(.28,round(size,2))
    labs[gi]=[round(lx,2),round(ly,2),round(ang,1),size]
out={'ext':list(EXT),'q':Q,
 'groups':[{'k':g[0],'n':g[1],'c':g[2],'lab':labs.get(i)} for i,g in enumerate(GROUPS)],
 'prov':[[p['g'],p['t'],p['s'],round(p['c'][0],2),round(p['c'][1],2),[e for r in rings(p['poly']) for e in [r]]] for p in PR],
 'states':[[s['g'],s['n'],round(s['c'][0],2),round(s['c'][1],2)] for s in ST],
 'stateLines':stateLines,'countryLines':countryLines,'coast':coast,'land':landR,'lakes':lakes,'rivers':rivers,
 'cities':[[c[0],round(c[1],2),round(c[2],2),(2 if c[3]>=3000000 else 1 if c[3]>=1000000 else 0),1 if c[4] else 0] for c in cities],
 'adj':adj}
js='window.MAP='+json.dumps(out,separators=(',',':'),ensure_ascii=False)+';\n'
open('/mnt/user-data/outputs/map-data.js','w',encoding='utf-8').write(js)
log('wrote map-data.js',round(len(js)/1e6,2),'MB; provinces',len(PR),'states',len(ST),'cities',len(cities))
# ---- quick preview image
try:
    from PIL import Image, ImageDraw
    S=18;W=int((EXT[2]-EXT[0])*S);H=int((EXT[3]-EXT[1])*S);im=Image.new('RGB',(W,H),(43,74,111));d=ImageDraw.Draw(im)
    tx=lambda x:(x-EXT[0])*S;ty=lambda y:(EXT[3]-y)*S
    for p in polys(LAND):d.polygon([(tx(x),ty(y)) for x,y in p.exterior.coords],fill=(120,118,104))
    TT=[1.0,.86,.92,.78,.9,1.1]
    for p in PR:
        c=GROUPS[p['g']][2];rgb=tuple(int(c[i:i+2],16) for i in (1,3,5));f=TT[p['t']]*(.94+.1*h01(PR.index(p)))
        col=tuple(min(255,int(v*f)) for v in rgb)
        for q in polys(p['poly']):d.polygon([(tx(x),ty(y)) for x,y in q.exterior.coords],fill=col,outline=tuple(int(v*.7) for v in col))
    for l in countryLines:
        pass
    im.save('/tmp/preview.png');log('preview saved',im.size)
except Exception as e:log('preview failed',e)
