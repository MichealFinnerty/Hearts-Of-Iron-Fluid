#!/usr/bin/env python3
"""INVICTA EUROPA - builds map-data.js: Europe, North Africa, East Africa down to Kenya and the Near East, with two
start dates (1 Aug 1939 and 1 Jan 1941), from Natural Earth public-domain data.

How the map is made
  1. PROVINCES: first-level administrative units from Natural Earth admin-1, re-cut to match 1939 borders
     (Prussian provinces and German states, 1939 Polish voivodeships, French departements with the Seine/Seine-et-Oise
     merges, Italian provinces without the post-1939 splits, Italian Istria/Fiume/Zara/Dodecanese, Hungarian Felvidek,
     Romanian Transylvania/Dobruja, Soviet-Finnish Karelia, and so on).
  2. TILES: every province is divided into tiles (Voronoi cells clipped to the province polygon), so tile edges always lie
     exactly on province borders and a tile never belongs to two provinces.
  2b. MERGING: county-sized provinces are joined with neighbours of the same country and owner until they are about the size of
     Munster (~24,600 km2); big provinces (Bavaria, Russian oblasts) and small special cases (Malta, Gibraltar, Berlin) stay as they are.
     Ireland is Munster, Leinster and Connacht. See MUNSTER / KEEP_ALONE / IRL_PROV in the merge block. 730 provinces result.
  3. Ownership for 1 Aug 1939 and 1 Jan 1941 is stored per province.
  4. NAMES: only country names and cities are drawn on the map, never province names, and only historical names: a country's colonies
     and exclaves carry the name of the country that holds them (Libya reads ITALY, Algeria reads FRANCE ...), as in HOI4.
     Each name is fitted (see 'COUNTRY NAME LABELS' / fit_label) as the biggest text box that never covers anybody else's land
     but may stretch out over the sea. Egypt belongs to the UK and Italy holds all of East Africa (British Somaliland from Aug 1940).

Usage: python3 build_map.py [density] [--preview] [--provinces-only]
       (density 1.0 = default, bigger = fewer, larger tiles)
Needs: pip install shapely numpy pillow.  Natural Earth files are downloaded once into ./ne_data (or $NE_DIR)."""
import json, math, os, random, sys, time, urllib.request
import numpy as np, shapely
from shapely.geometry import shape, Point, Polygon, MultiPolygon, MultiPoint, box
from shapely.ops import unary_union, voronoi_diagram
from shapely.prepared import prep
from shapely import STRtree
random.seed(11); np.random.seed(11)
T0=time.time()
def log(*a): print(round(time.time()-T0),*a,flush=True)
HERE=os.path.dirname(os.path.abspath(__file__)); DATA=os.environ.get('NE_DIR',os.path.join(HERE,'ne_data'))
EXT=(-18.5,-5.0,60.0,71.0); BOX=box(*EXT); Q=0.005
args=[a for a in sys.argv[1:] if not a.startswith('--')]; K=float(args[0]) if args else 1.0
RAW='https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/'
def feats(n):
    os.makedirs(DATA,exist_ok=True);p=os.path.join(DATA,n+'.geojson')
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
def P(*pts): return Polygon(pts)
GROUPS=[('GER','German Reich','#4a5058'),('SOV','Soviet Union','#7b1531'),('POL','Poland','#c05050'),('FRA','France','#4d77b3'),('VIC','Vichy France','#7f93b8'),('FFR','Free France','#3a5fa8'),
('GBR','United Kingdom','#c95a82'),('ITA','Italy','#3f8d5a'),('ESP','Spain','#c7a54a'),('POR','Portugal','#3f9296'),('SWE','Sweden','#6f9fd0'),('FIN','Finland','#a9bdd1'),('EST','Estonia','#7aa6c2'),
('LAT','Latvia','#a05a6a'),('LIT','Lithuania','#8aa84f'),('DEN','Denmark','#d07a6a'),('NOR','Norway','#7a5a8a'),('NED','Netherlands','#e0903a'),('BEL','Belgium','#d2b44a'),('LUX','Luxembourg','#8fb3c9'),
('SWI','Switzerland','#c35f5f'),('YUG','Yugoslavia','#8a6f4a'),('GRE','Greece','#6fa0d8'),('ROM','Romania','#c9ba62'),('HUN','Hungary','#e8a4a6'),('SVK','Slovakia','#8da3c4'),('BUL','Bulgaria','#6aa04c'),
('TUR','Turkey','#a3523d'),('IRL','Ireland','#3da16b'),('IRQ','Iraq','#7c8f50'),('IRN','Iran','#a58f60'),('SAU','Saudi Arabia','#b3a06c'),('YEM','Yemen','#a07f5c'),('LBR','Liberia','#5c8c6c'),
('DAN','Free City of Danzig','#8fb08a'),('VAT','Vatican City','#e6d27a'),('MON','Monaco','#d65a5a'),('SMR','San Marino','#6fb3e0'),('AND','Andorra','#d0a030'),('LIE','Liechtenstein','#5a6fc4')]
GI={g[0]:i for i,g in enumerate(GROUPS)}
# ---------------------------------------------------------------- default owners (1 Aug 1939, 1 Jan 1941) by modern country
DEF={}
def own(names,a,b=None):
    for n in names.split('|'): DEF[n]=(a,b or a)
own('Germany|Austria|Czech Republic','GER'); own('Slovakia','SVK'); own('Netherlands|Belgium|Luxembourg|Denmark|Norway','GER')
DEF['Netherlands']=('NED','GER');DEF['Belgium']=('BEL','GER');DEF['Luxembourg']=('LUX','GER');DEF['Denmark']=('DEN','GER');DEF['Norway']=('NOR','GER');DEF['France']=('FRA','GER');DEF['Monaco']=('FRA','GER')
own('Estonia','EST','SOV');own('Latvia','LAT','SOV');own('Lithuania','LIT','SOV')
own('Russia|Ukraine|Belarus|Georgia|Armenia|Azerbaijan|Kazakhstan|Uzbekistan|Turkmenistan','SOV'); own('Moldova','ROM','SOV'); own('Romania','ROM'); own('Hungary','HUN'); own('Bulgaria','BUL')
own('Finland|Aland','FIN'); own('Sweden','SWE'); own('Switzerland|Liechtenstein','SWI')
own('Italy|Albania|Libya|Eritrea|Ethiopia|Somalia|San Marino|Vatican','ITA'); own('Spain|Andorra|Western Sahara|Equatorial Guinea','ESP'); own('Portugal|Guinea Bissau|Angola|Sao Tome and Principe','POR')
own('United Kingdom|Malta|Cyprus|Northern Cyprus|Akrotiri Sovereign Base Area|Dhekelia Sovereign Base Area|Cyprus No Mans Area|Israel|West Bank|Gaza Strip|Jordan|Kuwait|Qatar|United Arab Emirates|Oman|Bahrain|Sudan|S. Sudan|Kenya|Uganda|United Republic of Tanzania|Nigeria|Ghana|Sierra Leone|Gambia|Isle of Man|Jersey|Guernsey|Gibraltar','GBR')
DEF['Poland']=('POL','GER'); DEF['Somaliland']=('GBR','ITA')   # British Somaliland fell to Italy in Aug 1940
own('Ireland','IRL'); own('Republic of Serbia|Croatia|Bosnia and Herzegovina|Slovenia|Montenegro|Macedonia|Kosovo','YUG'); own('Greece','GRE'); own('Turkey','TUR')
own('Egypt','GBR'); own('Iraq','IRQ'); own('Iran','IRN'); own('Saudi Arabia','SAU'); own('Yemen','YEM'); own('Liberia','LBR'); own('Rwanda|Burundi|Democratic Republic of the Congo','BEL')
own('Morocco|Algeria|Tunisia|Mauritania|Mali|Niger|Senegal|Guinea|Ivory Coast|Burkina Faso|Benin|Togo|Djibouti|Syria|Lebanon','FRA','VIC'); own('Chad|Republic of the Congo|Gabon|Central African Republic|Cameroon','FRA','FFR')
# ---------------------------------------------------------------- load admin-1, clip to the map, keep the shared borders consistent
A1=[]
for f in feats('ne_10m_admin_1_states_provinces'):
    p=f['properties']; g=mp(valid(shape(f['geometry'])).intersection(BOX))
    if not g.is_empty: A1.append((p,g))
try: simp=list(shapely.coverage_simplify([g for _,g in A1],0.012)); log('coverage_simplify ok',len(A1),'admin-1 units')
except Exception as e: log('coverage_simplify failed',e); simp=[g.simplify(0.006) for _,g in A1]
A1=[(p,mp(valid(g))) for (p,_),g in zip(A1,simp)]
# ---------------------------------------------------------------- the 1939 province table
REGS={}   # province key -> geometry
OWNP={}   # province key -> (owner 1939, owner 1941) as group codes
KEYADM={}
CUR=['']
def put(key,g,o):
    g=mp(g)
    if g.is_empty: return
    KEYADM.setdefault(key,set()).add(CUR[0])
    REGS[key]=mp(unary_union([REGS[key],g])) if key in REGS else g
    if o is not None: OWNP[key]=o
def move(src,key,o=None):
    if src in REGS: g=REGS.pop(src); oo=OWNP.pop(src,None); put(key,g,o or oo)
def carve(poly,key,o=None,frm=None):
    """cut poly out of the provinces in frm (a list of keys; default every province) and give it to key"""
    ks=[k for k in (frm if frm is not None else list(REGS)) if k in REGS]
    for k in ks:
        if k==key: continue
        t=mp(REGS[k].intersection(poly))
        if t.is_empty or t.area<1e-7: continue
        REGS[k]=mp(REGS[k].difference(poly)); put(key,t,o or OWNP.get(k))
        if REGS[k].is_empty: REGS.pop(k); OWNP.pop(k,None)
POL1939=Polygon([(26.6,55.7),(27.0,55.75),(27.45,55.82),(28.3,55.65),(28.3,55.4),(27.95,55.05),(27.9,54.75),(27.45,54.4),(27.15,54.1),(27.1,53.7),(27.3,52.9),(27.6,52.2),(27.7,51.5),(27.4,50.7),(26.9,50.1),(26.4,49.6),(26.3,49.1),(26.6,48.55),(25.8,48.25),(24.9,48.15),(23.8,48.0),(22.9,48.5),(22.2,49.05),(21.0,49.4),(20.0,49.25),(19.0,49.4),(18.6,49.7),(18.3,49.95),(18.25,50.12),(18.45,50.18),(18.7,50.22),(18.82,50.25),(18.83,50.45),(18.6,50.75),(18.35,50.9),(18.25,51.05),(17.95,51.2),(17.3,51.5),(16.8,51.68),(16.4,51.95),(15.85,52.2),(15.8,52.5),(16.1,52.85),(16.7,52.95),(16.9,53.4),(17.45,53.6),(17.65,53.9),(17.75,54.2),(18.0,54.5),(18.0,54.85),(18.9,54.85),(18.95,54.45),(18.88,53.95),(18.88,53.7),(19.2,53.55),(19.9,53.55),(20.05,53.5),(20.3,53.3),(20.55,53.2),(20.9,53.3),(21.1,53.4),(21.5,53.42),(21.95,53.5),(22.3,53.55),(22.45,53.75),(22.65,53.95),(22.75,54.15),(22.85,54.4),(23.0,54.3),(23.3,54.22),(23.7,54.17),(24.0,54.17),(24.3,54.3),(24.55,54.42),(24.75,54.6),(24.85,54.75),(25.0,54.85),(25.3,54.95),(25.7,55.12),(26.0,55.25),(26.25,55.32),(26.5,55.5)])
IT_MERGE={'Pordenone':'Udine','Isernia':'Campobasso','Verbano-Cusio-Ossola':'Novara','Biella':'Vercelli','Lecco':'Como','Lodi':'Milano','Monza e Brianza':'Milano','Prato':'Firenze','Crotene':'Catanzaro','Vibo Valentia':'Catanzaro',
 'Barletta-Andria Trani':'Bari','Fermo':'Ascoli Piceno','Carbonia-Iglesias':'Cagliari','Medio Campidano':'Cagliari','Oristrano':'Cagliari','Ogliastra':'Nuoro','Olbia-Tempio':'Sassari','Forlì-Cesena':'Forlì','Rimini':'Forlì',
 'Bozen':'Bolzano','Aoste':'Aosta','Turin':'Torino','Latina':'Littoria','Aoste ':'Aosta'}
FR_MERGE={'Paris':'Seine','Hauts-de-Seine':'Seine','Seine-Saint-Denis':'Seine','Val-de-Marne':'Seine','Yvelines':'Seine-et-Oise','Essonne':'Seine-et-Oise',"Val-d'Oise":'Seine-et-Oise','Seien-et-Marne':'Seine-et-Marne','Corse-du-Sud':'Corsica','Haute-Corse':'Corsica'}
FIN_MERGE={'Uusimaa':'Uusimaa','Finland Proper':'Turku and Pori','Satakunta':'Turku and Pori','Tavastia Proper':'Häme','Päijät-Häme':'Häme','Pirkanmaa':'Häme','Kymenlaakso':'Viipuri','South Karelia':'Viipuri',
 'Southern Savonia':'Mikkeli','Northern Savonia':'Kuopio','North Karelia':'Kuopio','Central Finland':'Vaasa','Ostrobothnia':'Vaasa','Central Ostrobothnia':'Vaasa','Southern Ostrobothnia':'Vaasa','Northern Ostrobothnia':'Oulu','Kainuu':'Oulu','Lapland':'Lapland'}
HR_MERGE={'Vukovarsko-Srijemska':'Slavonia','Osjecko-Baranjska':'Slavonia','Brodsko-Posavska':'Slavonia','Viroviticko-Podravska':'Slavonia','Pozesko-Slavonska':'Slavonia',
 'Grad Zagreb':'Croatia','Zagrebacka':'Croatia','Karlovacka':'Croatia','Sisacko-Moslavacka':'Croatia','Bjelovarsko-bilogorska':'Croatia','Koprivničko-Križevačka':'Croatia','Varaždinska':'Croatia','Medimurska':'Croatia','Krapinsko-Zagorska':'Croatia',
 'Licko-Senjska':'Lika','Zadarska':'Dalmatia','Šibensko-Kninska':'Dalmatia','Splitsko-Dalmatinska':'Dalmatia','Dubrovacko-Neretvanska':'Dalmatia','Istarska':'Pola','Primorsko-Goranska':'Fiume'}
SRB={'Severno-Backi':'Bačka','Južno-Backi':'Bačka','Zapadno-Backi':'Bačka','Severno-Banatski':'Banat','Srednje-Banatski':'Banat','Južno-Banatski':'Banat','Sremski':'Syrmia','Grad Beograd':'Belgrade'}
GR_NAMES={'Anatoliki Makedonia kai Thraki':'Eastern Macedonia and Thrace','Kentriki Makedonia':'Central Macedonia','Dytiki Makedonia':'Western Macedonia','Ayion Oros':'Central Macedonia','Dytiki Ellada':'Western Greece','Ionioi Nisoi':'Ionian Islands',
 'Ipeiros':'Epirus','Kriti':'Crete','Notio Aigaio':'South Aegean','Voreio Aigaio':'North Aegean','Peloponnisos':'Peloponnese','Stereá Elláda':'Central Greece','Thessalia':'Thessaly','Attiki':'Attica'}
PL_NAMES={'Greater Poland':'Poznań','Kuyavian-Pomeranian':'Pomorze','Pomeranian':'Pomorze','Lesser Poland':'Kraków','Subcarpathian':'Rzeszów','Lublin':'Lublin','Masovian':'Warsaw','Łódź':'Łódź','Świętokrzyskie':'Kielce','Podlachian':'Białystok','Silesian':'Silesia'}
CZ={'Prague':'Bohemia','Středočeský':'Bohemia','Jihočeský':'Bohemia','Plzeňský':'Bohemia','Karlovarský':'Bohemia','Ústecký':'Bohemia','Liberecký':'Bohemia','Královéhradecký':'Bohemia','Pardubický':'Bohemia',
 'Vysočina':'Moravia','Jihomoravský':'Moravia','Olomoucký':'Moravia','Zlínský':'Moravia','Moravskoslezský':'Moravia'}
WHOLE={'Latvia':'Latvia','Slovenia':'Drava Banovina','Macedonia':'Vardar Macedonia','Kosovo':'Kosovo','Montenegro':'Montenegro','Malta':'Malta','Liechtenstein':'Liechtenstein','Luxembourg':'Luxembourg','Northern Ireland':'Northern Ireland',
 'Azerbaijan':'Azerbaijan','Uganda':'Uganda','Burkina Faso':'Burkina Faso','Guinea':'Guinea','Moldova':'Bessarabia','Qatar':'Qatar','Bahrain':'Bahrain','Kuwait':'Kuwait','United Arab Emirates':'Trucial States','Cyprus':'Cyprus','Northern Cyprus':'Cyprus',
 'Akrotiri Sovereign Base Area':'Cyprus','Dhekelia Sovereign Base Area':'Cyprus','Cyprus No Mans Area':'Cyprus','Aland':'Åland','Guernsey':'Channel Islands','Jersey':'Channel Islands','Isle of Man':'Isle of Man','Gibraltar':'Gibraltar','Monaco':'Monaco',
 'Vatican':'Vatican','San Marino':'San Marino','Andorra':'Andorra','Kazakhstan':None}
SKIP={'Iceland','Faroe Islands','Seychelles','Svalbard','Jan Mayen'}
UK_REG={'Greater London':'London','North West':'North West England','North East':'North East England','Yorkshire and the Humber':'Yorkshire and the Humber','East Midlands':'East Midlands','West Midlands':'West Midlands','East':'East of England',
 'South East':'South East England','South West':'South West England','East Wales':'East Wales','West Wales and the Valleys':'West Wales','Highlands and Islands':'Highlands','South Western':'South West Scotland','Eastern':'East Scotland',
 'North Eastern':'North East Scotland','Southern':'Southern Scotland'}
LIB={'Tripolitania':'An Nuqat al Khams|Az Zawiyah|Tajura\' wa an Nawahi al Arba|Al Jifarah|Al Marqab|Misratah|Surt|Mizdah|Ghadamis','Cyrenaica':'Benghazi|Al Marj|Al Jabal al Akhdar|Al Butnan|Ajdabiya|Al Kufrah|Al Qubbah',
 'Fezzan':'Al Jufrah|Ash Shati\'|Sabha|Murzuq|Ghat|Wadi al Hayaa'}
LIBN={n:k for k,v in LIB.items() for n in v.split('|')}
DISAMB={'Bari','Plateaux','Nord','Northern','Eastern','Southern','Plateau','Ash Sharqiyah','Western','Jura','Savanes','Centre','Al Jawf','Littoral','Central'}   # same name in different countries
unmapped=set(); free_france=[]
for p,g in A1:
    ad=p['admin']; nm=p['name']; CUR[0]=ad
    if ad in SKIP: continue
    o=DEF.get(ad)
    if o is None: unmapped.add(ad); continue
    c=lambda k,gg=g,oo=o,ad=ad: put(k+(' ('+ad+')' if k in DISAMB else ''),gg,oo)
    if ad=='France':
        c(FR_MERGE.get(nm,nm))
    elif ad=='Italy':
        k=IT_MERGE.get(nm,nm); c(k)
    elif ad=='Spain': c({'Baleares':'Balearic Islands','Gerona':'Girona','Lérida':'Lleida','Orense':'Ourense','La Coruña':'A Coruña'}.get(nm,nm))
    elif ad=='Greece': c(GR_NAMES.get(nm,nm))
    elif ad=='Germany': c('DE:'+nm)
    elif ad=='Austria': c(nm if nm!='Wien' else 'Vienna')
    elif ad=='Czech Republic': c(CZ[nm])
    elif ad=='Belgium': c('Brabant' if nm in('Brussels','Flemish Brabant','Walloon Brabant') else nm+(' (BE)' if nm in('Limburg','Luxembourg') else ''))
    elif ad=='Netherlands': c('Overijssel' if nm=='Flevoland' else nm+(' (NL)' if nm=='Limburg' else ''))
    elif ad=='Denmark': c('DK:'+nm)
    elif ad=='Finland': c(FIN_MERGE[nm])
    elif ad=='Croatia': c(HR_MERGE.get(nm,nm))
    elif ad=='Republic of Serbia': c(SRB.get(nm,nm))
    elif ad=='Bosnia and Herzegovina': c('Herzegovina' if nm in('Herzegovina-Neretva','West Herzegovina','Trebinje') else 'Bosnia')
    elif ad=='Ireland': c('Dublin' if nm in('Dublin','Dún Laoghaire–Rathdown','Fingal','South Dublin') else 'Tipperary' if 'Tipperary' in nm else 'Laois' if nm=='Laoighis' else nm)
    elif ad=='United Kingdom':
        if p['type_en'] and 'Unitary District' in p['type_en'] or p['region'] in('Highlands and Islands','South Western','North Eastern','Southern') and p['type_en'] and 'Unitary' in p['type_en']:
            c('UK:S:'+(p['region'] or '?'))
        elif p['region']=='Northern Ireland': c('Northern Ireland')
        else: c('UK:'+(p['region'] or '?'))
    elif ad=='Libya': c(LIBN.get(nm,'Fezzan'))
    elif ad=='Algeria':
        cx=g.centroid; lat,lon=cx.y,cx.x
        south=nm in('Adrar','Béchar','Tindouf','Naâma','El Bayadh','Ghardaïa','Laghouat','Ouargla','El Oued','Illizi','Tamanrasset','Djelfa')
        c('Territoires du Sud' if south else 'Oran' if lon<1.25 else 'Alger' if lon<4.6 else 'Constantine')
    elif ad=='Poland':
        inn=mp(g.intersection(POL1939)); out=mp(g.difference(POL1939))
        if nm in PL_NAMES: put(PL_NAMES[nm],inn,('POL','GER' if PL_NAMES[nm] not in('Białystok',) else 'SOV'))
        if nm=='Warmian-Masurian': put('Warsaw',inn,('POL','GER')); put('East Prussia',out,('GER','GER'))
        elif nm=='Greater Poland': put('Pomerania' if out.centroid.y>=52.8 else 'Lower Silesia',out,('GER','GER'))
        elif nm in('Kuyavian-Pomeranian','Pomeranian'): put('Pomerania',out,('GER','GER'))
        elif nm=='West Pomeranian': put('Pomerania',g,('GER','GER'))
        elif nm=='Lubusz':
            put('Brandenburg',mp(g.intersection(box(0,52.0,40,60))),('GER','GER')); put('Lower Silesia',mp(g.intersection(box(0,0,40,52.0))),('GER','GER'))
        elif nm=='Lower Silesian': put('Lower Silesia',g,('GER','GER'))
        elif nm=='Opole': put('Upper Silesia',g,('GER','GER'))
        elif nm=='Silesian': put('Upper Silesia',out,('GER','GER'))
        elif nm in PL_NAMES: put(PL_NAMES[nm],out,('GER','GER'))
    elif ad=='Lithuania':
        inn=mp(g.intersection(POL1939)); out=mp(g.difference(POL1939))
        put('Wilno',inn,('POL','SOV')); put(nm.replace('s','') if False else 'LT:'+nm,out,o)
    elif ad=='Belarus':
        inn=mp(g.intersection(POL1939)); out=mp(g.difference(POL1939)); base='Minsk' if nm=='City of Minsk' else nm
        pk={'Brest':'Polesie','Grodno':'Nowogródek','Minsk':'Nowogródek','City of Minsk':'Nowogródek','Vitebsk':'Wilno'}.get(nm)
        if pk: put(pk,inn,('POL','SOV'))
        put('BY:'+base,out,o)
    elif ad=='Ukraine':
        pk={'Volyn':'Wołyń','Rivne':'Wołyń',"L'viv":'Lwów',"Ternopil'":'Tarnopol',"Ivano-Frankivs'k":'Stanisławów'}.get(nm)
        if pk: put(pk,g,('POL','SOV'))
        elif nm=='Chernivtsi': put('Bukovina',g,('ROM','SOV'))
        elif nm=='Transcarpathia': put('Carpatho-Ukraine',g,('HUN','HUN'))
        else: put('UA:'+nm,g,o)
    elif ad=='Romania':
        north=nm in('Satu Mare','Maramures','Salaj','Bistrita-Nasaud','Cluj','Mures','Harghita','Covasna','Bihor')
        c(nm if not north else nm,); OWNP[nm]=('ROM','HUN') if north else ('ROM','ROM')
    elif ad=='Hungary':
        c('HU:'+nm)
    elif ad=='Slovakia': c(nm)
    elif ad=='Bulgaria':
        if nm in('Dobrich','Silistra'): put('Southern Dobruja',g,('ROM','BUL'))
        else: c(nm)
    elif ad=='Russia':
        c('East Prussia' if nm=='Kaliningrad' else nm if nm not in('Moskva',) else 'Moscow')
        if nm=='Kaliningrad': OWNP['East Prussia']=('GER','GER')
    elif ad in WHOLE:
        if WHOLE[ad]: c(WHOLE[ad])
        else: c(nm)
    else: c(nm if ad not in('Yemen',) else nm)
if unmapped: log('no owner for (dropped):',sorted(unmapped))
log('keys from more than one country:',{k:sorted(v) for k,v in KEYADM.items() if len(v)>1})
# --- Germany 1939: Prussian provinces and the German states, cut out of the modern Laender
RHINE=P((6.75,51.87),(7.0,51.65),(7.05,51.5),(7.15,51.43),(7.3,51.36),(7.45,51.25),(7.7,51.1),(7.95,50.92),(8.05,50.7),(8.1,50.3),(5.5,50.3),(5.5,51.87))
RHEINHESSEN=box(7.85,49.6,8.5,50.06); PALATINATE=P((7.3,48.9),(8.6,48.9),(8.6,49.6),(8.3,49.66),(7.9,49.72),(7.5,49.77),(7.3,49.5))
BADEN=P((9.85,49.7),(9.62,49.52),(9.35,49.42),(9.1,49.3),(9.0,49.15),(8.9,49.0),(8.78,48.88),(8.75,48.7),(8.55,48.58),(8.45,48.38),(8.45,48.2),(8.55,48.1),(8.75,47.95),(9.0,47.88),(9.25,47.78),(9.45,47.7),(9.7,47.6),(9.7,47.0),(5.0,47.0),(5.0,50.0),(9.85,50.0))
OLDENBURG=P((7.95,52.65),(8.55,52.65),(8.75,52.95),(8.7,53.1),(8.65,53.2),(8.5,53.65),(8.0,53.75),(7.95,53.5))
VORPOMMERN=P((12.4,54.6),(12.4,54.25),(12.9,53.95),(13.3,53.45),(13.7,53.1),(14.3,53.1),(14.3,54.6))
HE='DE:Hessen'; ge=lambda k:'DE:'+k
carve(RHINE,'Rhine Province',('GER','GER'),[ge('Nordrhein-Westfalen')]); move(ge('Nordrhein-Westfalen'),'Westphalia')
carve(RHEINHESSEN,'Hesse',('GER','GER'),[ge('Rheinland-Pfalz')]); carve(PALATINATE,'Palatinate',('GER','GER'),[ge('Rheinland-Pfalz')]); move(ge('Rheinland-Pfalz'),'Rhine Province')
carve(box(8.3,49.4,9.3,50.0),'Hesse',('GER','GER'),[HE]); carve(box(8.5,50.28,9.45,50.85),'Hesse',('GER','GER'),[HE]); move(HE,'Hesse-Nassau')
carve(BADEN,'Baden',('GER','GER'),[ge('Baden-Württemberg')]); move(ge('Baden-Württemberg'),'Württemberg')
carve(VORPOMMERN,'Pomerania',('GER','GER'),[ge('Mecklenburg-Vorpommern')]); move(ge('Mecklenburg-Vorpommern'),'Mecklenburg')
carve(OLDENBURG,'Oldenburg',('GER','GER'),[ge('Niedersachsen')]); move(ge('Niedersachsen'),'Hanover')
carve(box(11.5,51.62,12.6,52.0),'Anhalt',('GER','GER'),[ge('Sachsen-Anhalt')]); move(ge('Sachsen-Anhalt'),'Province of Saxony')
for a,b in [('Bayern','Bavaria'),('Sachsen','Saxony'),('Thüringen','Thuringia'),('Brandenburg','Brandenburg'),('Berlin','Berlin'),('Hamburg','Hamburg'),('Bremen','Bremen'),('Schleswig-Holstein','Schleswig-Holstein'),('Saarland','Saarland')]: move(ge(a),b)
# --- Denmark: the historical parts
dk=lambda k:'DK:'+k
carve(box(14.5,54.9,15.5,55.4),'Bornholm',('DEN','GER'),[dk('Hovedstaden')]); carve(box(11.0,54.4,12.5,54.97),'Lolland-Falster',('DEN','GER'),[dk('Sjaælland')])
carve(box(9.75,54.8,11.0,55.8),'Funen',('DEN','GER'),[dk('Syddanmark')]); move(dk('Syddanmark'),'Southern Jutland'); move(dk('Hovedstaden'),'Zealand'); move(dk('Sjaælland'),'Zealand'); move(dk('Midtjylland'),'Central Jutland'); move(dk('Nordjylland'),'North Jutland')
# --- Latvia: Latgale, Kurzeme, Vidzeme, Zemgale
carve(box(26.0,55.0,29.0,58.2),'Latgale',('LAT','SOV'),['Latvia']); carve(box(20.0,55.0,23.2,58.2),'Kurzeme',('LAT','SOV'),['Latvia'])
carve(P((24.0,56.9),(24.35,56.78),(25.4,56.6),(26.0,56.9),(26.0,58.3),(24.0,58.3)),'Vidzeme',('LAT','SOV'),['Latvia']); move('Latvia','Zemgale')
# --- 1939 borders that differ from today's
MOR=['Guelmim - Es-Semara','Laâyoune - Boujdour - Sakia El Hamra','Oued el Dahab','Souss - Massa - Draâ','Tanger - Tétouan','Gharb - Chrarda - Béni Hssen','Taza - Al Hoceima - Taounate','Oriental','Fès - Boulemane']
SPA=('ESP','ESP')
carve(box(-18.0,20.5,-8.6,27.67),'Western Sahara',SPA,MOR)                                            # Spanish Sahara: the Moroccan-held part is in Natural Earth's Morocco
carve(box(-14.0,27.67,-8.67,28.1),'Cape Juby',SPA,MOR)                                                # Tarfaya strip (Cape Juby)
carve(box(-10.25,29.2,-9.7,29.65),'Ifni',SPA,MOR)                                                       # Ifni enclave
carve(P((-6.5,34.9),(-5.6,34.82),(-4.8,34.78),(-4.2,34.85),(-3.7,35.0),(-3.0,35.0),(-2.35,35.05),(-2.35,35.95),(-6.5,35.95)),'Spanish Morocco',SPA)   # northern protectorate, with Ceuta and Melilla
carve(box(-18.0,26.0,-8.6,27.67),'Saguia el-Hamra',SPA,['Western Sahara']); move('Western Sahara','Río de Oro',SPA)   # Spanish Sahara's two halves
carve(P((28.2,45.2),(29.8,45.2),(30.35,45.6),(30.55,46.05),(30.4,46.25),(29.95,46.5),(29.0,46.55),(28.2,46.1)),'Bessarabia',('ROM','SOV'),['UA:Odessa'])   # Southern Bessarabia (Budjak: Izmail, Cetatea Alba), incl. the Danube delta
carve(P((28.3,48.2),(28.6,48.15),(28.85,47.9),(29.0,47.75),(29.1,47.6),(29.12,47.4),(29.2,47.2),(29.45,47.0),(29.55,46.9),(29.57,46.85),(29.75,46.6),(29.95,46.5),(30.4,46.5),(30.4,48.3),(28.3,48.3)),'UA:Odessa',('SOV','SOV'),['Bessarabia'])   # Transnistria, east of the Dniester, was Soviet (Moldavian ASSR)
carve(P((27.3,57.55),(28.0,57.55),(28.0,57.95),(27.3,57.95)),'Võru',('EST','SOV'),['Pskov'])                                   # Petseri County, Estonian until 1945
carve(P((27.55,56.85),(28.25,56.85),(28.25,57.3),(27.55,57.3)),'Latgale',('LAT','SOV'),['Pskov'])                           # Abrene County (Pytalovo), Latvian until 1945
carve(P((18.45,54.17),(19.3,54.17),(19.3,54.46),(18.45,54.46)),'Danzig',('DAN','GER'))                                  # Free City of Danzig: a free city (own country) until Germany took it on 1 Sept 1939
carve(box(18.7,53.2,21.2,54.6),'East Prussia',('GER','GER'),['Pomerania'])                                               # east bank of the Vistula (Marienburg, Elbing, Marienwerder) belonged to East Prussia
carve(P((20.9,55.05),(22.05,55.05),(22.05,55.8),(20.9,55.8)),'Memelland',('GER','GER'),['LT:Klaipedos'])                     # returned to Germany, March 1939
carve(P((22.0,47.8),(24.75,47.8),(24.75,49.0),(22.0,49.0)),'Carpatho-Ukraine',('HUN','HUN'),['Bukovina','UA:Transcarpathia'])
carve(P((18.3,49.6),(18.95,49.6),(18.95,49.95),(18.3,49.95)),'Silesia',('POL','GER'),['Moravia'])                             # Zaolzie, Polish since Oct 1938
carve(P((28.9,69.0),(31.3,69.0),(31.3,70.4),(28.9,70.4)),'Petsamo',('FIN','SOV'),['Murmansk','Karelia'])
carve(P((27.9,60.0),(32.6,60.0),(32.6,62.3),(27.9,62.3)),'Karelian Isthmus',('FIN','SOV'),['Leningrad','Karelia','City of St. Petersburg'])
carve(P((16.8,47.7),(22.5,47.7),(22.5,48.75),(20.6,48.6),(19.6,48.25),(18.2,48.0),(17.3,48.0)),'Southern Slovakia',('HUN','HUN'),['Nitriansky','Banskobystrický','Košický','Bratislavský','Trnavský','Prešov'])
carve(box(26.2,35.2,29.9,37.9),'Dodecanese',('ITA','GRE'),['South Aegean'])
carve(box(15.12,44.06,15.32,44.18),'Zara',('ITA','ITA'),['Dalmatia'])
JM=P((13.4,46.6),(13.65,46.5),(13.75,46.3),(14.0,46.2),(14.1,46.05),(14.2,45.95),(14.3,45.8),(14.45,45.65),(14.5,45.45),(14.46,45.3),(14.52,44.95),(14.5,44.4),(13.0,44.4),(13.0,46.6))
carve(P((13.0,45.8),(14.6,45.8),(14.6,46.7),(13.0,46.7)).intersection(JM),'Gorizia',('ITA','ITA'),['Drava Banovina']); carve(JM,'Fiume',('ITA','ITA'),['Drava Banovina','Croatia'])
for k in('Pola','Fiume'): OWNP[k]=('ITA','ITA')
for k,c in {'Vatican':'VAT','San Marino':'SMR','Monaco':'MON','Andorra':'AND','Liechtenstein':'LIE'}.items():
    if k in OWNP: OWNP[k]=(c,c)       # independent micro states
carve(box(12.435,41.892,12.468,41.914),'Vatican',('VAT','VAT'))   # the real Vatican is 0.44 km2, too small to hold a tile: drawn about 20x too big so it can be seen and selected
carve(P((29.5,35.0),(31.0,35.0),(31.0,38.0),(29.5,38.0)),'Dodecanese',('ITA','GRE'),['South Aegean'])
# --- Yugoslav provinces: Banovina of Croatia and friends are not modelled; Herzegovina/Bosnia, Slavonia/Croatia/Dalmatia/Lika are the historical regions
for k in('Slavonia','Croatia','Lika','Dalmatia','Bosnia','Herzegovina','Bačka','Banat','Syrmia','Belgrade'): OWNP[k]=('YUG','YUG')
# --- France: the free zone in 1941 (departements mostly south and east of the demarcation line)
FREE=['Alpes-Maritimes','Alpes-de-Haute-Provence','Hautes-Alpes','Ardèche','Ariège','Aude','Aveyron','Bouches-du-Rhône','Cantal','Corrèze','Corsica','Drôme','Gard','Haute-Garonne','Gers','Hérault','Isère','Loire','Haute-Loire','Lot',
 'Lozère','Hautes-Pyrénées','Pyrénées-Orientales','Puy-de-Dôme','Rhône','Haute-Savoie','Savoie','Tarn','Tarn-et-Garonne','Var','Vaucluse','Haute-Vienne','Creuse','Allier','Lot-et-Garonne','Pyrénées-Atlantiques','Dordogne','Landes']
for k in FREE:
    if k in OWNP: OWNP[k]=('FRA','VIC')
# --- Italian East Africa and Libya names, and other tidy-ups
for k,o in list(OWNP.items()):
    if k in('Wilno','Nowogródek','Polesie','Wołyń','Lwów','Tarnopol','Stanisławów','Białystok'): OWNP[k]=('POL','SOV')
for k in('Poznań','Pomorze','Silesia','Kraków','Warsaw','Łódź','Kielce','Lublin','Rzeszów'): OWNP[k]=('POL','GER')
# unify Hungarian counties: fold the 'urban counties' into the county that surrounds them
def absorb_urban():
    hu=[k for k in REGS if k.startswith('HU:')]
    cap='HU:Budapest'
    urb=[]
    for p,g in A1:
        if p['admin']=='Hungary' and p['type_en'] in('Urban county',): urb.append('HU:'+p['name'])
    for u in urb:
        if u not in REGS: continue
        best=None;bl=0
        for k in REGS:
            if k==u or not k.startswith('HU:') or k in urb: continue
            l=REGS[u].boundary.intersection(REGS[k].buffer(.0005)).length
            if l>bl: best=k;bl=l
        if best: move(u,best)
absorb_urban()
for k in [k for k in REGS if k.startswith('HU:')]: move(k,k[3:].replace('Capital City','Budapest'),('HUN','HUN'))
for k in [k for k in REGS if k.startswith('BY:')]: move(k,k[3:] if k[3:] not in('Minsk','Vitebsk','Brest','Grodno') else k[3:],('SOV','SOV'))
for k in [k for k in REGS if k.startswith('UA:')]: move(k,k[3:],('SOV','SOV'))
for k in [k for k in REGS if k.startswith('LT:')]: move(k,k[3:],('LIT','SOV'))
for k in [k for k in REGS if k.startswith('UK:')]:
    if k.startswith('UK:S:'): move(k,UK_REG.get(k[5:],'Scotland'),('GBR','GBR'))
    else: move(k,UK_REG.get(k[3:],k[3:]),('GBR','GBR'))
# ---------------------------------------------------------------- merge county-sized provinces into Munster-sized ones
# Natural Earth admin-1 units are counties, departements, districts... far smaller than a HOI4 state. Here neighbouring provinces
# of the same country and the same owner (in BOTH scenarios) are joined until they are about the size of Munster (Ireland).
#   * a province that is already big (>= MERGE_BELOW: Bavaria, the Prussian provinces, Russian oblasts, ...) is left alone
#   * a merged province never grows past MERGE_CAP, so one huge neighbour does not swallow every small one
#   * only provinces that touch by land are joined, so islands, enclaves and city-states (Malta, Gibraltar, Berlin...) stay small
#   * a merged province takes the name of the member that holds its biggest city
MUNSTER=24600.0               # km2, the target size
MERGE_BELOW=.70*MUNSTER       # at or above this a province counts as 'done'
MERGE_CAP=1.40*MUNSTER        # hard ceiling for a merged province
KEEP_ALONE={'Berlin','Hamburg','Bremen','Vienna','London','Danzig','Memelland','Budapest','Belgrade','Malta','Gibraltar','Monaco','Liechtenstein','Luxembourg',
 'Vatican','San Marino','Andorra','Channel Islands','Isle of Man','Åland','Cyprus','Northern Ireland','Carpatho-Ukraine','Petsamo','Zara','Dodecanese'}
IRL_PROV={'Munster':['Cork','Kerry','Limerick','Tipperary','Waterford','Clare'],
 'Leinster':['Dublin','Wicklow','Wexford','Carlow','Kilkenny','Laois','Offaly','Kildare','Meath','Louth','Westmeath','Longford'],
 'Connacht':['Galway','Mayo','Roscommon','Sligo','Leitrim','Donegal','Cavan','Monaghan']}   # the three Ulster counties of the Free State go with Connacht
SCOTLAND={'Highlands','South West Scotland','East Scotland','North East Scotland','Southern Scotland','Scotland'}; WALES={'West Wales','East Wales'}   # England, Scotland and Wales are never joined
MERGE_RENAME={'Seine':'Île-de-France'}
def _km2(g): return sum(p.area*111.32*math.cos(math.radians(p.centroid.y))*110.57 for p in polys(g))
def merge_provinces():
    fixed=set(); n0=len(REGS)
    for nm,ms in IRL_PROV.items():
        ms=[m for m in ms if m in REGS]
        if not ms: continue
        o=OWNP[ms[0]]; g=mp(unary_union([REGS.pop(m) for m in ms]))
        for m in ms: OWNP.pop(m,None)
        REGS[nm]=g; OWNP[nm]=o; fixed.add(nm)
    keys=[k for k in REGS if k in OWNP and not REGS[k].is_empty]
    area={k:_km2(REGS[k]) for k in keys}
    a1t=STRtree([g for _,g in A1]); ctry={}; reg={}
    for k in keys:
        big=max(polys(REGS[k]),key=lambda p:p.area); idx=a1t.query(big.representative_point(),predicate='intersects')
        pp=A1[int(idx[0])][0] if len(idx) else {}
        ctry[k]=pp.get('admin','?'); reg[k]=(pp.get('admin'),pp.get('region')) if pp.get('region') else None
        if k in SCOTLAND: ctry[k]='Scotland'
        elif k in WALES: ctry[k]='Wales'
        elif ctry[k]=='United Kingdom' and k!='Northern Ireland': ctry[k]='England'
    ct=STRtree([REGS[k] for k in keys]); nb={k:{} for k in keys}
    for i,k in enumerate(keys):
        for j in ct.query(REGS[k].buffer(.003)):
            j=int(j)
            if j<=i: continue
            k2=keys[j]; L=REGS[k].boundary.intersection(REGS[k2].buffer(.003)).length
            if L>.05: nb[k][k2]=nb[k2][k]=L
    G={k:{'m':[k],'a':area[k],'o':OWNP[k],'c':ctry[k],'r':reg[k],'f':k in KEEP_ALONE or k in fixed} for k in keys}
    def run(below,cap,anyc=False):
        stuck=set()
        while True:
            cand=[g for g,v in G.items() if not v['f'] and v['a']<below and g not in stuck]
            if not cand: break
            g=min(cand,key=lambda x:G[x]['a']); me=G[g]; tot=sum(nb[g].values()) or 1; best=None; bs=-9
            for h,L in nb[g].items():
                H=G[h]
                if H['f'] or H['o']!=me['o'] or (H['c']!=me['c'] and not anyc) or me['a']+H['a']>cap: continue
                s=L/tot-.6*abs(me['a']+H['a']-MUNSTER)/MUNSTER+(.6 if me["r"] and me["r"]==H["r"] else 0)
                if s>bs: bs=s;best=h
            if best is None: stuck.add(g); continue
            h=best; H=G[h]; me['m']+=H['m']; me['a']+=H['a']; del G[h]
            if (H['a']>me['a']-H['a']): me['r']=H['r']
            for x,L in nb.pop(h).items():
                if x==g: continue
                nb[x].pop(h,None); nb[x][g]=nb[x].get(g,0)+L; nb[g][x]=nb[g].get(x,0)+L
            nb[g].pop(h,None)
    run(MERGE_BELOW,MERGE_CAP)            # pass 1: join neighbours up to ~Munster size
    run(.40*MUNSTER,1.65*MUNSTER)         # pass 2: stragglers under ~10,000 km2 may join a bit past the usual ceiling
    run(300,1.65*MUNSTER,True)            # pass 3: slivers at the map edge join any neighbour with the same owner
    # name each merged province after the member that holds its biggest city
    pop={}
    cities=[(f['properties']['pop_max'] or 0,Point(*f['geometry']['coordinates'][:2])) for f in feats('ne_10m_populated_places_simple')]
    mt=STRtree([REGS[k] for k in keys])
    for pm,pt in cities:
        if pm<100000: continue
        for j in mt.query(pt,predicate='intersects'):
            k=keys[int(j)]; pop[k]=max(pop.get(k,0),pm)
    report=[]
    for g,v in G.items():
        if len(v['m'])==1: continue
        name=max(v['m'],key=lambda k:(pop.get(k,0),area[k])); name=MERGE_RENAME.get(name,name)
        o=OWNP[g]; geom=mp(unary_union([REGS[m] for m in v['m']]))
        for m in v['m']: REGS.pop(m); OWNP.pop(m,None)
        REGS[name]=geom; OWNP[name]=o; report.append((name,round(v['a']),[m for m in v['m']]))
    log('merged provinces:',n0,'->',len(REGS))
    return report
MERGE_REPORT=merge_provinces()
REG=[k for k in REGS if k in OWNP and not REGS[k].is_empty]
miss=[k for k in REGS if k not in OWNP]
if miss: log('provinces without an owner (dropped):',miss)
RI={k:i for i,k in enumerate(REG)}
LAND=mp(unary_union([g for _,g in A1]))
log('provinces',len(REG))
if '--provinces-only' in sys.argv:
    for k in sorted(REG): print(k,OWNP[k],round(_km2(REGS[k])))
    for r in MERGE_REPORT: print('MERGED',r[0],r[1],r[2])
    sys.exit()
# ---------------------------------------------------------------- rivers (they are hard tile borders: see the bank split below)
# Class 1 = great rivers, 2 = major, 3 = minor. Every class cuts tiles; the class only sets how thick the line is drawn.
# Natural Earth has no Aisne/Oise/Sambre/Lys/Yser/Isonzo/Tagliamento/Tisza/Vah/Vardar and only part of the Neman and Saone, so those
# are drawn by hand from waypoints (lon,lat). A hand river only adds the stretches Natural Earth does not already have.
HAND={
 'Aisne':(2,[(5.17,48.95),(4.9,49.1),(4.75,49.3),(4.37,49.51),(4.0,49.45),(3.6,49.4),(3.32,49.38),(2.83,49.42)]),
 'Oise':(2,[(4.3,49.98),(4.08,49.92),(3.63,49.88),(3.35,49.75),(3.0,49.58),(2.83,49.42),(2.62,49.38),(2.48,49.26),(2.3,49.12),(2.1,49.05),(2.08,49.0)]),
 'Sambre':(2,[(3.7,49.97),(3.69,50.12),(3.97,50.28),(4.2,50.33),(4.44,50.41),(4.65,50.43),(4.87,50.46)]),
 'Lys':(2,[(2.25,50.45),(2.4,50.64),(2.65,50.65),(2.88,50.68),(3.12,50.78),(3.27,50.83),(3.53,50.98),(3.72,51.05)]),
 'Yser':(3,[(2.4,50.82),(2.55,50.92),(2.65,50.98),(2.78,51.05),(2.75,51.14),(2.73,51.15)]),
 'Dyle':(3,[(4.5,50.6),(4.6,50.72),(4.65,50.82),(4.7,50.88),(4.6,50.97),(4.5,51.03)]),
 'Albert Canal':(2,[(5.6,50.66),(5.67,50.75),(5.68,50.83),(5.6,50.9),(5.34,50.94),(5.05,51.08),(4.83,51.18),(4.55,51.23),(4.42,51.25)]),
 'Kiel Canal':(2,[(9.13,53.9),(9.3,53.95),(9.4,54.1),(9.6,54.27),(9.85,54.3),(10.05,54.33),(10.15,54.37)]),
 'Isonzo':(2,[(13.7,46.4),(13.75,46.25),(13.73,46.18),(13.65,46.0),(13.62,45.94),(13.5,45.85),(13.55,45.72)]),
 'Tagliamento':(2,[(12.5,46.43),(12.75,46.42),(13.0,46.4),(13.02,46.28),(12.95,46.18),(12.95,46.0),(13.0,45.85),(13.08,45.65)]),
 'Tisza':(1,[(24.2,48.1),(23.4,47.95),(22.5,48.1),(22.3,48.13),(21.41,48.12),(20.76,47.62),(20.2,47.17),(20.14,46.71),(20.15,46.25),(20.29,45.2)]),
 'Vah':(2,[(19.7,49.1),(19.3,49.07),(18.74,49.22),(18.4,49.1),(18.04,48.89),(17.83,48.59),(17.9,48.25),(18.12,47.76)]),
 'Vardar':(2,[(20.95,42.0),(21.43,42.0),(21.77,41.72),(22.0,41.5),(22.4,41.2),(22.55,41.0),(22.65,40.7),(22.72,40.5)]),
 'Neman':(1,[(26.5,53.2),(26.75,53.48),(26.1,53.55),(25.4,53.5),(24.5,53.4),(24.0,53.55),(23.83,53.68),(23.97,54.02),(24.05,54.4),(23.9,54.9),(23.2,55.0),(22.77,55.07),(22.2,55.08),(21.88,55.08),(21.5,55.2),(21.3,55.3)]),
 'Saone':(2,[(5.95,48.0),(5.9,47.7),(5.6,47.45),(5.38,47.2),(5.05,47.0),(4.85,46.78),(4.85,46.5),(4.83,46.3),(4.75,45.98),(4.83,45.76)]),
 'Sozh':(2,[(31.7,54.6),(31.3,53.9),(30.9,53.5),(30.8,52.9),(30.9,52.3),(30.8,51.95)]),
 'Orontes':(2,[(36.4,34.15),(36.6,34.5),(36.72,34.72),(36.75,35.13),(36.4,35.8),(36.2,36.2),(36.0,36.12),(35.95,36.05)]),
 'Medjerda':(2,[(7.95,36.28),(8.43,36.45),(8.78,36.5),(9.2,36.55),(9.6,36.65),(10.0,36.95),(10.25,37.09)]),
 'Mersey':(3,[(2.15,53.41),(2.4,53.4),(2.58,53.39),(2.75,53.33),(2.95,53.36),(3.02,53.43)]),
 'Clyde':(3,[(3.8,55.6),(3.78,55.67),(4.0,55.8),(4.25,55.86),(4.5,55.93),(4.7,55.95)]),
 'Moselle':(2,[(6.85,47.9),(6.59,48.02),(6.45,48.17),(6.29,48.37),(6.1,48.55),(5.89,48.68),(6.13,48.76),(6.05,48.91),(6.18,49.12),(6.16,49.36),(6.37,49.47),(6.64,49.75),(7.07,49.92),(7.12,49.95),(7.17,50.15),(7.59,50.36)]),
 'Meuse':(1,[(5.5,47.97),(5.6,48.15),(5.69,48.35),(5.67,48.55),(5.59,48.76),(5.54,48.89),(5.42,49.0),(5.38,49.16),(5.28,49.3),(5.19,49.49),(5.08,49.6),(4.95,49.7),(4.72,49.76),(4.65,49.92),(4.7,50.0),(4.83,50.13),(4.91,50.26),(4.87,50.46),(5.05,50.48),(5.23,50.52),(5.4,50.58),(5.57,50.64),(5.69,50.74),(5.69,50.85),(5.75,51.0),(5.87,51.1),(5.99,51.19),(6.1,51.28),(6.17,51.37),(6.1,51.5),(6.0,51.62),(5.97,51.7)]),
 'Somme':(2,[(3.29,49.85),(2.93,49.93),(2.3,49.89),(1.83,50.1),(1.58,50.2)]),
}
from shapely.ops import nearest_points
from shapely.geometry import LineString
def lines_of(g):
    if g is None or g.is_empty: return []
    if g.geom_type=='LineString': return [g]
    if g.geom_type in('MultiLineString','GeometryCollection'): return [l for x in g.geoms for l in lines_of(x)]
    return []
def dedup(lines,acc,tol=.012,snap=.06,minlen=.03):
    """keep only the parts of `lines` farther than tol from what is already in acc, and join their ends to the river they flow into"""
    out=[]
    for l in sorted(lines,key=lambda l:-l.length):
        pieces=[l]
        if acc:
            near=[a for a in acc if a.bounds[0]<l.bounds[2]+tol and a.bounds[2]>l.bounds[0]-tol and a.bounds[1]<l.bounds[3]+tol and a.bounds[3]>l.bounds[1]-tol]
            if near:
                U=unary_union(near);pieces=lines_of(l.difference(U.buffer(tol)))
                fixed=[]
                for p in pieces:
                    if p.length<minlen:continue
                    cs=list(p.coords)
                    for end in (0,-1):
                        pt=Point(cs[end]);d=U.distance(pt)
                        if 0<d<snap:
                            q=nearest_points(U,pt)[0];cs=[(q.x,q.y)]+cs if end==0 else cs+[(q.x,q.y)]
                    fixed.append(LineString(cs))
                pieces=fixed
        for p in pieces:
            if p.length>=minlen:out.append(p);acc.append(p)
    return out
def river_lines():
    cl=[];eu=[];ranked=[]
    for f in feats('ne_10m_rivers_lake_centerlines'):
        p=f['properties']
        if p.get('featurecla')!='River':continue
        g=shape(f['geometry']).intersection(BOX)
        for l in lines_of(g):
            r=p.get('scalerank') or 99;y=l.centroid.y
            if r>(9 if y>=37 else 7):continue          # south of 37N only the big rivers: deserts are full of dry wadis
            ranked.append((1 if r<=4 else 2 if r<=7 else 3,r,l))
    for f in feats('ne_10m_rivers_europe'):
        p=f['properties']
        if p.get('featurecla')=='Lake Centerline' or not p.get('name'):continue   # unnamed streams are too small to matter
        for l in lines_of(shape(f['geometry']).intersection(BOX)): ranked.append((3,p.get('scalerank') or 12,l))
    acc=[];RV={1:[],2:[],3:[]};added=[]
    for c in (1,2,3):    # great rivers claim their course first, minor ones only add what is missing
        RV[c]=dedup([l.simplify(.004) for cc,r,l in sorted(ranked,key=lambda t:t[1]) if cc==c],acc)
    U0=unary_union(acc)
    for nm,(hc,pts) in HAND.items():   # hand-drawn last: add only the continuous stretches Natural Earth does not cover
        hl=LineString(pts);n=max(2,int(hl.length/.02));sm=[hl.interpolate(i/n,normalized=True) for i in range(n+1)]
        unc=[U0.distance(p)>.05 for p in sm];new=[];i=0
        while i<len(sm):
            if not unc[i]: i+=1;continue
            j=i
            while j+1<len(sm) and unc[j+1]: j+=1
            if (j-i+1)*.02>=.25:    # a gap shorter than ~25 km is left alone
                a=max(0,i-1);b=min(len(sm)-1,j+1);new.append(LineString([(p.x,p.y) for p in sm[a:b+1]]))
            i=j+1
        RV[hc]+=new;acc+=new
        if new: added.append((nm,round(sum(l.length for l in new),2),round(hl.length,2)))
    return RV,added
RV,HAND_ADDED=river_lines()
log('rivers',{c:(len(v),round(sum(l.length for l in v))) for c,v in RV.items()})
log('hand-drawn rivers adding something (name, new deg, full deg):',HAND_ADDED)
RIVU=unary_union([l for v in RV.values() for l in v]).intersection(LAND)
# ---- rivers that serve as borders: the admin-1 border runs a kilometre or two off the river, so move the border vertices onto the river.
RIVL=lines_of(RIVU);RT=STRtree(RIVL);SNAPD=.035
def _snap_pt(x,y):
    p=Point(x,y);best=None;bd=SNAPD
    for i in RT.query(p.buffer(SNAPD)):
        l=RIVL[int(i)];d=l.distance(p)
        if d<bd:bd=d;best=l
    if best is None:return (x,y)
    q=best.interpolate(best.project(p));return (q.x,q.y)
def _snap_ring(cs):
    out=[]
    for a,b in zip(cs[:-1],cs[1:]):
        L=math.hypot(b[0]-a[0],b[1]-a[1]);pts=[a]
        if L>.007 and len(RT.query(LineString([a,b]).buffer(SNAPD)))>0:
            n=int(L/.005)+1;pts=[(a[0]+(b[0]-a[0])*i/n,a[1]+(b[1]-a[1])*i/n) for i in range(n)]
        out+=pts
    out.append(cs[-1]);sn=[_snap_pt(*p) for p in out];res=[sn[0]]
    for p in sn[1:]:
        if p!=res[-1]:res.append(p)
    return res
def snap_poly(pg):
    out=[]
    for p in polys(pg):
        ext=_snap_ring(list(p.exterior.coords))
        if len(ext)<4:continue
        hol=[h for h in (_snap_ring(list(r.coords)) for r in p.interiors) if len(h)>=4]
        out.append(valid(Polygon(ext,hol)))
    return mp(unary_union(out)) if out else pg
_moved=0
for _k in REG:
    _o=REGS[_k];_n=snap_poly(_o)
    if not _n.is_empty and abs(_n.area-_o.area)<.05*_o.area+1e-4:REGS[_k]=_n;_moved+=1
log('province borders snapped to rivers for',_moved,'of',len(REG),'provinces')
try:    # snapping added a vertex every ~0.5 km along river borders; thin them out again without opening gaps between neighbours
    _simp=list(shapely.coverage_simplify([REGS[k] for k in REG],.0025))
    _v0=sum(len(r.coords) for k in REG for p in polys(REGS[k]) for r in [p.exterior]+list(p.interiors))
    for _k,_g in zip(REG,_simp):
        if not _g.is_empty:REGS[_k]=mp(_g)
    _v1=sum(len(r.coords) for k in REG for p in polys(REGS[k]) for r in [p.exterior]+list(p.interiors))
    log('coverage_simplify after snapping: vertices',_v0,'->',_v1)
except Exception as e: log('coverage_simplify after snapping failed',e)
if '--rivers-only' in sys.argv:
    json.dump({'rv':{c:[list(l.coords) for l in v] for c,v in RV.items()},'land':[list(p.exterior.coords) for p in polys(LAND)]},open('/tmp/rv.json','w'))
    sys.exit()
# ---- seeds (dart throwing; smaller tiles in the crowded west, bigger in deserts and the far east)
def rf(x,y):
    if y>=36.5: r=.30+.011*max(0,x-12)+.02*max(0,y-58)
    else: r=min(1.25,.6+.018*(37-y))
    return K*min(1.3,r)
def seeds_for(pg):
    pts=[];grid={};cs=1.4
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
    if not pts:   # a province too small for the loop above still gets one tile
        big=max(polys(pg),key=lambda p:p.area);rp=big.representative_point();pts.append((rp.x,rp.y))
    return pts
def cells(pg,pts):
    if len(pts)==1: return {0:pg}
    vd=voronoi_diagram(MultiPoint(pts),envelope=box(*pg.bounds).buffer(2));tree=STRtree([Point(p) for p in pts]);out={}
    for c in vd.geoms:
        idx=tree.query(c,predicate='contains')
        if len(idx)==0:continue
        try: piece=mp(c.intersection(valid(shapely.clip_by_rect(pg,*c.bounds))))
        except Exception: continue
        if not piece.is_empty: out[int(idx[0])]=piece
    return out
# ---- bank split: a province is cut along its rivers into "banks" and every tile is made inside one bank, so no tile can straddle a river.
from shapely.ops import polygonize,linemerge
RIVL=lines_of(RIVU);RT=STRtree(RIVL)
def amin(f):
    c=f.representative_point();return .25*rf(c.x,c.y)**2     # a bank smaller than ~1/4 of a tile is a sliver, not a bank
def banks_of(pg):
    cand=[RIVL[int(i)] for i in RT.query(pg)]
    pb=pg.boundary;pbb=pb.buffer(3e-7)
    Ls=[l for l in lines_of(unary_union(cand).intersection(pg).difference(pb.buffer(.004))) if l.length>.01] if cand else []   # stretches on the border are the border itself
    if not Ls: return [pg],0,0
    def classify(lines):
        bd=[];rv=[]
        for e in lines_of(shapely.union_all(lines,grid_size=1e-7)):    # grid noding: a river end that is a hair off the border still meets it
            (bd if pbb.contains(e.interpolate(.5,normalized=True)) else rv).append(e)
        return bd,rv
    bd,rv=classify([pb]+Ls)
    # a river that ends inside the province (a source) cannot close a bank: continue it in a straight line to the nearest border or river,
    # but only up to ~2.5 tiles. The continuation is an ordinary tile border (not a river: it gets no crossing penalty).
    cnt={}
    for e in bd+rv:
        for c in (e.coords[0],e.coords[-1]):
            k=(round(c[0],6),round(c[1],6));cnt[k]=cnt.get(k,0)+1
    ext=[];brg=[]
    for e in rv:
        for end in (0,-1):
            c=e.coords[end]
            if cnt[(round(c[0],6),round(c[1],6))]>1 or pbb.contains(Point(c)):continue
            if pb.distance(Point(c))<=.0075:   # the river stops just short of the border (it was trimmed there): join it to the border
                h=nearest_points(pb,Point(c))[0];seg=LineString([c,(h.x,h.y)])
                if seg.length>1e-6:brg.append(seg)
                continue
            d=min(e.length*.9,.05);ln=e.length
            q=e.interpolate(d if end==0 else ln-d);dx,dy=c[0]-q.x,c[1]-q.y;n=math.hypot(dx,dy)
            if n<1e-9:continue
            cap=3.5*rf(c[0],c[1]);ray=LineString([c,(c[0]+dx/n*cap,c[1]+dy/n*cap)])
            others=unary_union([pb]+[x for x in rv if x is not e]);hit=ray.intersection(others)
            pts=[g for x in ([hit] if hit.geom_type!='GeometryCollection' else list(hit.geoms)) for g in ([x] if x.geom_type=='Point' else [Point(z) for z in x.coords][:1] if x.geom_type=='LineString' else list(getattr(x,'geoms',[])))]
            pts=[g for g in pts if g.geom_type=='Point' and Point(c).distance(g)>1e-5]
            if pts:
                h=min(pts,key=lambda g:Point(c).distance(g));seg=LineString([c,(h.x,h.y)])
                if pg.buffer(1e-7).contains(seg):ext.append(seg)
    if ext or brg: bd,rv=classify([pb]+Ls+ext+brg)
    EU=unary_union(ext) if ext else None;LU=unary_union(Ls+brg)     # bridges join a river end to the border: they are part of the river cut and are never merged away
    isx=lambda e:(EU is not None) and e.interpolate(.5,normalized=True).distance(EU)<4e-7 and e.interpolate(.5,normalized=True).distance(LU)>4e-7
    keep=list(rv);xf=[isx(e) for e in keep];dropped=0
    if os.environ.get('DBGB'):print('   Ls',[(round(l.length,3),round(l.coords[0][0],3),round(l.coords[0][1],3),round(l.coords[-1][0],3),round(l.coords[-1][1],3)) for l in Ls],'ext',len(ext),'bd',len(bd),'rv',len(rv),'faces0',[round(f.area,4) for f in polygonize(bd+keep) if pg.contains(f.representative_point())])
    mids=lambda e:e.interpolate(.5,normalized=True)
    for it in range(80):
        faces=[f for f in polygonize(bd+keep) if pg.contains(f.representative_point())]
        small=sorted([f for f in faces if f.area<amin(f)],key=lambda f:f.area)
        act=False
        for f in small:
            fb=f.boundary.buffer(3e-7)
            adjk=[q for q,e in enumerate(keep) if fb.contains(mids(e))]
            xs=[q for q in adjk if xf[q]]
            if xs:    # merge with the neighbour across a continuation line (it is not a river, so nothing is lost)
                q=max(xs,key=lambda q:keep[q].length);keep.pop(q);xf.pop(q);act=True;break
            if f.area<max(.0015,.008*rf(*f.representative_point().coords[0])**2) and adjk:     # a true sliver: the river is within metres of the border
                if any(fb.contains(mids(e)) for e in bd): drop=adjk
                else: drop=[max(adjk,key=lambda q:keep[q].length)]
                dropped+=len(drop);keep=[e for q,e in enumerate(keep) if q not in drop];xf=[x for q,x in enumerate(xf) if q not in drop];act=True;break
        if not act:break
    else:
        faces=[f for f in polygonize(bd+keep) if pg.contains(f.representative_point())]
    if not faces: return [pg],0,0
    return faces,len(keep),dropped
for _a in sys.argv:
    if _a.startswith('--debug=') or _a.startswith('--debugxy='):
        _ks=_a[8:].split(',') if _a.startswith('--debug=') else [k for k in REG if REGS[k].contains(Point(*map(float,_a[10:].split(','))))]
        for _k in _ks:
            _pg=REGS[_k];cand=[RIVL[int(i)] for i in RT.query(_pg)];Ls=lines_of(unary_union(cand).intersection(_pg)) if cand else []
            for l in Ls: print('  piece',[ (round(x,3),round(y,3)) for x,y in (l.coords[0],l.coords[-1])],'len %.2f'%l.length,'ends to border %.5f %.5f'%(_pg.boundary.distance(Point(l.coords[0])),_pg.boundary.distance(Point(l.coords[-1]))))
            _pb=_pg.boundary;_pbb=_pb.buffer(1e-7);_ed=lines_of(unary_union([_pb]+Ls));_bd=[e for e in _ed if _pbb.contains(e.interpolate(.5,normalized=True))];_rv=[e for e in _ed if e not in _bd]
            for e in _ed: print('   E',len(e.coords),[(round(x,4),round(y,4)) for x,y in (e.coords[0],e.coords[-1])],'bd' if e in _bd else 'rv','len %.3f'%e.length)
            print('  edges',len(_ed),'bd',len(_bd),'rv',len(_rv),'polygonize',[round(f.area,3) for f in polygonize(_bd+_rv)],'type',_pg.geom_type)
            _f,_nk,_nd=banks_of(_pg)
            print('DEBUG',_k,'area %.3f'%_pg.area,'river pieces',len(Ls),'len %.2f'%sum(l.length for l in Ls),'faces',len(_f),[round(f.area,3) for f in _f],'kept edges',_nk,'dropped',_nd,'amin %.4f'%amin(_pg))
        sys.exit()
PR=[];NBANK=0;NCUT=0;NDROP=0;NKEEP=0
for ri,k in enumerate(REG):
    pg=REGS[k];faces,nk,nd=banks_of(pg);NBANK+=len(faces);NDROP+=nd;NKEEP+=nk;NCUT+=(len(faces)>1);got=0
    for fc in faces:
        pts=seeds_for(fc);c1=cells(fc,pts);new=[]
        for i,p in enumerate(pts):
            pc=c1.get(i)
            if pc is None:new.append(p);continue
            cen=pc.centroid;new.append((cen.x,cen.y) if pc.contains(cen) else p)
        gf=0
        for i,pc in cells(fc,new).items():
            ps=[p for p in polys(pc) if p.area>=.0015]
            if ps:
                pc=MultiPolygon(ps);rp=pc.representative_point();PR.append({'r':ri,'poly':pc,'c':(rp.x,rp.y)});got+=1;gf+=1
        if not gf:   # tiny bank: keep it as a single tile
            pc=mp(fc);rp=pc.representative_point();PR.append({'r':ri,'poly':pc,'c':(rp.x,rp.y)});got+=1
    if not got:
        pc=mp(pg);rp=pc.representative_point();PR.append({'r':ri,'poly':pc,'c':(rp.x,rp.y)})
log('banks',NBANK,'provinces cut by a river',NCUT,'| river edges kept',NKEEP,'dropped as slivers',NDROP)

# ---- second pass: a river that ends inside a bank (a tributary's source) never closes a bank, so it still runs through tile interiors.
# Any tile with river in its interior is split along that river; a free river end is extended in a straight line to the tile edge.
def river_inside(part):
    cand=[RIVL[int(i)] for i in RT.query(part)]
    if not cand:return []
    U=unary_union(cand).intersection(part)
    if U.is_empty or U.intersection(part.buffer(-6e-4)).length<1e-3:return []      # nothing in the interior: the river is on the tile edge
    return [l for l in lines_of(U) if l.length>1e-3]
def split_part(part,force=False):
    ls=river_inside(part)
    if not ls:return [part]
    global NTRY;NTRY+=1
    cut=[];pb=part.boundary
    for l in ls:
        cs=list(l.coords)
        for end in (0,-1):
            p=Point(cs[end])
            if pb.distance(p)<1e-6 or any(o is not l and o.distance(p)<1e-6 for o in ls):continue
            a=cs[end];b=cs[1] if end==0 else cs[-2];dx,dy=a[0]-b[0],a[1]-b[1];n=math.hypot(dx,dy)
            if n==0:continue
            far=(a[0]+dx/n*1.5,a[1]+dy/n*1.5);ray=LineString([a,far]).intersection(part)
            seg=[g for g in lines_of(ray) if g.distance(p)<1e-6]
            if seg:
                ex=list(seg[0].coords);ex=ex if Point(ex[0]).distance(p)<1e-6 else ex[::-1]
                cs=(ex[::-1][:-1]+cs) if end==0 else (cs+ex[1:])
        cut.append(LineString(cs))
    faces=[f for f in polygonize(lines_of(shapely.union_all([pb]+cut,grid_size=1e-7))) if part.buffer(1e-6).contains(f.representative_point())]
    if force:return faces if len(faces)>=2 else [part]
    if len(faces)<2 or any(f.area<max(.02*part.area,.002) for f in faces):
        global NSKIP;NSKIP+=1
        if len(faces)<2:
            global NF1;NF1+=1
            if NF1<=6:print('NOSPLIT',[(round(x,3),round(y,3)) for x,y in list(part.representative_point().coords)],'river pieces',len(ls),'cut',[(round(c.coords[0][0],3),round(c.coords[0][1],3),round(c.coords[-1][0],3),round(c.coords[-1][1],3),round(pb.distance(Point(c.coords[0])),5),round(pb.distance(Point(c.coords[-1])),5)) for c in cut][:3],'faces',len(faces))
        return [part]     # skims the tile edge: leave it, the river is on the border anyway
    return faces
NSPLIT=0;NMERGE=0;NTRY=0;NSKIP=0;NF1=0;RS={}
def mkt(r,poly):
    rp=poly.representative_point();return {'r':r,'poly':poly,'c':(rp.x,rp.y)}
for _pass in range(3):
    for i in range(len(PR)):      # a tile made of several pieces becomes one tile per piece if a river runs through any of them
        pr=PR[i];ps=polys(pr['poly'])
        if len(ps)>1 and any(river_inside(p) for p in ps if p.area>.003):
            PR[i]=mkt(pr['r'],MultiPolygon([ps[0]]))
            for p in ps[1:]:PR.append(mkt(pr['r'],MultiPolygon([p])))
    tree=STRtree([p['poly'] if p else Point(0,0) for p in PR]);changed=0
    for i in range(len(PR)):
        pr=PR[i]
        if pr is None:continue
        parts=polys(pr['poly'])
        if len(parts)!=1 or parts[0].area<=.003:continue
        part=parts[0];faces=split_part(part,force=True)
        if len(faces)<2:RS['nofaces']=RS.get('nofaces',0)+1;continue
        lim=max(.02*part.area,.002);small=[f for f in faces if f.area<lim];big=[f for f in faces if f.area>=lim]
        if not big:RS['nobig']=RS.get('nobig',0)+1;continue
        extra=[];merges=[];ok=True
        for f in small:      # a sliver between the river and the tile edge joins the neighbouring tile of the same province on that side
            best=None;bl=0
            for j in tree.query(f.buffer(1e-4)):
                j=int(j)
                if j==i or PR[j] is None or PR[j]['r']!=pr['r']:continue
                l=f.boundary.intersection(PR[j]['poly'].buffer(1e-6)).length
                if l>bl:best=j;bl=l
            if best is not None and bl>.4*f.boundary.length*.3:
                u=unary_union([PR[best]['poly'],f])
                if u.geom_type=='Polygon' or (u.geom_type=='MultiPolygon' and len(u.geoms)==len(polys(PR[best]['poly']))):merges.append((best,mp(u)));continue
            if f.area>=.0015:extra.append(f)
            else:ok=False;break
        if not ok:RS['tinysliver']=RS.get('tinysliver',0)+1;continue
        for j,u in merges:PR[j]=mkt(PR[j]['r'],u);NMERGE+=1
        PR[i]=mkt(pr['r'],MultiPolygon([big[0]]))
        for f in big[1:]+extra:PR.append(mkt(pr['r'],MultiPolygon([f])))
        NSPLIT+=1;changed+=1
    log('river split pass',_pass,'tiles split',changed,RS)
    if not changed:break
PR=[p for p in PR if p is not None]
log('tiles split by a river',NSPLIT,'slivers merged into a neighbour',NMERGE)
if '--tilesonly' in sys.argv:
    print('NF1',NF1,'NSKIP',NSKIP);sys.exit()
log('tiles',len(PR))
def _inside_len(pr):
    t=0
    for pt in polys(pr['poly']):
        cand=[RIVL[int(i)] for i in RT.query(pt)]
        if cand:t+=unary_union(cand).intersection(pt.buffer(-6e-4)).length
    return t
_in=sum(_inside_len(p) for p in PR);_tot=RIVU.length
import pickle
pickle.dump({'pr':[(p['r'],shapely.to_wkb(p['poly'])) for p in PR],'rv':{c:[shapely.to_wkb(l) for l in v] for c,v in RV.items()},'reg':REG},open('/tmp/dump.pkl','wb'))
log('RIVER CHECK: river length %.1f deg, still inside tile interiors %.2f deg (%.2f%%)'%(_tot,_in,100*_in/_tot))
# ---- terrain
reg=feats('ne_10m_geography_regions_polys')
mnt=unary_union([valid(shape(f['geometry'])) for f in reg if f['properties'].get('FEATURECLA')=='Range/mtn']).intersection(BOX)
plat=unary_union([valid(shape(f['geometry'])) for f in reg if f['properties'].get('FEATURECLA')=='Plateau']).intersection(BOX)
marsh=unary_union([box(26,51.2,31.8,52.9),box(3.5,51.0,7.3,53.6),box(28.5,44.3,30.0,45.5),box(21,53.3,24,54.2),box(30.0,29.5,32.2,31.6),box(31.0,6.5,32.5,9.5)])
desert=unary_union([box(-18.5,15,35,31),box(34,12,45,28),box(38,16,60,30),box(44,-5,51,11)])
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
    elif poly.intersection(desert).area>.5*a:t=6
    else:
        pf=.85 if (y>=56 and x>=10) else .4 if y>=56 else .55 if (50<=y<56 and x>=22) else .3 if (y>=46 and x>=18) else .14
        if y>=65:pf=.5
        if y<12 and x>8:pf=.55
        if y<37 and y>=12:pf=.02
        t=1 if h01(i)<pf else 0
    idx=ctree.query(poly,predicate='contains')
    if len(idx) and max(cities[j][3] for j in idx)>=1000000:t=5
    pr['t']=t
# ---- provinces: outline (union of their tiles, so it follows the tile edges exactly), label spot
tiles_of={}
for i,p in enumerate(PR): tiles_of.setdefault(p['r'],[]).append(i)
PV=[]
for ri,k in enumerate(REG):
    u=mp(unary_union([PR[i]['poly'] for i in tiles_of.get(ri,[])]))
    parts=sorted(polys(u),key=lambda p:-p.area);big=parts[0];rp=big.representative_point();cen=big.centroid
    lx,ly=(cen.x,cen.y) if big.contains(cen) else (rp.x,rp.y)
    b=big.bounds;size=min((b[2]-b[0])*.8/(len(k)*.78+.1),(b[3]-b[1])*.45,1.2)
    PV.append({'k':k,'u':u,'lx':lx,'ly':ly,'size':size})
log('province outlines')
pt=STRtree([p['poly'] for p in PR]);adj=[]
for i,p in enumerate(PR):
    for j in pt.query(p['poly'].buffer(.012)):
        j=int(j)
        if j<=i:continue
        q=PR[j]['poly']
        if p['poly'].distance(q)<.006 and p['poly'].buffer(.012).intersection(q).area>4e-5:adj+=[i,j]
log('adjacency pairs',len(adj)//2)
# tile pairs that touch ALONG a river: crossing between them is a river crossing
from shapely.prepared import prep as _prep
radj=[]
for q in range(0,len(adj),2):
    i,j=adj[q],adj[q+1];a=PR[i]['poly'];b=PR[j]['poly']
    sh=a.boundary.intersection(b.buffer(.002))
    if sh.is_empty or sh.length==0:continue
    cand=[RIVL[int(t)] for t in RT.query(sh.buffer(.004))]
    if not cand:continue
    if sh.intersection(unary_union(cand).buffer(.0035)).length>.5*sh.length:radj+=[i,j]
log('river-crossing pairs',len(radj)//2,'of',len(adj)//2)
# ---- encoding and outlines
CB=LAND.boundary.buffer(.007)
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
def lns(g,tol=.004,sub=False):
    out=[];b=(g.simplify(tol) if tol else g).boundary
    if sub:b=b.difference(CB)
    for l in ([b] if b.geom_type=='LineString' else list(getattr(b,'geoms',[]))):
        if l.geom_type=='LineString' and len(l.coords)>=2:
            e=enc(list(l.coords),False)
            if e:out.append(e)
    return out
provLines=[]
for v in PV:provLines+=lns(v['u'],.004,True)
import shapely.affinity
from PIL import Image, ImageDraw
# ---------------------------------------------------------------------------
# COUNTRY NAME LABELS
# One rule, as in HOI4: a country's colonies and exclaves carry the name of the country that holds them (Libya and the Italian
# East Africa read ITALY, Algeria / Tunisia / Morocco / French West Africa read FRANCE, Sudan and Kenya read UNITED KINGDOM ...).
# Only occupied Norway (1941) keeps its own historical name.  No modern or invented names, and no province names - ever.
#
# Placement (as in HOI4 the name FOLLOWS the land): a spine is run down the middle of the landmass in each of several directions,
# smoothed, and the lettering is laid along it - diagonal or gently curved - as large as it will go, with the letters spread out to fill the
# run.  The strip the letters occupy
#   * never covers land held by anybody else (water does NOT block it, so a name may overhang the coast or run on out over the sea), and
#   * keeps its middle on the landmass it names and at least MINOWN of its area on that country's own land.
# (A straight tilted box, fit_label, is kept as the fallback.)
# A country's capital landmass gets one name; every other big landmass (colonies, Libya, Sudan ...) gets a repeat of it, and a huge
# colonial block gets two.
# ---------------------------------------------------------------------------
BEBAS={'A':.401,'B':.404,'C':.383,'D':.406,'E':.363,'F':.344,'G':.391,'H':.42,'I':.192,'J':.265,'K':.414,'L':.344,'M':.538,'N':.427,'O':.4,'P':.386,'Q':.4,'R':.403,'S':.372,'T':.364,'U':.402,'V':.382,'W':.557,'X':.406,'Y':.394,'Z':.362,' ':.16,'-':.27,'.':.188,'&':.42}
LAB_FONT=1.3; LAB_SP=.36; LAB_CAP=.70*1.3      # game draws caps at size*1.3 px with size*.36 px between letters
def text_units(t): return sum(BEBAS.get(c,.4) for c in t)*LAB_FONT+LAB_SP*(len(t)-1)     # text width in units of "size"
LAB_PADX=.40; LAB_PADY=.34                                                                # breathing room (units of size) so lettering never touches a border
MINOWN=.75                                                                                # share of a name's box that must lie on the country's own land
LABEL_CITIES=[]    # (name,lon,lat,zoom it first appears at) - filled in by the builder

def _raster(g,x0,y0,res,nx,ny):
    im=Image.new('L',(nx,ny),0); d=ImageDraw.Draw(im)
    for p in polys(g):
        d.polygon([((x-x0)/res,(y-y0)/res) for x,y in p.exterior.coords],fill=1)
        for r in p.interiors: d.polygon([((x-x0)/res,(y-y0)/res) for x,y in r.coords],fill=0)
    return np.asarray(im,dtype=bool)

def label_rect(cx,cy,s,ang,text):
    w=(text_units(text)+LAB_PADX)*s/2; h=(LAB_CAP+LAB_PADY)*s/2
    return shapely.affinity.rotate(box(cx-w,cy-h,cx+w,cy+h),ang,origin=(cx,cy))

def fit_label(own,anchor,text,cap,blocks=(),smin=.07,angs=range(-55,56,5),rel=1.0):
    """Biggest text box that lies on own land or water (never on anyone else's land), centred on `anchor`, >=MINOWN on own land.
    Returns (size,cx,cy,angle_deg_ccw) or None."""
    W=text_units(text)+LAB_PADX; H=LAB_CAP+LAB_PADY
    ax0,ay0,ax1,ay1=anchor.bounds; ext=max(ax1-ax0,ay1-ay0)
    cap=min(cap,rel*ext/W)
    if cap<smin: return None
    pad=min(7.0,max(1.5,.35*ext))
    X0,Y0,X1,Y1=max(EXT[0],ax0-pad),max(EXT[1],ay0-pad),min(EXT[2],ax1+pad),min(EXT[3],ay1+pad)
    res=min(.07,max(.0035,math.sqrt(anchor.area)/80))
    while ((X1-X0)/res)*((Y1-Y0)/res)>160000: res*=1.2
    nx=int((X1-X0)/res)+1; ny=int((Y1-Y0)/res)+1
    roi=box(X0,Y0,X1,Y1)
    m_own=_raster(own.intersection(roi),X0,Y0,res,nx,ny)
    m_land=_raster(LAND.intersection(roi),X0,Y0,res,nx,ny)
    m_anc=_raster(anchor.intersection(roi),X0,Y0,res,nx,ny)
    d=m_own.copy()                                                  # land that is not ours, with a one-cell skirt round our own coast (the coastlines differ by a sliver)
    for dy in (-1,0,1):
        for dx in (-1,0,1): d|=np.roll(np.roll(m_own,dy,0),dx,1)
    m_ok=m_own|~(m_land&~d)
    for b in blocks: m_ok&=~_raster(b.intersection(roi),X0,Y0,res,nx,ny)
    cx0,cy0=(X0+X1)/2,(Y0+Y1)/2; hw,hh_=(X1-X0)/2,(Y1-Y0)/2
    best=None
    for a in angs:
        ca,sa=math.cos(math.radians(a)),math.sin(math.radians(a))
        ru=abs(ca)*hw+abs(sa)*hh_; rv=abs(sa)*hw+abs(ca)*hh_
        mx=int(2*ru/res)+1; my=int(2*rv/res)+1
        U,V=np.meshgrid(-ru+(np.arange(mx)+.5)*res,-rv+(np.arange(my)+.5)*res)
        ix=((cx0+U*ca-V*sa-X0)/res).astype(int); iy=((cy0+U*sa+V*ca-Y0)/res).astype(int)
        v=(ix>=0)&(ix<nx)&(iy>=0)&(iy<ny); ix=np.clip(ix,0,nx-1); iy=np.clip(iy,0,ny-1)
        ok=m_ok[iy,ix]&v; ow=m_own[iy,ix]&v; an=m_anc[iy,ix]&v
        So=np.zeros((my+1,mx+1),np.int32); So[1:,1:]=ok.cumsum(0).cumsum(1)
        Sw=np.zeros((my+1,mx+1),np.int32); Sw[1:,1:]=ow.cumsum(0).cumsum(1)
        def feas(s):
            ww=int(math.ceil(W*s/res)); hh=int(math.ceil(H*s/res))
            if ww>=mx or hh>=my: return None
            full=(So[hh:,ww:]-So[:-hh,ww:]-So[hh:,:-ww]+So[:-hh,:-ww])==ww*hh
            if not full.any(): return None
            own=(Sw[hh:,ww:]-Sw[:-hh,ww:]-Sw[hh:,:-ww]+Sw[:-hh,:-ww])>=MINOWN*ww*hh
            full&=own
            if not full.any(): return None
            pos=np.argwhere(full); pos=pos[an[pos[:,0]+hh//2,pos[:,1]+ww//2]]    # centre must sit on the landmass being named
            return (pos,ww,hh) if len(pos) else None
        if feas(smin) is None: continue
        lo,hi=smin,cap
        if feas(hi) is not None: lo=hi
        else:
            for _ in range(12):
                mid=(lo+hi)/2
                if feas(mid) is not None: lo=mid
                else: hi=mid
        pos,ww,hh=feas(lo)
        sc=lo*(1-.011*abs(a))
        if best is None or sc>best[0]:
            ccx,ccy=anchor.centroid.x-cx0,anchor.centroid.y-cy0          # of all the places the box fits, the one nearest the middle of the landmass
            uc=(ccx*ca+ccy*sa-(-ru))/res-ww/2; vc=(-ccx*sa+ccy*ca-(-rv))/res-hh/2
            j=int(np.argmin((pos[:,1]-uc)**2+(pos[:,0]-vc)**2)); iy_,ix_=pos[j]
            u_=-ru+(ix_+ww/2)*res; v_=-rv+(iy_+hh/2)*res
            best=(sc,lo*.97,cx0+u_*ca-v_*sa,cy0+u_*sa+v_*ca,a)
    return None if best is None else best[1:]


def fit_curve(own,anchor,text,cap,blocks=(),smin=.07,angs=range(-75,91,15)):
    """Lay `text` along a smoothed spine through `anchor`.  Returns (size, path[(lon,lat)...]) or None."""
    Wb=text_units(text); Hh=LAB_CAP+LAB_PADY
    ax0,ay0,ax1,ay1=anchor.bounds; ext=max(ax1-ax0,ay1-ay0)
    cap=min(cap,(.30 if anchor.area>=300 else .17)*math.sqrt(anchor.area))                       # lettering about a fifth of the territory's size, as HOI4 does
    if cap<smin: return None
    pad=min(7.0,max(1.5,.35*ext))
    X0,Y0,X1,Y1=max(EXT[0],ax0-pad),max(EXT[1],ay0-pad),min(EXT[2],ax1+pad),min(EXT[3],ay1+pad)
    res=min(.07,max(.0035,math.sqrt(anchor.area)/80))
    while ((X1-X0)/res)*((Y1-Y0)/res)>160000: res*=1.2
    nx=int((X1-X0)/res)+1; ny=int((Y1-Y0)/res)+1
    roi=box(X0,Y0,X1,Y1)
    m_own=_raster(own.intersection(roi),X0,Y0,res,nx,ny)
    m_land=_raster(LAND.intersection(roi),X0,Y0,res,nx,ny)
    m_anc=_raster(anchor.intersection(roi),X0,Y0,res,nx,ny)
    d=m_own.copy()
    for dy in (-1,0,1):
        for dx in (-1,0,1): d|=np.roll(np.roll(m_own,dy,0),dx,1)
    m_ok=m_own|~(m_land&~d)
    for b in blocks: m_ok&=~_raster(b.intersection(roi),X0,Y0,res,nx,ny)
    cx0,cy0=(X0+X1)/2,(Y0+Y1)/2; hw,hh_=(X1-X0)/2,(Y1-Y0)/2
    acx,acy=anchor.centroid.x-cx0,anchor.centroid.y-cy0
    best=None
    for a in angs:
        ca,sa=math.cos(math.radians(a)),math.sin(math.radians(a))
        ru=abs(ca)*hw+abs(sa)*hh_; rv=abs(sa)*hw+abs(ca)*hh_
        mx=int(2*ru/res)+1; my=int(2*rv/res)+1
        U,V=np.meshgrid(-ru+(np.arange(mx)+.5)*res,-rv+(np.arange(my)+.5)*res)
        ix=((cx0+U*ca-V*sa-X0)/res).astype(int); iy=((cy0+U*sa+V*ca-Y0)/res).astype(int)
        v=(ix>=0)&(ix<nx)&(iy>=0)&(iy<ny); ix=np.clip(ix,0,nx-1); iy=np.clip(iy,0,ny-1)
        ok=m_ok[iy,ix]&v; ow=m_own[iy,ix]&v; an=m_anc[iy,ix]&v
        # spine: the middle of the longest run of land in every column
        mid=np.full(mx,np.nan)
        for c in np.flatnonzero(an.any(0)):
            idx=np.flatnonzero(an[:,c]); br=np.flatnonzero(np.diff(idx)>1)
            st=np.r_[0,br+1]; en=np.r_[br,len(idx)-1]; k=int(np.argmax(en-st)); mid[c]=(idx[st[k]]+idx[en[k]])/2
        cv=np.flatnonzero(~np.isnan(mid))
        if len(cv)<4: continue
        c0,c1=cv[0],cv[-1]; n=c1-c0+1
        m=np.interp(np.arange(mx),cv,mid[cv])
        win=max(3,int(.22*n))|1; kern=np.ones(win)/win
        seg=np.pad(m[c0:c1+1],win//2,mode='edge'); m[c0:c1+1]=np.convolve(seg,kern,mode='valid')
        kk=max(3,int(.15*n))                                          # run the spine straight on past both coasts
        sl0=np.clip((m[c0+kk]-m[c0])/kk,-1.2,1.2); sl1=np.clip((m[c1]-m[c1-kk])/kk,-1.2,1.2)
        m[:c0]=m[c0]+sl0*(np.arange(c0)-c0); m[c1+1:]=m[c1]+sl1*(np.arange(c1+1,mx)-c1)
        m=np.clip(m,0,my-1)
        dm=np.gradient(m); g=np.sqrt(1+dm*dm); cum=np.concatenate([[0],np.cumsum(g*res)])   # cum[i]: arc length up to column i
        cso=np.zeros((my+1,mx),np.int32); cso[1:]=ok.cumsum(0); csw=np.zeros((my+1,mx),np.int32); csw[1:]=ow.cumsum(0)
        cols=np.arange(mx)
        lc=an[np.clip(np.round(m).astype(int),0,my-1),cols]              # is the spine over the landmass at this column?
        lcs=np.concatenate([[0],np.cumsum(lc)])
        def strips(s):
            half=(Hh*s/res)/2*g
            r0=np.floor(m-half).astype(int); r1=np.ceil(m+half).astype(int)
            inb=(r0>=0)&(r1<my); r0=np.clip(r0,0,my-1); r1=np.clip(r1,0,my-1)
            cnt=r1-r0+1
            full=inb&((cso[r1+1,cols]-cso[r0,cols])==cnt)
            fo=(csw[r1+1,cols]-csw[r0,cols])/cnt
            return full,fo
        def feas(s):
            Lb=Wb*s
            if Lb>=cum[-1]: return None
            full,fo=strips(s)
            i=np.arange(mx); j=np.searchsorted(cum,cum[i]+Lb)-1
            good=j<mx; i=i[good]; j=j[good]
            if not len(i): return None
            okc=np.concatenate([[0],np.cumsum(full)]); fc=np.concatenate([[0],np.cumsum(fo)])
            L=j-i+1
            good=((okc[j+1]-okc[i])==L)&(((fc[j+1]-fc[i])/L)>=MINOWN)&(((lcs[j+1]-lcs[i])/L)>=.85)
            i=i[good]; j=j[good]
            if not len(i): return None
            mc=(i+j)//2; good=an[np.clip(np.round(m[mc]).astype(int),0,my-1),mc]    # the middle of the name must sit on the landmass being named
            i=i[good]; j=j[good]
            return (i,j,full) if len(i) else None
        if feas(smin) is None: continue
        lo,hi=smin,min(cap,cum[-1]/Wb*.999)
        if feas(hi) is not None: lo=hi
        else:
            for _ in range(11):
                mid_=(lo+hi)/2
                if feas(mid_) is not None: lo=mid_
                else: hi=mid_
        i_,j_,full=feas(lo)
        cc=(acx*ca+acy*sa+ru)/res; q=int(np.argmin(np.abs((i_+j_)/2-cc)))
        i2,j2=int(i_[q]),int(j_[q]); Lb=Wb*lo
        while True:                                                    # spread the lettering out over the land it follows, up to 1.2x its natural length (tight lettering, as HOI4 has it)
            grew=False
            if cum[j2+1]-cum[i2]<1.2*Lb:
                if j2+1<mx-1 and full[j2+1] and lc[j2+1]: j2+=1; grew=True
                if i2>0 and full[i2-1] and lc[i2-1] and cum[j2+1]-cum[i2]<1.2*Lb: i2-=1; grew=True
            if not grew: break
        tot=cum[j2+1]-cum[i2]
        # biggest wins; at equal size the direction with the longest run (the land's long axis) wins, so the name follows the shape of the country
        sc=lo*(1+.5*min(1.0,(tot/Lb-1)/.2))*(1-.0012*abs(a))
        if best is None or sc>best[0]:
            arc=cum[i2:j2+2]-cum[i2]
            K=int(max(4,min(14,round(tot/(1.2*lo)))))
            tt=np.linspace(0,tot,K); cc_=np.interp(tt,arc,np.arange(i2,j2+2))
            mm=np.interp(cc_,np.arange(mx),m)
            Uc=-ru+(cc_+.5)*res; Vc=-rv+(mm+.5)*res
            px=cx0+Uc*ca-Vc*sa; py=cy0+Uc*sa+Vc*ca
            best=(sc,lo*.97,list(zip(px.tolist(),py.tolist())))
    return None if best is None else best[1:]

NE0=None
def ne_geom(codes):
    global NE0
    if NE0 is None:
        NE0={}
        for f in feats('ne_10m_admin_0_countries'):
            NE0.setdefault(f['properties'].get('ADM0_A3'),[]).append(valid(shape(f['geometry'])))
    return unary_union([g for c in codes for g in NE0.get(c,[])])
# the one place a state's land keeps its own historical name: Norway under German occupation (1 Jan 1941)
LABEL_OCCUPIED=[('NORWAY','GER',['NOR'])]
# the landmass that carries each state's main name (else its biggest landmass)
LABEL_HOME={'GER':(13.4,52.5),'SOV':(37.6,55.75),'FRA':(2.35,47.0),'VIC':(2.35,46.5),'GBR':(-1.5,52.8),'ITA':(12.5,42.5),'ESP':(-3.7,40.4),'POR':(-8.0,39.5),'BEL':(4.6,50.6)}

def _footprints(s):
    return [box(cx-10/max(zc,9/s),cy-9/max(zc,9/s),cx+(16+5.6*len(cn))/max(zc,9/s),cy+17/max(zc,9/s)) for (cn,cx,cy,zc) in LABEL_CITIES]

def _place(own,comps,home,name,cap,gi,out):
    placed=[]
    for comp in sorted(comps,key=lambda p:(p is not home,-p.area)):
        is_home=comp is home
        if not is_home and comp.area<3: continue                               # islands and tiny exclaves do not repeat the name
        if not is_home and any(r.intersects(comp) for r in placed): continue   # already lying under a name
        reps=1                                                                  # one name per landmass, however big (HOI4 draws one large FRANCE over French Africa)
        blocks=[]; s1=None
        for k in range(reps):
            cp=cap if is_home else min(cap,4.5 if comp.area>=300 else 2.0)
            r=fit_curve(own,comp,name,cp,blocks)
            if r is None: break
            s,path=r
            if not is_home and s<.5: break
            if k>0 and s<.65*s1: break
            band=LineString(path).buffer((LAB_CAP+LAB_PADY)*s/2)
            fp=unary_union(_footprints(s)) if LABEL_CITIES else None      # keep clear of city names where that costs little
            if fp is not None and band.intersects(fp):
                r2=fit_curve(own,comp,name,cp,blocks+[fp])
                if r2 and r2[0]>=.75*s:
                    m1=LineString(path).interpolate(.5,normalized=True); m2=LineString(r2[1]).interpolate(.5,normalized=True)
                    if m1.distance(m2)<=.2*text_units(name)*s: s,path=r2; band=LineString(path).buffer((LAB_CAP+LAB_PADY)*s/2)
            s1=s1 or s; placed.append(band); blocks.append(band.buffer(1.0*s))
            mp_=LineString(path).interpolate(.5,normalized=True)
            out.append([gi,round(mp_.x,2),round(mp_.y,2),0,round(float(s),2),name,[round(float(c),2) for xy in path for c in xy]])

def make_labels(gi,code,gname,u,cap=3.2):
    """u: everything group gi holds.  Returns [group, lon, lat, angle, size, NAME] label records."""
    out=[]; occ=None
    for nm,oc,codes in LABEL_OCCUPIED:
        if oc==code:
            zg=ne_geom(codes); reg=mp(u.intersection(zg))
            if not reg.is_empty:
                pcs=[p for p in polys(reg) if p.area>=3]
                if pcs: _place(reg,pcs,max(pcs,key=lambda p:p.area),nm,2.4,gi,out)
                u=mp(u.difference(zg))
    comps=[p for p in polys(u) if p.area>=.05]
    if not comps: return out
    hp=LABEL_HOME.get(code); home=None
    if hp:
        for p in comps:
            if p.buffer(.05).contains(Point(*hp)): home=p;break
    if home is None: home=max(comps,key=lambda p:p.area)
    _place(u,comps,home,gname.upper(),cap,gi,out)
    return out
LABEL_CITIES[:]=[(c[0],c[1],c[2],26 if c[4] else 22 if c[3]>=3000000 else 34 if c[3]>=1000000 else 60) for c in cities]
SCN=[('1939','1 August 1939',1939,8,1),('1941','1 January 1941',1941,1,1)]
scn=[]
for si,(sid,date,y,m,d) in enumerate(SCN):
    own=[GI[OWNP[k][si]] for k in REG];GU={}
    for ri,k in enumerate(REG): GU.setdefault(own[ri],[]).append(PV[ri]['u'])
    lines=[];labs=[]
    for g,lst in GU.items():
        u=mp(unary_union(lst).buffer(.004).buffer(-.004));lines+=lns(u,.004,True)
        labs+=make_labels(g,GROUPS[g][0],GROUPS[g][1],u,3.3 if GROUPS[g][0]=='SOV' else 3.2)
    kept=sorted(labs,key=lambda l:-l[4])
    scn.append({'id':sid,'date':date,'own':own,'lines':lines,'labs':kept})
    log('scenario',sid,len(GU),'owners')
LANDR=rings(LAND,.006)
lakes=[]
for f in feats('ne_10m_lakes'):
    for p in polys(valid(shape(f['geometry'])).intersection(BOX)):
        if p.area>.12:lakes+=rings(p,.01)
rivers=[]      # [class1 lines, class2 lines, class3 lines]: the same lines the tiles were cut along
for c in (1,2,3):
    g=unary_union(RV[c]).intersection(LAND);lst=[]
    for l in lines_of(linemerge(g) if g.geom_type!='LineString' else g):
        e=enc(list(l.coords),False)
        if e:lst.append(e)
    rivers.append(lst)
# ---- country ranks and victory points
# Six ranks (the last is the Micro State). Every country gets a rank in each scenario (RANK_BY_SCN, default 4 = Minor Power).
# Victory points follow the LAND (like HOI4): each 1939 country has a VP budget set by its 1939 rank, split 30% capital,
# 50% biggest cities, 20% smaller cities, shared out by population (Natural Earth: mean of city-proper and metro figures, modern values: no 1939 census is available).
RANKS=[('Hyperpower','Unchallenged global dominance',150),('Superpower','Can project power across the entire world and compete for global dominance',100),
 ('Great Power','One of the world\'s strongest countries, with major military and economic influence',70),('Regional Power','Dominates or heavily influences its geographical region',40),
 ('Minor Power','Has a functioning military and economy but limited international influence',20),
 ('Micro State','A tiny sovereign state or free city, too small to field a real army',5)]
MICRO={'DAN':5,'VAT':5,'MON':5,'SMR':5,'AND':5,'LIE':5}    # micro states: capital-only VP, no army to speak of
RANK_BY_SCN=[{**MICRO,'GBR':1,'GER':2,'SOV':2,'FRA':2,'ITA':2,'POL':3,'TUR':3,'YUG':3},   # 1 Aug 1939 (nobody on this map is a Hyperpower)
             {**MICRO,'GBR':1,'GER':2,'SOV':2,'FRA':2,'ITA':2,'POL':3,'TUR':3,'YUG':3}]   # 1 Jan 1941 (VIC and FFR fall back to Minor Power)
BIGN=[8,6,5,4,3,0]; SMALLN=[30,20,14,8,5,0]       # how many 'big' and 'smaller' cities each rank can have
CAPS={'GER':'Berlin','SOV':'Moscow','POL':'Warsaw','FRA':'Paris','GBR':'London','ITA':'Rome','ESP':'Madrid','POR':'Lisbon','SWE':'Stockholm','FIN':'Helsinki','EST':'Tallinn','LAT':'Riga','LIT':'Kaunas','DEN':'Kobenhavn',
 'NOR':'Oslo','NED':'Amsterdam','BEL':'Brussels','LUX':'Luxembourg','SWI':'Bern','YUG':'Belgrade','GRE':'Athens','ROM':'Bucharest','HUN':'Budapest','SVK':'Bratislava','BUL':'Sofia','TUR':'Ankara','IRL':'Dublin',
'IRQ':'Baghdad','IRN':'Tehran','SAU':'Riyadh','YEM':'Sanaa','LBR':'Monrovia',
 'DAN':'Gdansk','VAT':'Vatican City','MON':'Monaco','SMR':'San Marino','AND':'Andorra la Vella','LIE':'Vaduz'}
METRO={'GBR':{'United Kingdom','Isle of Man','Jersey','Guernsey'},'FRA':{'France','Monaco'},'ITA':{'Italy','San Marino','Vatican'},'ESP':{'Spain','Andorra'},'POR':{'Portugal'},'BEL':{'Belgium'}}  # colonies only fill 'smaller city' slots
ptree=STRtree([v['u'] for v in PV])
def prov_of(x,y):
    pt=Point(x,y);idx=ptree.query(pt,predicate='intersects')
    if len(idx):return int(idx[0])
    j=int(ptree.nearest(pt));return j if PV[j]['u'].distance(pt)<.35 else -1
places={}
for f in feats('ne_10m_populated_places_simple'):
    p=f['properties'];x,y=f['geometry']['coordinates'][:2]
    if not(EXT[0]<x<EXT[2] and EXT[1]<y<EXT[3]) or (p['pop_max'] or 0)<50000:continue
    pi=prov_of(x,y)
    if pi<0:continue
    places.setdefault(OWNP[REG[pi]][0],[]).append({'n':p['name'],'a':p['nameascii'],'x':x,'y':y,'pop':((p['pop_min'] or 0)+(p['pop_max'] or 0))/2,'adm':p['adm0name'],'pi':pi})
# micro states are too small for any 50,000+ city to fall inside them, so their capital is placed by hand (capital-only VP)
MICRO_CAPS={'DAN':('Gdańsk',18.65,54.35),'VAT':('Vatican City',12.45,41.9),'MON':('Monaco',7.42,43.74),'SMR':('San Marino',12.45,43.94),'AND':('Andorra la Vella',1.52,42.51),'LIE':('Vaduz',9.52,47.14)}
for _c,(_n,_x,_y) in MICRO_CAPS.items():
    _pi=prov_of(_x,_y)
    if _pi>=0 and OWNP[REG[_pi]][0]==_c: places[_c]=[{'n':_n,'a':_n,'x':_x,'y':_y,'pop':50000,'adm':'micro','pi':_pi}]
    else: log('WARNING micro state capital not inside its province:',_c,_n)
def share(total,items):
    n=len(items);w=[max(c['pop'],1) for c in items];S=sum(w);rest=total-n;raw=[rest*wi/S for wi in w];fl=[int(r) for r in raw];left=rest-sum(fl)
    for i in sorted(range(n),key=lambda i:-(raw[i]-fl[i]))[:left]:fl[i]+=1
    return [1+v for v in fl]
VPL=[]
for gi,(code,gname,_) in enumerate(GROUPS):
    cs=places.get(code)
    if not cs:continue
    r=RANK_BY_SCN[0].get(code,4);B=RANKS[r][2];metro=METRO.get(code);ism=lambda c:metro is None or c['adm'] in metro
    cap=next((c for c in cs if c['a']==CAPS.get(code) and ism(c)),None) or max([c for c in cs if ism(c)] or cs,key=lambda c:c['pop'])
    capvp=round(.30*B);pb=round(.50*B);ps=B-capvp-pb
    rest=sorted([c for c in cs if c is not cap],key=lambda c:-c['pop'])
    bigs=[c for c in rest if ism(c)][:min(BIGN[r],pb)];rest=[c for c in rest if c not in bigs]
    smalls=[];ncol=0;maxcol=SMALLN[r]//3 if metro else SMALLN[r]
    for c in rest:
        if len(smalls)>=min(SMALLN[r],ps):break
        if not ism(c):
            if ncol>=maxcol:continue
            ncol+=1
        smalls.append(c)
    if not smalls:pb+=ps;ps=0
    if not bigs:capvp+=pb;pb=0
    ent=[(cap,capvp,2)]+[(c,v,1) for c,v in zip(bigs,share(pb,bigs) if bigs else [])]+[(c,v,0) for c,v in zip(smalls,share(ps,smalls) if smalls else [])]
    assert sum(v for _,v,_ in ent)==B,(code,sum(v for _,v,_ in ent),B)
    for c,v,k in ent:VPL.append([c['n'],round(c['x'],2),round(c['y'],2),v,k,c['pi'],gi])
    log('VP',code,RANKS[r][0],B,'capital',cap['n'],capvp,'| big',len(bigs),'| small',len(smalls))
for si in range(len(scn)):
    present=set(scn[si]['own']);scn[si]['rank']=[(RANK_BY_SCN[si].get(g[0],4) if gi in present else -1) for gi,g in enumerate(GROUPS)]
out={'name':'Invicta Europa','ext':list(EXT),'q':Q,'groups':[{'k':g[0],'n':g[1],'c':g[2]} for g in GROUPS],
 'provs':[[v['k'],round(v['lx'],2),round(v['ly'],2),round(v['size'],2),rings(v['u'])] for v in PV],
 'tiles':[[p['r'],p['t'],round(p['c'][0],2),round(p['c'][1],2),rings(p['poly'])] for p in PR],
 'ranks':[{'n':a,'d':b,'vp':c} for a,b,c in RANKS],'vp':VPL,'scn':scn,'provLines':provLines,'coast':lns(LAND,.006),'land':LANDR,'lakes':lakes,'rivers':rivers,'radj':radj,
 'radj':radj,'cities':[[c[0],round(c[1],2),round(c[2],2),(2 if c[3]>=3000000 else 1 if c[3]>=1000000 else 0),1 if c[4] else 0] for c in cities],'adj':adj}
js='window.MAP='+json.dumps(out,separators=(',',':'),ensure_ascii=False)+';\n'
open(os.path.join(HERE,'map-data.js'),'w',encoding='utf-8').write(js)
log('wrote map-data.js',round(len(js)/1e6,2),'MB; tiles',len(PR),'provinces',len(PV),'cities',len(cities))
if '--preview' in sys.argv:
    from PIL import Image, ImageDraw
    S=9;W=int((EXT[2]-EXT[0])*S);H=int((EXT[3]-EXT[1])*S);im=Image.new('RGB',(W*2,H),(43,74,111));d=ImageDraw.Draw(im)
    for si in (0,1):
        ox=si*W;tx=lambda x:(x-EXT[0])*S+ox;ty=lambda y:(EXT[3]-y)*S
        for p in polys(LAND):d.polygon([(tx(x),ty(y)) for x,y in p.exterior.coords],fill=(120,118,104))
        for p in PR:
            o=scn[si]['own'][p['r']]
            c=GROUPS[o][2];rgb=tuple(int(c[i:i+2],16) for i in (1,3,5));col=tuple(min(255,int(v*(.94+.1*h01(p['r'])))) for v in rgb)
            for q in polys(p['poly']):d.polygon([(tx(x),ty(y)) for x,y in q.exterior.coords],fill=col,outline=tuple(int(v*.75) for v in col))
    im.save('/tmp/preview.png');log('preview saved')
