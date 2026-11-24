# Invicta Europa

Browser grand strategy game. Europe and North Africa, starting 1 August 1939 (or 1 January 1941). Runs entirely in the browser, so it can be hosted free on GitHub Pages.

Open `index.html` (or the GitHub Pages link) to play. No install or build step is needed.

## Controls
- **Left click**: select. Click a tile to outline its whole province in white. Hold and drag to draw a selection box (for units, once units exist).
- **Middle mouse drag**: pan the map.
- **Mouse wheel**: zoom. **WASD / arrow keys** also pan. **Esc** clears the selection.

## Choosing a nation
- **Main menu:** the two start dates sit right of the title and the map sits under it. Clicking a start date switches the map to that scenario. Hovering one opens a panel of hanging flags for that scenario's great powers (rank Great Power or higher), strung on one thin thread and grouped by faction (Axis, Allies, Comintern). Click a banner to choose that nation (hovering one also lights it up on the menu map); double-click a banner to begin straight away. **Clear** removes the choice.
- **No flag chosen:** Begin Campaign opens the full map in "choose your nation" mode. Hover to light up a whole country, left-click it, then press **Play** (or Enter). Esc cancels.
- Factions are in `FACTION_OF`, and the rank cut-off for "great powers" is `MAJOR_MAX_RANK`, both in `index.html`. Flags are drawn in code: `flagSVG` for the small flat flags (HUD, confirm bar) and `flatArt` + `silkURL` for the hanging flags (the flat flag is painted on a canvas, then lit like folded silk and cached); countries without a drawn flag show their colour and code. Germany uses the 1939 national flag.

## Map
- A **province** is a historical first-level division of 1939 (Prussian provinces, 1939 Polish voivodeships, French departements, Italian provinces, Spanish provinces and so on).
- Counties and departements are then **merged into provinces about the size of Munster** (~24,600 km2): neighbours of the same country and owner are joined, big provinces (Bavaria, Russian oblasts) stay as they are, and islands, enclaves and city-states (Malta, Gibraltar, Berlin) stay small. Ireland is Munster, Leinster and Connacht.
- The **Sahara** (the area marked by hand, treated as unplayable desert) is one big province per country: French Sahara, Eastern Sahara and Libyan Desert. French Sahara also takes Chad's desert and Mauritania's Adrar and Tiris Zemmour; because a province has one owner, Chad's desert is Vichy in 1941 (the rest of Chad stays Free French). The northern half of Constantine (the coast) is part of Alger, so the Tunisia-Morocco coast road stays open. Tiles, coasts and country borders are unchanged. 710 provinces in total.
- A **tile** is a small cell inside a province. Tiles never cross a province border.
- Geography comes from [Natural Earth](https://www.naturalearthdata.com/) (public domain). Province borders are approximations of 1939 borders built from modern admin-1 data plus hand-drawn cuts.

## Rebuilding the map
`build_map.py` regenerates `map-data.js`:

```
pip install shapely numpy pillow
python3 build_map.py            # add --preview for a PNG, bigger number = bigger tiles
```

## Country ranks and Victory Points

Six ranks, each with a Victory Point (VP) budget. A nation's VP are spread over its cities in a 30 / 50 / 20 split
(capital / biggest cities / smaller cities), weighted by city population.

| Rank | Description | VP |
|---|---|---|
| Hyperpower | Unchallenged global dominance | 150 |
| Superpower | Can project power across the entire world and compete for global dominance | 100 |
| Great Power | One of the world's strongest countries, with major military and economic influence | 70 |
| Regional Power | Dominates or heavily influences its geographical region | 40 |
| Minor Power | Has a functioning military and economy but limited international influence | 20 |
| Micro State | A tiny sovereign state or free city, too small to field a real army | 5 (all in the capital) |

Ranks in both scenarios: Superpower = Britain; Great Power = Germany, Soviet Union, France, Italy;
Regional Power = Poland, Turkey, Yugoslavia; the micro states (Free City of Danzig, Vatican City, Monaco, San Marino, Andorra, Liechtenstein) are Micro States; everyone else is a Minor Power. Danzig is its own country on 1 August 1939 and part of Germany on 1 January 1941 (annexed 1 September 1939). Nobody on this map is a Hyperpower.
Change them in `RANK_BY_SCN` in `build_map.py` (one table per scenario), and the VP budgets in `RANKS`.

VP follow the land, as in HOI4: if a province is conquered, its VP go with it. City populations come from Natural Earth
(modern figures, as no 1939 census is available), so a few cities rank differently from 1939.
Colonial cities only fill small-city slots, so empires don't get all their VP from overseas.

## Names on the map

Only country names and cities are drawn on the map. Province names never appear (a province's name is in the hover tooltip and the info panel), and only
historical names are used: a country's colonies and exclaves carry the name of the state that holds them, as in HOI4 (Libya and East Africa read ITALY,
Algeria, Tunisia, Morocco and French West Africa read FRANCE, Sudan and Egypt read UNITED KINGDOM). Occupied Norway (1941) keeps its own name.
Names follow the land the way HOI4's do: each is fitted by `build_map.py` along ONE gentle curve (an arc or a straight line, never an S-bend) through the middle of its
territory, using the real Bebas Neue letter widths, so ITALY runs down the boot and UNITED KINGDOM sits diagonally on Britain. A name takes at most ~60% of the land along its
line, so about a fifth of the territory is left clear at each end, and it is scaled down whenever it would not fit. It never covers another country's land. The curve is stored
as an optional 7th element of each label record (`[group,lon,lat,angle,size,name,[lon,lat,...]]`) and drawn letter by letter in `index.html` (`pathText`). It also steers clear of city names where that costs little.
The capital landmass gets one name and each other sizeable landmass one more (a colonial block like French Africa gets a single large name, never two); names on tiny territories (Luxembourg) are tiny and only show once you zoom in.

Egypt is part of the United Kingdom's territory, and Italy holds all of East Africa (British Somaliland from August 1940), so there are fewer nations to draw.

Fonts: Bebas Neue (country names), Oswald (headings and buttons), Libre Baskerville (dates and historical text),
IBM Plex Mono (statistics). All are self-hosted in `fonts/`.

## Rivers are tile borders

No tile ever straddles a river: rivers are the borders between tiles, so every river crossing is a crossing between two tiles.

How it works (`build_map.py`):
1. **River set.** Natural Earth's river lines plus its Europe supplement (named rivers only), in three classes (1 great, 2 major, 3 minor) that set the drawn thickness. South of 37°N only the big rivers are used, because the desert is full of dry wadis. Duplicate branches are removed. The Aisne, Sambre, Lys, Yser, Dyle, Albert Canal, Kiel Canal, Isonzo, Tagliamento, Tisza, Váh, Vardar, Sozh, Orontes, Medjerda, Mersey and Clyde, plus parts of the Neman and Saône, are missing from Natural Earth, so they are hand-drawn from waypoints in `HAND` (good to a few km, not survey-grade). A hand river only adds the stretches Natural Earth lacks.
2. **Borders follow rivers.** Province borders that run within about 3.5 km of a river are moved onto it.
3. **Banks.** Each province is cut along its rivers into banks. A river that ends inside a province (a source) is continued in a straight line to the nearest border or river, at most ~3.5 tile radii; the continuation is an ordinary tile border, not a river.
4. **Tiles are made inside one bank.** A second pass splits any tile that still has river in its interior.

Result: about 97.7% of all river length lies exactly on tile edges (98.9% within 3 km); great rivers 97.9%. The rest is short stretches, mostly where a source ends inside a very large tile.

Data: `rivers` = three lists of lines (class 1, 2, 3). `radj` = pairs of tile indices that touch **along a river** (flat list `[i,j,i,j…]`, the river-crossing pairs; `adj` stays the plain neighbour list). Selecting a tile shows how many of its neighbours are across a river.

Debugging the cutting for one province: `python3 build_map.py --debug=Aisne` or `--debugxy=4.4,51.0`.

## Performance

The map has 13,429 tiles inside 731 provinces. Border line work is not stroked per frame: zoomed out (zoom <= 30) it is baked into transparent rasters (`ovl`, pre-built in the background by `prebake`), and zoomed in it is drawn as spatially culled chunks (`chunkLines` / `strokeCh`), so only borders on screen cost anything. The geometry is the same as in `map-data.js`; only the way it is drawn changed.
