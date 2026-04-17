import ee

# -------------------------------
# Authenticate & initialise
# -------------------------------
ee.Authenticate()   # Only needed first time
ee.Initialize(project='twdtw-paper-project')

# -------------------------------
# LOAD YOUR AOI (VERY IMPORTANT)
# -------------------------------
# Replace this with your asset path
knepp = ee.FeatureCollection("projects/twdtw-paper-project/assets/Knepp_Estate/knepp_boundary_lcc_pretty")

# -------------------------------
# Step 0 — Buffer AOI
# -------------------------------
buffer_distance = 5000  # metres
knepp_buffered = knepp.geometry().buffer(buffer_distance)

# -------------------------------
# Step 1 — Load Landsat collections
# -------------------------------
L5 = ee.ImageCollection("LANDSAT/LT05/C02/T1_L2")
L7 = ee.ImageCollection("LANDSAT/LE07/C02/T1_L2")
L8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
L9 = ee.ImageCollection("LANDSAT/LC09/C02/T1_L2")

# -------------------------------
# Step 2 — Preprocess
# -------------------------------
def preprocess(image):
    optical = image.select("SR_B.*").multiply(0.0000275).add(-0.2)

    qa = image.select("QA_PIXEL")
    mask = qa.bitwiseAnd(1 << 3).eq(0).And(
           qa.bitwiseAnd(1 << 4).eq(0))

    return image.addBands(optical, None, True).updateMask(mask)

# -------------------------------
# Step 3 — Rename bands
# -------------------------------
def rename_L57(image):
    return image.select(
        ["SR_B1", "SR_B2", "SR_B3", "SR_B4"],
        ["Blue", "Green", "Red", "NIR"]
    )

def rename_L89(image):
    return image.select(
        ["SR_B2", "SR_B3", "SR_B4", "SR_B5"],
        ["Blue", "Green", "Red", "NIR"]
    )

# -------------------------------
# Step 4 — Merge collections
# -------------------------------
landsat = (
    L5.map(preprocess).map(rename_L57)
    .merge(L7.map(preprocess).map(rename_L57))
    .merge(L8.map(preprocess).map(rename_L89))
    .merge(L9.map(preprocess).map(rename_L89))
    .filterBounds(knepp_buffered)
    .filterDate("2000-01-01", "2025-12-31")
)

# -------------------------------
# Step 5 — Add NDVI
# -------------------------------
def add_ndvi(image):
    ndvi = image.normalizedDifference(["NIR", "Red"]).rename("NDVI")
    return image.addBands(ndvi).toFloat()

with_indices = landsat.map(add_ndvi)

# -------------------------------
# Step 6 — Monthly composites (LIGHT)
# -------------------------------
def make_month(m):
    start = ee.Date("2000-01-01").advance(m, "month")
    end = start.advance(1, "month")

    filtered = with_indices.filterDate(start, end)

    composite = filtered.median()

    return ee.Algorithms.If(
        filtered.size().gt(0),
        composite
        .set("system:time_start", start.millis())
        .clip(knepp_buffered)
        .toFloat(),
        None
    )

months = ee.List.sequence(0, 12 * (2025 - 2000))

monthly = ee.ImageCollection.fromImages(
    months.map(make_month)
).filter(ee.Filter.notNull(["system:time_start"]))

# Convert to list for iteration
img_list = monthly.toList(monthly.size())

# -------------------------------
# Step 7 — EXPORT LOOP
# -------------------------------
n = img_list.size().getInfo()

tasks = []

for i in range(n):
    img = ee.Image(img_list.get(i))

    date = ee.Date(img.get("system:time_start")).format("YYYY_MM").getInfo()

    start = ee.Date(img.get("system:time_start"))
    end = start.advance(1, "month")

    filtered = with_indices.filterDate(start, end)

    # ---- Solar angles ----
    mean_sun_el = ee.Number(filtered.aggregate_mean("SUN_ELEVATION"))
    mean_sun_az = ee.Number(filtered.aggregate_mean("SUN_AZIMUTH"))

    sza = ee.Number(90).subtract(mean_sun_el)

    sza_img = ee.Image.constant(sza).rename("SZA")
    saa_img = ee.Image.constant(mean_sun_az).rename("SAA")

    final_img = img.addBands([sza_img, saa_img]).toFloat()

    task = ee.batch.Export.image.toDrive(
        image=final_img,
        description=f"Knepp_Estate_with_Buffer_{date}",
        folder="Knepp_Estate_Buffered",
        region=knepp_buffered,
        scale=30,
        maxPixels=1e13
    )

    task.start()
    tasks.append(task)

print(f"Submitted {len(tasks)} tasks!")
