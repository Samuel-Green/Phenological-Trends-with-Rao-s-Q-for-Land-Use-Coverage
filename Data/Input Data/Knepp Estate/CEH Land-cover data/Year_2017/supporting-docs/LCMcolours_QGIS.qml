<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>
<qgis maxScale="0" styleCategories="AllStyleCategories" hasScaleBasedVisibilityFlag="0" version="3.18.1-Zürich" minScale="1e+08">
  <flags>
    <Identifiable>1</Identifiable>
    <Removable>1</Removable>
    <Searchable>1</Searchable>
    <Private>0</Private>
  </flags>
  <temporal enabled="0" fetchMode="0" mode="0">
    <fixedRange>
      <start></start>
      <end></end>
    </fixedRange>
  </temporal>
  <customproperties>
    <property key="WMSBackgroundLayer" value="false"/>
    <property key="WMSPublishDataSourceUrl" value="false"/>
    <property key="embeddedWidgets/count" value="0"/>
    <property key="identify/format" value="Value"/>
  </customproperties>
  <pipe>
    <provider>
      <resampling zoomedOutResamplingMethod="nearestNeighbour" enabled="false" zoomedInResamplingMethod="nearestNeighbour" maxOversampling="2"/>
    </provider>
    <rasterrenderer band="1" nodataColor="" alphaBand="-1" type="paletted" opacity="1">
      <rasterTransparency/>
      <minMaxOrigin>
        <limits>None</limits>
        <extent>WholeRaster</extent>
        <statAccuracy>Estimated</statAccuracy>
        <cumulativeCutLower>0.02</cumulativeCutLower>
        <cumulativeCutUpper>0.98</cumulativeCutUpper>
        <stdDevFactor>2</stdDevFactor>
      </minMaxOrigin>
      <colorPalette>
        <paletteEntry color="#e10000" alpha="255" value="1" label="Broadleaved woodland"/>
        <paletteEntry color="#006600" alpha="255" value="2" label="Coniferous woodland"/>
        <paletteEntry color="#732600" alpha="255" value="3" label="Arable and horticulture"/>
        <paletteEntry color="#00ff00" alpha="255" value="4" label="Improved grassland"/>
        <paletteEntry color="#7fe57f" alpha="255" value="5" label="Neutral grassland"/>
        <paletteEntry color="#70a800" alpha="255" value="6" label="Calcareous grassland"/>
        <paletteEntry color="#998100" alpha="255" value="7" label="Acid grassland"/>
        <paletteEntry color="#ffff00" alpha="255" value="8" label="Fen, marsh and swamp"/>
        <paletteEntry color="#801a80" alpha="255" value="9" label="Heather"/>
        <paletteEntry color="#e68ca6" alpha="255" value="10" label="Heather grassland"/>
        <paletteEntry color="#008073" alpha="255" value="11" label="Bog"/>
        <paletteEntry color="#d2d2ff" alpha="255" value="12" label="Inland rock"/>
        <paletteEntry color="#000080" alpha="255" value="13" label="Saltwater"/>
        <paletteEntry color="#0000ff" alpha="255" value="14" label="Freshwater"/>
        <paletteEntry color="#ccaa00" alpha="255" value="15" label="Supralittoral rock"/>
        <paletteEntry color="#ccb300" alpha="255" value="16" label="Supralittoral sediment"/>
        <paletteEntry color="#ffff80" alpha="255" value="17" label="Littoral rock"/>
        <paletteEntry color="#ffff80" alpha="255" value="18" label="Littoral sediment"/>
        <paletteEntry color="#8080ff" alpha="255" value="19" label="Saltmarsh"/>
        <paletteEntry color="#000000" alpha="255" value="20" label="Urban"/>
        <paletteEntry color="#808080" alpha="255" value="21" label="Suburban"/>
      </colorPalette>
      <colorramp name="[source]" type="randomcolors">
        <Option/>
      </colorramp>
    </rasterrenderer>
    <brightnesscontrast contrast="0" brightness="0" gamma="1"/>
    <huesaturation colorizeRed="255" colorizeStrength="100" saturation="0" colorizeOn="0" colorizeBlue="128" grayscaleMode="0" colorizeGreen="128"/>
    <rasterresampler maxOversampling="2"/>
    <resamplingStage>resamplingFilter</resamplingStage>
  </pipe>
  <blendMode>0</blendMode>
</qgis>
