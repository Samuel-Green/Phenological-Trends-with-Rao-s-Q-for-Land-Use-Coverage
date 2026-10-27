<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>
<qgis styleCategories="AllStyleCategories" maxScale="0" minScale="1e+08" hasScaleBasedVisibilityFlag="0" version="3.4.11-Madeira">
  <flags>
    <Identifiable>1</Identifiable>
    <Removable>1</Removable>
    <Searchable>1</Searchable>
  </flags>
  <customproperties>
    <property value="false" key="WMSBackgroundLayer"/>
    <property value="false" key="WMSPublishDataSourceUrl"/>
    <property value="0" key="embeddedWidgets/count"/>
    <property value="Value" key="identify/format"/>
  </customproperties>
  <pipe>
    <rasterrenderer type="paletted" opacity="1" alphaBand="-1" band="1">
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
        <paletteEntry value="1" alpha="255" color="#ff0000" label="Broadleaved woodland"/>
        <paletteEntry value="2" alpha="255" color="#006600" label="Coniferous woodland"/>
        <paletteEntry value="3" alpha="255" color="#732600" label="Arable"/>
        <paletteEntry value="4" alpha="255" color="#00ff00" label="Improved grassland"/>
        <paletteEntry value="5" alpha="255" color="#7fe57f" label="Neutral grassland"/>
        <paletteEntry value="6" alpha="255" color="#70a800" label="Calcareous grassland"/>
        <paletteEntry value="7" alpha="255" color="#998100" label="Acid grassland"/>
        <paletteEntry value="8" alpha="255" color="#ffff00" label="Fen, Marsh and Swamp"/>
        <paletteEntry value="9" alpha="255" color="#801a80" label="Heather and shrub"/>
        <paletteEntry value="10" alpha="255" color="#e68ca6" label="Heather grassland"/>
        <paletteEntry value="11" alpha="255" color="#008073" label="Bog"/>
        <paletteEntry value="12" alpha="255" color="#d2d2ff" label="Inland rock"/>
        <paletteEntry value="13" alpha="255" color="#000080" label="Saltwater"/>
        <paletteEntry value="14" alpha="255" color="#0000ff" label="Freshwater"/>
        <paletteEntry value="15" alpha="255" color="#ccb300" label="Supralittoral rock"/>
        <paletteEntry value="16" alpha="255" color="#ccb300" label="Supralittoral sediment"/>
        <paletteEntry value="17" alpha="255" color="#ffff80" label="Littoral rock"/>
        <paletteEntry value="18" alpha="255" color="#ffff80" label="Littoral sediment"/>
        <paletteEntry value="19" alpha="255" color="#8080ff" label="Saltmarsh"/>
        <paletteEntry value="20" alpha="255" color="#000000" label="Urban"/>
        <paletteEntry value="21" alpha="255" color="#808080" label="Suburban"/>
      </colorPalette>
      <colorramp type="randomcolors" name="[source]"/>
    </rasterrenderer>
    <brightnesscontrast brightness="0" contrast="0"/>
    <huesaturation colorizeGreen="128" colorizeBlue="128" colorizeOn="0" colorizeRed="255" colorizeStrength="100" grayscaleMode="0" saturation="0"/>
    <rasterresampler maxOversampling="2"/>
  </pipe>
  <blendMode>0</blendMode>
</qgis>
