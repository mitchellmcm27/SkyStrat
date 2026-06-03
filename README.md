# SkyStrat

SkyStrat is a GIS platform for mapping stratigraphic height by fitting cylindrical fold geometry to structural measurements.
Structural measurements (strikes/dips) can be measured in the field or calculated by SkyStrat solving a three-point problem.

For more information, see the provided video tutorials.

## Get the code

Clone this repository, or download an archive from the [Releases](https://github.com/mitchellmcm27/SkyStrat/releases). On the releases page, go to the latest version, click the arrow next to Assets, and click the source code in your preferred format to download. Then, unzip and use the instructions below to install.

## Arc GIS Plugin installation

Open the Catalog pane in ArcGIS Pro. Under "Folders", navigate to the **BuskMethod.pyt** file where you saved the download.
Right click, and choose *Add to Project*.

Open the Geoprocessing Pane, click Toolboxes between Favorites and Portal at the top of the pane, and the added toolboxes will appear under the Project heading. Expand the **Busk Method Toolbox** menu item, and click **Busk Downplunge View** to launch the tool.

To do three-point-problem calculations, a similar sequence can be followed with the **ThreePointProblemToolbox.pyt** file.

## QGIS Plugin installation

Installing a local plugin requires copying the files to the QGIS **python\plugins** directory.
The easist way to develop a plugin is to create a symlink, which allows changes to be automatically synced to the correct directory without manually copying files.

After copying or symlinking, restart QGIS, and then enable the "SkyStrat" plugin in this QGIS plugins manager.

![enable SkyStrat](enable_skystrat.png "Enable SkyStrat")

The plugins will be available in the **Processing Toolbox** under the *SkyStrat* menu item.

To uninstall, delete the `SkyStrat` symlink in the **QGIS3/profiles/default/python/plugins** folder.

### Windows symlink

cmd prompt (run as Admin):

```cmd
mklink /D "C:\Users\<username>\AppData\Roaming\QGIS\QGIS3\profiles\default\python\plugins\SkyStrat" "C:\Path\To\This\Repo\SkyStrat\qgis\SkyStrat"
```

Powershell (run as Admin):

```powershell
New-Item -ItemType SymbolicLink -Path "C:\Users\<username>\AppData\Roaming\QGIS\QGIS3\profiles\default\python\plugins\SkyStrat" -Target "C:\Path\To\This\Repo\SkyStrat\qgis\SkyStrat"
```

### Mac symlink (not yet tested)

```bash
ln -s /path/to/this/repo/SkyStrat/qgis/SkyStrat ~/Library/Application\ Support/QGIS/QGIS3/profiles/default/python/plugins/SkyStrat
```

### Linux symlink (not yet tested)

```bash
ln -s /path/to/this/repo/SkyStrat/qgis/SkyStrat ~/.local/share/QGIS/QGIS3/profiles/default/python/plugins/SkyStrat
```
