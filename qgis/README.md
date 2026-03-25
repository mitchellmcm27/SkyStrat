# QGIS Plugin installation

Installing a local plugin requires copying the files to the QGIS `python\plugins` directory.
The easist way to develop a plugin is to create a symlink, which allows changes to be automatically synced to the correct directory without manually copying files.

## Windows

cmd prompt (run as Admin):

```cmd
mklink /D "C:\Users\<username>\AppData\Roaming\QGIS\QGIS3\profiles\default\python\plugins\sky_strat" "C:\Path\To\This\Repo\SkyStrat\qgis\sky_strat"
```

Powershell (run as Admin):

```powershell
New-Item -ItemType SymbolicLink -Path "C:\Users\<username>\AppData\Roaming\QGIS\QGIS3\profiles\default\python\plugins\sky_strat" -Target "C:\Path\To\This\Repo\SkyStrat\qgis\sky_strat"
```

## Mac (TODO)

## Linux (TODO)