# Run FEBio simulations
Information on how to run the FEBio simulations.
For a given problem an FEBio definition file is used. The simulation can then be executed via:

Executable in `/home/mkoenig/git/buildroadrunner/FEBio/cbuild/bin/febio3`

```
export PATH=${PATH}:/home/mkoenig/git/porous_media/buildroadrunner/FEBio/cbuild/bin
```

```bash
febio3 "/home/mkoenig/git/porous_media/febio/lobule_BCflux.feb"
```

This creates a result file
```bash
/home/mkoenig/git/porous_media/febio/lobule_BCflux.xplt
```

To create the VTK files use FEBIO Studio >= 2.1.
https://github.com/febiosoftware/FEBioStudio/issues/48
https://febio.org/downloads/
https://febio.org/knowledgebase/getting-started/updating-febio-and-febio-studio/

To update to the latest development version, please follow these steps:

- Open FEBio Studio and go to the menu Tools\Options (on MacOS this is FEBio Studio\Preferences).
- In the dialog that appears, select the Auto-Update settings from the left-side panel.
- On the right-side panel, click the “Update to development version” and follow the instructions on screen.
