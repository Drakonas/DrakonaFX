# DrakonaFX

DrakonaFX is a set of customizable Shader presets that can be used with [Dolphin Ishiiruka](https://forums.dolphin-emu.org/Thread-unofficial-ishiiruka-dolphin-custom-version) for the purpose of extended effects originating from ReShade and DolphinFX.

## History

These are the custom shader presets originally created by the visual effects team of Hamasaki Productions, originally designed for use with our Twilight Princess Gamecube Mod project.

While these shaders were originally designed for our project, we soon realized they worked quite well on many other games depending on what the user liked. So we decided to release them as a pack separately.

This pack also includes our later released Modern shader presets, which are very versatile.

Thanks go to Tino for Ishiiruka and creating a custom GLSL that allows more capability with the shaders.

Credit for piecing together the code goes to Drakonas, who packaged the GLSL's for release, and brought many effects from ReShade for use with Ishiiruka with help from Tino.

I should note that the Gaussian Anamorphic LensFlare now found in the official builds of Ishiiruka were reworked by Tino from my original code 2 years back.

We hope you like these shaders, and now on with the installation!

## Requirements
- [Ishiiruka](https://forums.dolphin-emu.org/Thread-unofficial-ishiiruka-dolphin-custom-version) - The latest was recommended previously, but now stable build 758 is recommended until I update the code. Ishiiruka is required to use DrakonaFX. These shaders won't work with normal Dolphin.
- A gaming laptop/PC within the last 8 years - Depending on how decent your rig is, we have made multiple versions of each preset to allow for high and low performing machines alike.

As a basic guide on what hardware is required I can say from our tests that the minimum needed for running the
Ultra preset is dependant solely on the graphics card. From what we have seen, an AMD R9 280X/290 or Nvidia GTX 780 Ti/970
or higher should be enough to guarantee 30 FPS solid throughout the game with no visual defects or framerate dips below 30.

As far as the lower presets go, performance increases dramatically. In retrospect, a 970 with proper graphics settings runs Ultra
at 1080p at around 30-35 FPS with an unlocked framerate, while High runs at about 45-60 FPS, and Medium runs past 70-90 FPS, and Low past 90-100 FPS.

We've seen that even a GTX 640 1GB can run High at 720p, and a GTX 670 or AMD Radeon 7950 can run High at 1080p consistently, and beyond that the
presets run on most gaming computers with a dedicated graphics card that isn't more than 8 years old.

## Installation

Installation instructions differ depending on if there is a portable.txt in the folder where you extracted Ishiiruka.

1. Portable.txt exists: Place GLSL files in User/Shaders/Postprocessing, residing in the same folder as Dolphin.exe
2. Portable.txt missing: Place GLSL files in Documents/Dolphin Emulator/Shaders/Postprocessing.

FYI: When updating Dolphin/Ishiiruka, it is highly recommended to start with a fresh Sys directory or new folder altogether.
When updating Dolphin or Ishiiruka, if you do not erase the Sys directory, there could be files left over that are not intended for the new build of the emulator.

## Usage

The shaders will now show up in Options -> Graphics Settings -> Post-Processing. You can add the shader you want there. The options, if they have never been changed within your installation before, will be set to the defaults as we designed them.

Some games require different Post-Processing Triggers to apply the shaders properly. The most common one that works is On EFB Copy, but for Twilight Princess On Swap or On Project is recommended, for example. Also, sometimes it depends on the scene being shown. On Projection is the only reported trigger to apply the shaders to the particle effects in the Goron Mines for example.

You can change the settings for the preset if you wish to by selecting the highlighting the shader that has already been added and clicking "Options..". If you wish to revert to the defaults we designed for the preset, you can click "Restore Defaults" in the "Options.." window.

Warning: Opening Graphics options while the game is running can produce unwanted results. We recommend changing other settings while the game is not running, but if you wish to edit shaders while the game is running, make sure to have Options -> Graphics Settings -> Advanced -> Compile Shaders on Startup turned off, or else your changes may not propogate. But make sure to turn the option back on when you are done. The emulator runs slower with the shaders not completely compiled on startup.

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
[GNU GPL v2](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)