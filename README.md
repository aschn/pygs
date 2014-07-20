pygs
====

A visualizer and related tools for Anna's thylakoid simulation data.


Install
-----

1. Clone this repo and navigate to its top-level directory.
2. Install wx for your platform: http://www.wxpython.org/download.php
3. Install other dependencies: `pip install -r requirements.txt`


Quickstart
----

On the command line, run `./pygs/visualizer.py`

A GUI should appear.

Click "open files" and select particle coordinate files.
A sample configuration file is at `tests/data/particle_config_L400.dat`.

Clicking on a particle will toggle its tagged state, and write the particle identity to the terminal.

You can quit by clicking the close button, or by control-C on the command line.

Useful controls
------

* play: Play a movie starting from the current timepoint.
* draw PS frame: Draw the current screen as a postscript file.
* layers: Toggle the visibility of the bottom, top, and overlay views.
* movie controls > zoom: Set the zoom level of the view.
* initialize data > open files: Select particle coordinate files to load.
* initialize data > width / height: Set the width and height of the box, in nanometers. (This may be automatically set from the configuration file name.)
* colors: Set the base colors of tagged and untagged LHCII and PSII particles.
* color schemes: Toggle other color schemes. Play around with this one!
* tag within a rectangle: Tag all particles within a rectangle enclosed by xmin, xmax, ymin, and ymax, in nanometers.
* tag within a circle: Tag all particles within a circle centered on x, y with radius radius, in nanometers.
