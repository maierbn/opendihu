# uni-stuttgart-beamer-template

This is an unofficial L<sup>A</sup>T<sub>E</sub>X template for Beamer
presentations by and for members of the University of Stuttgart,
Germany, following the new corporate design the university has given
itself in 2016.

The template is private work by members of the University of Stuttgart.
It is not endorsed by the University of Stuttgart or any of its
institutions, institutes, or departments.
The template is made available to the public as we think it might
be useful for other people.
Neither the authors nor the University of Stuttgart provide support
for the template going beyond this `README`.

The license for the template is located in `LICENSE`.

If you happen to improve the template, it would be nice if you merged
back your work using pull requests.

## Usage

The directory `tex` contains L<sup>A</sup>T<sub>E</sub>X source code
and resulting PDFs, while `gfx` contains graphics.
The main component of the template is the new Beamer theme "Stuttgart",
located in `tex/beamerthemeStuttgart.sty`.
Also included is an example presentation in `tex/talk.pdf`
(L<sup>A</sup>T<sub>E</sub>X source in
`tex/talk.tex` and `tex/slides.tex`).

### SCons

The preferred way to compile the example presentation is using SCons.
SCons is a build system tool (equivalent to CMake, Autotools, etc.),
which can also be used to build L<sup>A</sup>T<sub>E</sub>X
documents.

On Ubuntu, install SCons via
```
sudo apt install scons
```
Then, download and extract this repository to a folder.
Change to this folder and execute
```
scons
```
and you're done.
This not only compiles the presentation using the right amount of
L<sup>A</sup>T<sub>E</sub>X calls (very much similar to `latexmk`),
but generates also a separate handout in `tex/handout.pdf`,
which contains the slides in Beamer's `handout` mode
(e.g., ignoring `pause` commands) compressed with 4 pages per sheet.
You can clean up all generated files using
```
scons -c
```

### PDFL<sup>A</sup>T<sub>E</sub>X

Of course, you can also compile the example presentation in the
traditional way using
```
cd tex
pdflatex talk.tex
```
