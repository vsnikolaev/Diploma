# Prerequisites

SMASH (Simulating Many Accelerated Strongly-interacting Hadrons) is a
relativistic hadronic transport approach for the dynamical description of heavy
ion reactions.
Please see [Phys. Rev. C 94, 054905
(2016)](https://arxiv.org/abs/1606.06642) for details and cite this reference,
if you are using SMASH.

To get started, install the application from the site [SMASH](https://github.com/smash-transport/smash). 
Users guide [here] (http://theory.gsi.de/~smash/userguide/1.6/index.html)

# Running SMASH

To run SMASH iprut `./smash` at ...[Path]/SMASH/build
The default files are called: configuration file `config.yaml`, particles file `particles.txt` and decaymodes file `decaymodes.txt`.
To rut with files from other locations input:

  ./smash -i [Path]/config.yaml
  ./smash -p [Path]/particles.txt
  ./smash -d [Path]/decaymode.txt

Note, that files can be named as you like.
To change the output directory use option `-o`
 
  ./smash -i [Path]/config.yaml -o [Other_Path]/

All line options can be viewed with `-h`
To rut SMASH silently redirect stdout to /dev/null. Warning and error will be displayed

./smash > /dev/null
