  I tried to use potentials in SMASH for Snn  energies 7.7 GeV and higher.
  But the same error occured during event generation.
  With an energy of 5 GeV there was no error.
  In this folder there are files with logs, run file and config file.

   6 log files are made with different potentials and energies 7.7 and 11.5 GeV Au+Au
   Potential parameters are taken from the article: 
http://inspirehep.net/record/1664230/files/PoS(EPS-HEP2017)178.pdf
   
   With fermi motion "off" or "frothen" all runs well (potentials off). When I "on" fermi motion (at User Guide it
 is recommended to also activate potentials) an error occurs. The error always is the same: 



FATAL         Fpe        : Floating point trap was raised: Invalid (domain error occurred)
stack trace:
  /lib/x86_64-linux-gnu/libpthread.so.0 : ()+0xf890
  /lib/x86_64-linux-gnu/libm.so.6 : ()+0x15543
  /lib/x86_64-linux-gnu/libm.so.6 : pow()+0x1c
  /lustre/nyx/hades/user/parfenov/Soft/Models/SMASH/smash/build/smash : smash::Potentials::skyrme_force(double, smash::ThreeVector, smash::ThreeVector, smash::ThreeVector) const+0x36
 ....
At log files full list.
