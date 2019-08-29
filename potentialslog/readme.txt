  I tried to use potentials in SMASH for Snn  energies 7.7 GeV and higher.
  But the same error occured during event generation.
  With an energy of 5 GeV there was no error.
  In this folder there are files with logs, run file and config file.

   6 log files are made with different potentials and energies 7.7 and 11.5 GeV Au+Au:
#           a [MeV]   b [MeV]     t    S_pot [MeV]
#
#default    -209.2    156.4     1.35    18
#
#soft       -356.0    303.0     1.17    18
#
#hard       -124.0    71.0      2.0     18
#Energy Snn 7.7 or 11.5 
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
  /lustre/nyx/hades/user/parfenov/Soft/Models/SMASH/smash/build/smash : smash::Potentials::all_forces(smash::ThreeVector const&, std::vector<smash::ParticleData, std::allocator<smash::ParticleData> > const&) const+0x244
  /lustre/nyx/hades/user/parfenov/Soft/Models/SMASH/smash/build/smash : smash::update_momenta(smash::Particles*, double, smash::Potentials const&, smash::RectangularLattice<std::pair<smash::ThreeVector, smash::ThreeVector> >*, smash::RectangularLattice<std::pair<smash::ThreeVector, smash::ThreeVector> >*)+0x5fb
  /lustre/nyx/hades/user/parfenov/Soft/Models/SMASH/smash/build/smash : smash::Experiment<smash::ColliderModus>::run_time_evolution()+0x396
  /lustre/nyx/hades/user/parfenov/Soft/Models/SMASH/smash/build/smash : smash::Experiment<smash::ColliderModus>::run()+0x191
  /lustre/nyx/hades/user/parfenov/Soft/Models/SMASH/smash/build/smash : main()+0x1952
  /lib/x86_64-linux-gnu/libc.so.6 : __libc_start_main()+0xf5
  /lustre/nyx/hades/user/parfenov/Soft/Models/SMASH/smash/build/smash() [0x44a997]

