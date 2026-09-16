# README #
This Simple.h, Simple.cpp are LAMMPS scripts for the simplified lubrication force. The theory (or rather the algebra) behind the simplification is given in the pdf file.

### What is this repository for? ###

* Calculating lubrication force between two particles using [LAMMPS](https://github.com/lammps).
* Calculating Stokes drag for every particle using LAMMPS.

### How do I get set up? ###

* Copy Simple.h and Simple.cpp in LAMMPS installation's src folder and recreate the LAMMPS executable.
* It depends on pair_lubricate.h that is included with LAMMPS.
* All the flags are the same as [lubricate/poly](http://lammps.sandia.gov/doc/pair_lubricate.html). flagfld now stands for Stokes drag. flagVF is left for legacy reasons and will not affect the results when set to 0 or 1.
* Ensure that the inner and outer cutoffs are greater than sphere diameter, and use pair_modify mix arithmetic for spheres of different sizes. See collision.lmp example.
* Modify the two particle collison code: collision.lmp, and compare it using the matlab scripts provided.

### Contribution guidelines ###

* Please email Ranga for updates or features needed. You are welcome to modify and use this code to your specs.

### Who do I talk to? ###

* Ranga (r.radhakrishnan@ed.ac.uk)


### List of things to work on for future updates ###

1. Make the class to be independent of pair_lubricate files.
2. Make the cutoffs with arithmetic mean by default.
3. When Stokes drag is included using flagfld=1, change v_stream velocities based on the bottom edge of the box. Currently v_stream assumes the box bottom edges to be at 0, it needs to depend on the actual box dimensions. Note that the v_stream in both the cases will only differ by a constant, and will not have much impact on the steady state results.

### Acknowledgements ###
Tim Najuch (University of Edinburgh) for leading the discussion on this topic and helping me develop this code.
Dr. Chris Ness (Cambridge University)for verifying the theory behind the forces.
Prof. Jin Sun (University of Edinburgh) for his kind support and encouragement.
EPSRC for funding under Future Formulation of Complex Products.