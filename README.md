# Free energy of hairpin formation

## Introduction

### Force Field

This tutorial is designed to introduce free energy calculation using the oxDNA2 force field. The oxDNA force field consists of DNA nucleotides, each represented by an anisotropic bead. Each nucleotide interacts with others through:

1. excluded volume
2. backbone connection 
3. phosphate repulsion
4. stacking
5. cross stacking
6. hydrogen bonding 

The values of the interaction parameters were taken from a variety of experiments, including: structural information (i.e. B form DNA crystals), thermodynamic information (agreement with the Santa-Lucia 1998 parameters for nearest neighbours DNA thermodynamics), and some information relating to the angle of DNA origami.  

The dynamics of the model are *confusing*, but since we will be running only Monte Carlo simulations, we don't have to worry about them at all.

NOTE: I'm using the oxDNA force field here. But there's actually an oxDNA2 force field that's a bit more realistic. To run using the oxDNA2 force field, use `interaction_type=DNA2`, and `salt_concentration=0.5` in your input file.

### Implementation

oxDNA is implemented in custom C++ with CUDA support for GPUs. Installation of the CUDA support is a tiny bit tricky, but we'll (probably) be doing just Monte Carlo simulations, so we won't need GPU support. 

oxDNA is also available as an extension of LAMMPS, which supports MPI. You could also implement it in HOOMD if you wanted (for example) multi-GPU acceleration.

### Example problem

A very simple but relevant example problem is to evaluate the free energy of a DNA hairpin loop closing as the temperature varies. 

Imagine you had a sequence consisting of 

`5'-ATGCTTTTTTTTGCAT-3'`

This could be either in the 'open' form above, or could self-hybridize, forming loop (i.e. the 5'-ATGC-3' region could hybridize to the 5'-CGAT-3' region). We might (for example) be interested in the free energy difference between those two states as a function of temperature.

This will introduce the following ideas:

1. Initialization
2. Trajectory inspection
3. Reaction coordinates
4. Biased Monte-Carlo sampling

We can also compare the oxDNA predictions to those of a simple model: beads on a Gaussian spring. 

I'm estimating this will take about 90 minutes.

## Instructions

Follow the oxDNA instructions to install oxDNA.

Also you'll be more efficient when typing if you make sure you can call oxDNA from the terminal. So either put an alias in your `~/.bashrc`, or append the oxDNA installation directory to your PATH variable.

You should also make sure the `UTILS` folder of oxDNA is on your path, and that you have python2.7 installed.

### Initialisation

Create a file called seq, and add the line:

`TAGCTTTTTTTTGCTA`

Then run `generate-sa.py 50 seq`.

This means generate a starting configuration in a box side of 50 oxDNA units (each about 0.85 nm), using the sequence in seq. Due to some historical problems, oxDNA likes sequence files written 3' to 5'.

This will generate a file called generated.dat and generated.top.

Use oxDNA viewer to have a look at it (just drag into the oxDNA viewer window).

Then have a look at the .dat file.

It has a line for time step, box size, energy (the triple of internal energy, kinetic energy, and total).
Then it has a load of lines of 15 space separated floats.
These correspond to:

1. x,y,z of center of mass
2. phosphate-to-base vector (they call this a versor)
3. vector which points along the spine of the DNA helix (i.e. normal to the base)
4. vx, vy, vz center of mass velocity
5. angular momentum as a three vector

### Running a basic simulation

Move to the `2.raw-monte-carlo directory`. Copy the initially generated files from earlier (generated.dat and generated.top), and move them to this directory.

Have a look at the input file. I've annotated this one so you can see what the features mean.

This input file is special because it prints out the configuration of hydrogen bonds in the system at any point in time.

When we run the simulation we get a lot of things printed to the terminal. These are:

Time step; internal energy; translational acceptance prob. ; rotational acceptance prob; don't know

### Viewing oxDNA trajectories

To view the trajectory the best viewer is the Javascript one developed at Petr Sulc's lab.

`https://sulcgroup.github.io/oxdna-viewer/`

### Biased Monte Carlo

We want to find out the free energy difference between the fully bonded and fully unbonded states. Ideally we would just run the simulation for ages, and then look at the occupancy of these two states: maybe it is unbound 90% of the time, or maybe it is unbound 20% of the time.

However, the "wait and see" approach (i.e. unbiased Monte-Carlo) ends up being very inefficient. It's inefficient because the transition state between the unbound and bound states is very unfavorable in terms of free energy. The transition state of only one or two bonds being formed pays the entropic price of constraining the loop, but only gets limited internal energy benefit. So transitions between the two free energy minima are rare. So we don't see many transitions and therefore we collect poor statistics about relative free energy.

At the extreme level, the transition state is so unfavorable that we may never see transitions! So the bonded form would live forever being bonded, and the unbonded form would never bind. So obviously it would be impossible to calculate the free energy differences.

To solve this, we use biased Monte Carlo. We add an additional bias, and correct for it later. Here we want to bias the existence of the transition state. When the transition state is favorably biased, we will see lots and lots of transitions between the two free energy minima.

The reaction coordinate I have chosen to bias along is:

when there are no hydrogen bonds:

(1) binned minimum distance between all complementary nucleotides

when there is at least one hydrogen bond:

(2) the number of hydrogen bonds in the system

This works for the same reason you can work out the free energy as a function of length of a *macroscopic spring*, by just applying a force.

### Biased Monte Carlo in oxDNA

oxDNA supports biased Monte-Carlo using square well potentials. To do biased Monte Carlo, move into directory `2.2-biased-monte-carlo`. 

Have a look at the input file. You can see that there are some additional input commands compared to the normal Monte Carlo simulation. 

This commands are:

1. `umbrella_sampling = 1` turns umbrella sampling on
2. `op_file = op.txt`  #the file that defines the order parameter
3. `weights_file = wfile.txt` #the file that defines the artificial biases.
4. `safe_weights = 0` # I will talk about this later.
5. `default_weight = 0` # same

In the directory is a file called `op.txt`, which contains the order parmeter information. For the bonding, I have added all the different complementary bonds in the system I care about. In the mindistance, I have added the same, but also added some interfaces. Since the oxDNA biasing support just uses square wells, we need to fine some way of discretizing the space. Note that mindistance means the minimum distance over all pairs mentioned here. Similarly the number of bonds means the total number -- obviously there are four different ways to get three bonds (so you may think about what that means  for the configurational entropy as a function of number of bonds...)

In the file `wfile.txt`, each line consists of two space separated ints, and one float. 

The line:

`0 3 70`

means: when the number of bonds is 0, and the minimum distance is in bin 3 (i.e. greater than 4), bias with a Bolzmann factor of 70.

The line:

`4 0 57`

means: when the number of bonds is 4, and the mininum distance is in bin 0 (i.e. less than 1.5), bias with a Bolzmann factor of 57.

(notion of true reaction coordinate)

When we run a biased Monte Carlo sampling run in oxDNA, we also print additional columns correponding to the order parameter(s), and the bias.

### Running biased Monte Carlo

In `2.2-biased-monte-carlo/ENERGY0`, I have set up the biased Monte Carlo simulation. You can run it from there. Here I just chose some approximate guesses about what the bias should be to encourage transitions. 

When practically doing biased Monte-Carlo, you obviously need to start *without* a knowledge of a decent bias. Instead you guess, run the simulation, then use this inefficiently biased simulation to update the biases, and then iterate on this repetedly. In `ENERGY1` through `ENERGY5`, I have run (fairly brief) simulations in parallel (i.e. trivially in parallel -- just rerun with different seeds). Then in the analysis folder I have a Jupyter notebook where I have evaulated the predictions of free energy, as a first pass.

### Improving the bias

The first bias was a bit rough, and led to probably inefficient sampling. With the updated free energies, I can bias the next set of simulations much more efficiently! I have done that in `2.3-biased-monte-carlo-optimal`. In the analysis folder there I also have a Jupyter notebook where I calculate a much better free energy estimate. That is nominally the end of the tutorial.

## Automated approaches

It seems like the job of choosing the correct bias could be automated. There are various techniques for doing this -- the grown up way to do this is either Wang-Landau sampling, or Metadynamics. If you have a very big and complicated transition path, you could also do umbrella sampling, where you look at free energy differences between tiny bits of the reaction coordinate with separate simulations, then stich them all together. There's a book by Dan Frenkel on my desk you're welcome to borrow. This goes over everything.

## Extension

There are lots of extra question you could now ask:

1. Can I extract the entropy and internal energy profiles along the reaction coordinate?
2. How does this distribution change with temperature?
3. How does this distribution change with length of the single stranded polyT region?

LDM told me that you were going to look at vesicle fusion mediated by DNA tethers (like SNARE protein mimics). But he didn't tell me what he wanted your simulation project to be on. I hope this gives you a decent base to work off. I'm always available if you need a hand.

Cheers,

Will








