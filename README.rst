Installation
------------

`Stack`_ is the easiest way to build and run the simulation.
First you need to install `stack`_. Installation depends on the platform.
Most linux distributions should provide an easy way to get `stack`_ from
repositories. Example installation on Debian 8 x64::

    sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 575159689BEFB442
    echo 'deb http://download.fpcomplete.com/debian jessie main'|sudo tee /etc/apt/sources.list.d/fpco.list
    sudo apt-get update && sudo apt-get install stack -y

Detailed instructions can be found in `stack's documentation`_.

When you have stack you can proceed with building the simulation.
The following commands do everything for you::

    git clone https://github.com/Motions/prototype.git
    cd prototype
    stack setup
    stack build

Some dependencies require a fair amount of RAM when being built. Our experiences
tell that at least 3 GB is necessary (use swap if needed).

.. _stack: http://docs.haskellstack.org/en/stable/README.html
.. _stack's documentation: http://docs.haskellstack.org/en/stable/README.html#how-to-install

Running the simulation
----------------------

You can either use stack (if you used stack build)::

    stack exec -- prototype arguments

Or simply execute it from the directory containing the executable (the path may vary a little)::

    cd .stack-work/dist/x86_64-linux/Cabal-1.22.4.0/build/prototype/
    ./prototype arguments
    
To run the simulation we need an initial state. To initialize the simulation run it with the "init" command::

    ./prototype init 

The detailed usage will be shown. Example initialization could be::

    ./prototype init lbsites.txt rbsites.txt -r 20 -b 100 -l 100 -o init-state
    
Then you can use the output file as an input file to the "run" command. To show its usage::

    ./prototype run
    
An example run::

    ./prototype run -p results.pdb -i init-state -o end-state -s 10000
    
The "-p" argument file (in this example results.pdb) will contain the final state of the simulation in the PDB format. If you want it to contain all the intermediate states, use the "-m" flag. The "-o" argument file can be used as an input file to another run if you want to continue the simulation later.
