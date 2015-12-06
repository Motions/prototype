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
