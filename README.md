# The Photometry for Orbits, Eclipses, and Transists (POET)
POET pipeline for reduction of Spitzer data

Requirements

- Python >3.5

# Quick Start Guide:
1. git clone url
2. Change project_name to the desired name of your project.
3. In the project_name folder, create two directories named analyses and data. 
The data directory is where the downloaded exoplanet data will exist. ex. Wa077bo21/S19.2.0/files
The analyses directory is where commands will be executed. ex. Wa077bo21/run/poet.py
  
4. Place the following into your bash profile:
  export PYTHONPATH=$PYTHONPATH:/home/username/....../projectname/code/lib/python
  export PYTHONPATH=$PYTHONPATH:/home/username/...../projectname/code/python/pipeline/trunk
  
5. To run the first half of the POET pipeline, run the following code from withing the run directory. 
6.  python poet.py p1
7.  python poet.py p2
8.  python poet.py p3
9.  python poet.py p4 fgc
10. python poet.py p5 fgc/ap2500715/
