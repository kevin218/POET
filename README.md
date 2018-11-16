# The Photometry for Orbits, Eclipses, and Transists (POET)
POET pipeline for reduction of Spitzer data

Requirements

- Python >3.5

# Quick Start Guide:
1. git clone url
2. Place the following into your bash profile: 
    `export PYTHONPATH=$PYTHONPATH:/home/username/.../code/lib`
3. Copy .../code/run and its contents into your analysis folder and modify the POET control files (.pcf)
4. Read the documentation in .../code/doc
5. To run the first half of the POET pipeline, run the following sample code from within the run directory:
    * `python poet.py p1`
    * `python poet.py p2`
    * `python poet.py p3`
    * `python poet.py p4 fgc`
    * `python poet.py p5 fgc/ap2500715/`
6. To run the second half of the POET pipeline, copy `run_interactive.py` into .../fgc/ap2500715/{newFolder}
7. Create the 'ancil' directory in .../fgc/ap2500715/{newFolder} and copy `eg00-initvals.txt` and `eg00_params.py` into the directory
8. Open a Python3 session and follow the instructions inside `run_interactive.py`
