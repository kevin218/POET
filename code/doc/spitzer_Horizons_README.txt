The following is an explaination of how js41_spitzer.vec was created.
You can follow the same procedure to create a similar file for any
other data set.  

In order to create js41_spitzer.vec or a similar file:
1. Go to the horizons webpage: http://ssd.jpl.nasa.gov/?horizons
2. You will see that the HORIZONS system can be accessed using any of the following methods: 
 telnet  
 email
 web-interface
Chose web-interface. 
3. You will see the default settings: 

Current Settings
Ephemeris Type [change] :  OBSERVER 
Target Body [change] :  Mars [499] 
Observer Location [change] :  Geocentric [500] 
Time Span [change] :  Start=2007-03-19, Stop=2007-04-18, Step=1 d 
Table Settings [change] :  defaults 
Display/Output [change] :  default (formatted HTML) 

(Note: Time span will change to reflect current date)

You want to change the settings so that they look like this (for js41 times are as stated, for other observations change the times accordining):


Current Settings
Ephemeris Type [change] :  VECTORS 
Target Body [change] :  Spitzer Space Telescope [SST SIRTF (Spacecraft)] [-79] 
Coordinate Origin [change] :  Solar System Barycenter (SSB) [500@0]
Time Span [change] :  Start=2006-12-28 00:00, Stop=2006-12-29 00:00, Step=1 h 
Table Settings [change] :  output units=KM-S; reference plane=FRAME; labels=YES
Display/Output [change] :  download/save (plain text file) 


Do this by: 
a) Changing the ephermeris type settings.  
Your chooses are:  
- Observer Table (Use this table type to generate a table of observer quantities (such as R.A./Dec.) for any object with respect to a geocentric or topocentric observer. )
 - Vector Table (Use this table type to generate a Cartesian state vector table of any object with respect to any major body. )
 - Orbital Elements (Use this table type to generate a table of osculating orbital elements for any object with respect to an appropriate major body.)
CHOOSE Vector Table
b) Changing the Taget Body settings.  
You will be given a search box, type in "Spitzer" and search either all bodies or major bodies, it should come up with "Spitzer Space Telescope [SST SIRTF (Spacecraft)] [-79]" 
c)Changing the Coordinate Origin.
in the 'Lookup Named Location' search box, type @0 
This should give you "Solar System Barycenter (SSB) [500@0]"
d) change the Time Span Settings
JS41 was taken on 2006-12-28 centered around 15:00 UT, so I chose to get horizons data for the entire day (you need atleast a 3-4 hours on either side) in 1 hour intervals.
e) Changing the Table Settings.
-change the output units to "km & km/s"
-change the reference plane to "Earth mean equator and equinox of referece epoch
-find labels : (box)  -- enable labeling of each vector component, and check the box
-find object page : (box)  -- include object information/data page on output, the box should be checked by default , but if not checked, check it.        
f) Display/Outpit settings.
Either change to download/save (in which case it will prompt you name the file when you go to "Generate Ephemeris") 
Or leave as is and copy and paste the information from the browser. 

When all the settings have been promperly changed Press the "Generate Ephemeris" button.  


*** Update 6/3/07:
Formerly we used the setting
Coordinate Origin [change] :  Sun (body center) [500@10]


If you have any additional questions contact Karen Horning (karen.horning@gmail.com)
