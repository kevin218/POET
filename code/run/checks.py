
#------------------------------------------------------------------
# Checkings

import sys, os
r = os.getcwd().split("/")
maindir = "/".join(r[:r.index("run")])
sys.path.append(maindir + '/lib/')
import manageevent  as me

# Load to see and check the results:
# change evtname to the name of your event.
event = me.loadevent('evtname_ini', load=['data','uncd','bdmskd'])
event = me.loadevent('evtname_bpm', load=['data','uncd','mask'])
event = me.loadevent('evtname_den', load=['data','uncd','mask'])
event = me.loadevent('evtname_ctr', load=['data','uncd','mask'])
event = me.loadevent('evtname_pht')

# Check visually the centering results:
import frameviewer as fv
event = me.loadevent('evtname_pht')
fv.frameviewer(event, zoom=True)   # zoom-in around the target.
fv.frameviewer(event, zoom=False)


# If by any chance you want to run the pipeline from an interactive session:
import poet_1event  as p1
import poet_2badpix as p2
import poet_3center as p3
import poet_4photom as p4
import poet_5checks as p5
import poet_denoise as pd

# Change evtname accordingly to your event name.
p1.Event(         'evtname.pcf')
p2.badpix(        'evtname_ini')
pd.run_denoising( 'evtname_bpm', 'denoise.pcf')
#p3.run_centering( 'evtname_bpm', 'center.pcf')
p3.run_centering( 'evtname_den', 'center.pcf')
p4.run_photometry('evtname_ctr', 'photom.pcf')
p5.checks(        'evtname_pht')


