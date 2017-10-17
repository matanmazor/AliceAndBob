#experiment 
#last updated on May 8th, 2017

from psychopy import core, visual, monitors, event
from psychopy.visual import ShapeStim # for drawing the arrows
from lib.preRNG import preRNG # pre-RNG of protocol folder
from random import shuffle # set block order
import os, sys  # handy system and path functions

def closeAll(win, log, premature = 0):
    win.close()
    log.write("%-18s%f\n" % ('end', start.getTime()))
    log.close()
    if premature:
        print os.path.join(logDir, 'run'+str(runNum)+'terminatedByUser')
        os.rename(os.path.join(logDir, 'run'+str(runNum)),
                  os.path.join(logDir, 'run'+str(runNum)+'terminatedByUser'))
    core.quit()
    sys.exit(0)  

# Timing-
initialBlank = 10 #seconds
blockDuration = 8 #seconds
restDuration = 10 #seconds
cycleDuration = blockDuration + restDuration #seconds

# Set order of events
blocks = ["narf"] #a place keeper
preRNG("path to protocol folder")
for runIdx in [1,2,3,4]:        
    blockOrder = [0]*10+[1]*10
    shuffle(blockOrder)
    blocks.append(blockOrder)

# Prompt user for run number
runNum = input("please enter the run number: ")

_thisDir = os.getcwd()
logDir = os.path.join(_thisDir,"log")
if not os.path.isdir(logDir):
    os.mkdir(logDir)
log = open(os.path.join(logDir, 'run'+str(runNum)),'w')
        
# Get monitor parameters
mon = monitors.Monitor('testMonitor')
screenSize = mon.getSizePix()

# Create the main windows
win = visual.Window(monitor='testMonitor', size=screenSize, units='cm',
                    color='black', fullscr = True)
win.mouseVisible = False

# Create shapes
rightArrow = [(-.2,.05),(-.2,-.05),(0,-.05),(0,-.1),(.2,0),(0,.1),(0,.05)]
rightArrow = ShapeStim(win, vertices=rightArrow, fillColor='white', size=4)

leftArrow = [(.2,.05),(.2,-.05),(0,-.05),(0,-.1),(-.2,0),(0,.1),(0,.05)]
leftArrow = ShapeStim(win, vertices=leftArrow, fillColor='white', size=4)

fixation = [(0, -0.5), (0, 0.5), (0,0), (-0.5,0), (0.5, 0)]
fixation = ShapeStim(win,vertices = fixation, lineWidth=5, closeShape = False,
                     lineColor="white")

arrows = (leftArrow, rightArrow)
arrowNames = ('left','right')

# Wait for t
template = "The %s scan is just about to start"
serialNumbers = ("Narf","first","second","third","fourth")
message = visual.TextStim(win, text=template % (serialNumbers[runNum]))
message.draw()
win.flip()
while not event.getKeys(['t','T']):
    pass

# EXPERIMENT STARTS NOW
start = core.Clock()
fixation.draw()
win.flip()
while start.getTime()<initialBlank:
    if event.getKeys(['escape']): 
        closeAll(win, log, 1)
for eventIdx, eventType in enumerate(blocks[runNum]):
    arrows[eventType].draw()
    win.flip()
    log.write("%-18s%f\n" % (arrowNames[eventType], start.getTime()))
    while start.getTime() < initialBlank + eventIdx*cycleDuration + blockDuration:
        if event.getKeys(['escape']): 
            closeAll(win, log, 1)
    fixation.draw()
    win.flip()
    log.write("%-18s%f\n" % ('rest', start.getTime()))
    while start.getTime() < initialBlank + (eventIdx+1)*cycleDuration:
        if event.getKeys(['escape']): 
            closeAll(win, log, 1)
closeAll(win, log)

