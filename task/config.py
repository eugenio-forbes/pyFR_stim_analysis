# pyFR Configuration

# for the stimulation session
isPrepSession = False
if isPrepSession:
    # NO STIM
    #numTrials   = 25
    numTrials = 20
    stimType = [0]
    #numSessions = 4
    numSessions = 5
#listLen     = [2,3,4,5,6]
    listLen     = [9,10,11,12,13]
else:
    # STIM
    #numTrials=50
    numTrials   = 20
    numSessions = 5
    stimType = [0]
    #listLen     = [3]
    listLen     = [10]

# Pause+Jitter after orienting stim before first word
orient_Duration                = 1000
PauseBeforeWords               = 1000
JitterBeforeWords              = 1000
PauseBeforeRecall              = 500
JitterBeforeRecall             = 200
PauseBtwnStimAndFirstWordStart = 500

# Word Font size (percentage of vertical screen)
wordHeight = .1

# Duration word is on the screen
#word_Duration = 750
word_Duration = 1500
#word_ISI      = 200
word_ISI      = 800
word_Jitter   = 50

# Duration of recall in ms
recallDuration     = 20000
terminateRecallKey = 'SPACE'

# Beep at start and end of recording (freq,dur,rise/fall)
startBeepFreq     = 800
startBeepDur      = 500
startBeepRiseFall = 100
stopBeepFreq      = 400
stopBeepDur       = 500
stopBeepRiseFall  = 100

# Orienting Stimulus text
orientText      = '+'
recallStartText = '*******'

# Math distractor options
doMathDistract = True
MATH_numVars = 3
MATH_maxNum  = 9
MATH_minNum  = 1
MATH_maxProbs = 20
MATH_plusAndMinus = False
MATH_minDuration_Practice = 30000
MATH_minDuration = 10000
MATH_textSize    = .1
MATH_correctBeepDur = 500
MATH_correctBeepFreq = 400
MATH_correctBeepRF = 50
MATH_correctSndFile = None
MATH_incorrectBeepDur = 500
MATH_incorrectBeepFreq = 200
MATH_incorrectBeepRF = 50
MATH_incorrectSndFile = None

# Word pool to use
wp                    = 'pools/word-pool.txt'
wp_imagery            = 'pools/tor_imag_only.txt'
wp_medianimagVal      = 5.4
presentationType      = 'text'  # image, sound, text
presentationAttribute = 'name'


# Instructions text file
practice_instruct = 'text/fr_practice_instruct.txt'
instruct          = 'text/fr_instruct.txt'

# Default font
defaultFont = 'fonts/Verdana.ttf'

# Realtime configuration
# ONLY MODIFY IF YOU KNOW WHAT YOU ARE DOING!
# HOWEVER, IT SHOULD BE TWEAKED FOR EACH MACHINE
doRealtime = False
rtPeriod   = 120
rtComputation = 9600
rtConstraint  = 1200

# stim parameters 
pulseSignal = "1111111111111111"

pulsestim_length = PauseBtwnStimAndFirstWordStart
pulsestim_label  = 'ecog_stim'
pulsestim_signal = "1111111111111111"

pulseWidth = 10 # ms, probably shouldn't be changed
pulseFrequency = 100 # Hz
#pulseDuration = 3000 # ms
pulseDuration = 23000







