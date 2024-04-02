#!/usr/bin/python

from pyepl.locals import *
from pyepl.hardware import addPollCallback

# other modules
import random
import time
import os
import sys
import shutil


# Set the current version
VERSION = '0.1.1'
MIN_PYEPL_VERSION = '1.0.0'

class TupleText(Text):
    """
    Subclass Text so we can get offset time from present().
    """
    def __init__(self, text, font = None, size = None, color = None):
	Text.__init__(self, text, font, size, color)
	
    def present(self, clk = None, duration = None, jitter = None, bc = None, minDuration = None):
	v = VideoTrack.lastInstance()
	
	# get the clock if needed
	if clk is None:
	    clk = exputils.PresentationClock()
	    
	# show the image
	t = v.showCentered(self)
	onset = v.updateScreen(clk)
	
	if bc:
	    # wait for button press
	    button,bc_time = bc.waitWithTime(minDuration,duration,clk)
	else:
	    clk.delay(duration,jitter)

	v.unshow(t)
	offset = v.updateScreen(clk)
	
	if bc:
	    return onset,button,bc_time
	else:
	    return onset, offset

def textInput(screenText,video,keyboard,clock):

    # set up keyboard entry
    ans_but = keyboard.keyChooser()

    done = False
    while not done:

        rstr = ''
        field = video.showProportional(Text(screenText),.4,.5)
        input = video.showRelative(Text(rstr),RIGHT,field,20)
        video.updateScreen(clock)

        kret,timestamp = ans_but.waitWithTime(maxDuration = None, clock=clock)
        while kret:
            
            # process the response
            if kret.name == 'BACKSPACE':
                
                # remove last char
                if len(rstr) > 0:
                    rstr = rstr[:-1]

                # update text
                input = video.replace(input,Text(rstr))
                video.updateScreen(clock)
                kret,timestamp = ans_but.waitWithTime(maxDuration = None,clock=clock)
                
            elif kret.name not in ['BACKSPACE','RETURN','ENTER']:
                newstr = kret.name.strip('[]')
                rstr = rstr + newstr

                # update the text
                input = video.replace(input,Text(rstr))
                video.updateScreen(clock)
                kret,timestamp = ans_but.waitWithTime(maxDuration = None,clock=clock)
                
            elif kret.name in ['RETURN','ENTER']:
                video.clear("black")
                return rstr

def listWAScheck(wordInds, WASthresh):
    """
    Check if similarity between any two words in a list exceeds some
    threshold.
    """
    # check to make sure no two words in the list are too similar
    listGood = True
    for word1 in wordInds:
        for word2 in wordInds:
           val =  semMat[word1][word2]
           if val >= WASthresh and val < 1:
               listGood = False
               return listGood
    return listGood

def verifyFiles(config):
    """
    Verify that all the files specified in the config are there so
    that there is no random failure in the middle of the experiment.

    This will call sys.exit(1) if any of the files are missing.
    """
    
    # make the list of files from the config vars
    files = (config.wp,             
             config.instruct,
             config.defaultFont)
    
    for f in files:
        if not os.path.exists(f):
            print "\nERROR:\nPath/File does not exist: %s\n\nPlease verify the config.\n" % f
            sys.exit(1)

def prepare(exp, config):
    """
    Prepare the trials...
    """
    
    # verify that we have all the files
    verifyFiles(config)
    
    # get the state
    state = exp.restoreState()
    
    # load the word pool
    wp      = Pool(config.wp)
    wp_imag = Pool(config.wp_imagery)
    wp_hiImag = []
    wp_loImag = []
    count=0
    for allW in wp:
        if float(wp_imag[count].name)<config.wp_medianimagVal:
            wp_loImag.append(str(allW.name))
        else:
            wp_hiImag.append(str(allW.name))
        count=count+1

    print wp_loImag,'\n\n\n'
    print wp_hiImag,'\n\n\n'

    # copy the word pool to the sessions dir
    shutil.copy(config.wp,        exp.session.fullPath())    
    shutil.copy(config.wp_imagery,exp.session.fullPath())    
    sessionList_loImag = []
    sessionList_hiImag = []
    sessionLL   = []

    for sessionNum in xrange(config.numSessions):
        # Set the session so that we know the correct directory
        exp.setSession(sessionNum)
        
        # Get the session-specific config
        sessionconfig = config.sequence(sessionNum)
        
        # Shuffle the word pool
        random.shuffle(wp_loImag)
        random.shuffle(wp_hiImag)
        sliceStart = 0
        trialList_loImag  = []
        trialList_hiImag  = []

        # make the list length vector
        LL_vect = []        
        for LL_vect_count in xrange(sessionconfig.numTrials):
            LL_vect.append(sessionconfig.listLen[LL_vect_count%len(sessionconfig.listLen)])
        random.shuffle(LL_vect)
        
        for trialNum in xrange(sessionconfig.numTrials):
            # Get the trial specific config
            trialconfig = sessionconfig.sequence(trialNum)
            
            # slice the list
	    sliceEnd = sliceStart+LL_vect[trialNum]
	    if sliceEnd>len(wp_loImag) or sliceEnd>len(wp_hiImag):
		# reshuffle the word pool and reset the indices
		random.shuffle(wp_loImag)
                random.shuffle(wp_hiImag)
		sliceStart = 0
		sliceEnd = sliceStart+LL_vect[trialNum]
                print "RAN OUT OF WORD AT LIST %d\n\n"%trialNum
            currentList_loImag = wp_loImag[sliceStart:sliceEnd]
            currentList_hiImag = wp_hiImag[sliceStart:sliceEnd]
            trialList_loImag.append(currentList_loImag)
            trialList_hiImag.append(currentList_hiImag)

            # update the slice start
            sliceStart += LL_vect[trialNum]

            # write the list out to file
            listFile_loImag = exp.session.createFile("%d_loImag.txt"  % trialNum)
            listFile_hiImag = exp.session.createFile("%d_hiImag.txt"  % trialNum)

            for eachWord in currentList_loImag:
                listFile_loImag.write("%s\n" % eachWord)
            listFile_loImag.close()

            for eachWord in currentList_hiImag:
                listFile_hiImag.write("%s\n" % eachWord)
            listFile_hiImag.close()

        # append the trials to the session
        sessionList_loImag.append(trialList_loImag)
        sessionList_hiImag.append(trialList_hiImag)
        sessionLL.append(LL_vect)

    # save the prepared data
    exp.saveState(state,
                  trialNum=0,
                  sessionList_loImag=sessionList_loImag,
                  sessionList_hiImag=sessionList_hiImag,
                  sessionLL=sessionLL,
                  sessionNum=0)

def run(exp, config):
    """
    Run a session of free recall.
    """

    # set priority
    if config.doRealtime:
        setRealtime(config.rtPeriod,config.rtComputation,config.rtConstraint)
    
    # verify that we have all the files
    verifyFiles(config)
    
    # get the state
    state = exp.restoreState()
    
    # set up the session
    # have we run all the sessions
    if state.sessionNum >= len(state.sessionList_hiImag):
        print "No more sessions!"
        return

    # set the session number
    exp.setSession(state.sessionNum)
    
    # get session specific config
    sessionconfig = config.sequence(state.sessionNum)

    # create tracks...
    video = VideoTrack("video")
    audio = AudioTrack("audio")
    keyboard = KeyTrack("keyboard")
    log = LogTrack("session")
    mathlog = LogTrack("math")
    if config.isPrepSession:
        eeg = EEGTrack("eeg") # by doing this, we prevent startLogging()
    else:
        eeg = EEGTrack("eeg", autoStart=False) # by doing this, we prevent startLogging()
        eeg.startService() #  from being called, which prevents the intermittent 
        # callback-triggered pulses from being sent
        eeg.logall = True    
    # set the default font
    setDefaultFont(Font(sessionconfig.defaultFont))
    
    # get a presentation clock
    clock = PresentationClock()

    # log start
    timestamp = clock.get()
    log.logMessage('SESS_START\t%d' % (state.sessionNum + 1), timestamp)
    bc = ButtonChooser(Key(sessionconfig.terminateRecallKey))

    # do instructions on first trial of each session
    if state.trialNum == 0:
        # do miketest
        video.clear("black")
        soundgood = micTest(2000,1.0)
        if not soundgood:
            # quit
            return
	
	# Show main instructions
	video.clear("black")
        instruct(open(sessionconfig.instruct,'r').read())

    # present each trial in the session
    didThisOnce = False
    while state.trialNum < len(state.sessionList_hiImag[state.sessionNum]):
        # load trial specific config
        trialconfig = sessionconfig.sequence(state.trialNum)

	# create the beeps
	startBeep = Beep(trialconfig.startBeepFreq,
			 trialconfig.startBeepDur,
			 trialconfig.startBeepRiseFall)
	stopBeep = Beep(trialconfig.stopBeepFreq,
			trialconfig.stopBeepDur,
			trialconfig.stopBeepRiseFall)
        
        # show the current trial and wait for keypress
        while True:
            video.clear("black")                    
            readyText = video.showCentered(Text("Press any key for Trial #%d" % (state.trialNum + 1)))    
            video.updateScreen(clock)

            # press 'T' to gain control
            settings_bc = keyboard.keyChooser().waitWithTime(minDuration=None, 
                                                             maxDuration=None, clock=clock)
            if (settings_bc[0].name == 'T' or not didThisOnce):
                didThisOnce = True
                video.unshow(readyText)
                video.updateScreen(clock)
                while True:
                    doImag = textInput('HI or LO Imagery words? [H]/[L]:',video,keyboard,clock)
                    if doImag == 'H':
                        dohiImageryWords = True
                        logImagText='HiImag'
                        break
                    elif doImag == 'L':
                        dohiImageryWords = False
                        logImagText='LoImag'    
                        break
                if not trialconfig.isPrepSession:
                    while True:
                        doStim = textInput('Stimulation? [Y]/[N]:',video,keyboard,clock)
                        if doStim == 'Y':                
                            thisIsAStimTrial = True
                            elec        = textInput('Electrode Number:',video,keyboard,clock)
                            mAmps       = textInput('mAmps:',video,keyboard,clock)
                            freqHZ      = textInput('Frequency (Hz):',video,keyboard,clock)
                            duration    = textInput('Stim duration (sec):',video,keyboard,clock)
                            while True:
                                StimTimeIn  = textInput('Encoding or Recalled [E]/[R]:',video,keyboard,clock)
                                if StimTimeIn == 'E':
                                    doEncodingStim = True
                                    break
                                elif StimTimeIn == 'R':
                                    doEncodingStim = False
                                    break
                            logStimText = "STIM\t%s\t%s\t%s\t%s\t%s_stim"%(elec,mAmps,freqHZ,duration,StimTimeIn)
                            break
                        elif doStim == 'N':
                            logStimText      = 'NO_STIM'
                            thisIsAStimTrial = False
                            break
                else:
                    thisIsAStimTrial = False
                    logStimText      = 'NO_STIM'                
            else:
                break
                
        # remove wait for any key and start list
        video.unshow(readyText)
        video.updateScreen(clock) 
        timestamp = clock.get()
        log.logMessage('TRIAL\t%d\t%s\t%s' % (state.trialNum + 1,logImagText,logStimText),timestamp)

        # which list (hi or low imagery?)
        if dohiImageryWords:
            thisSessList=state.sessionList_hiImag[state.sessionNum][state.trialNum]
        else:
            thisSessList=state.sessionList_loImag[state.sessionNum][state.trialNum]
        
        # write the list out to file
        listFile = exp.session.createFile("%d.lst"  % state.trialNum)
        for eachWord in thisSessList:
            listFile.write("%s\n" % eachWord)
        listFile.close()

        # display the "cross-hairs"
        timestamp = flashStimulus(Text(trialconfig.orientText, size = trialconfig.wordHeight),
                                      clk=clock,duration=trialconfig.orient_Duration)
        clock.delay(trialconfig.PauseBeforeWords, trialconfig.JitterBeforeWords)
        log.logMessage('ORIENT', timestamp)

        # record at the start fo the study period
        (rec,timestamp) = audio.startRecording("studyperiod_"+str(state.trialNum),t=clock.get())
        log.logMessage('REC_STUDY', timestamp)
        
        # start STIM HERE
        if thisIsAStimTrial and doEncodingStim:
            ts_stim = eeg.timedPulse(trialconfig.pulsestim_length,
                                     trialconfig.pulsestim_label, 
                                     trialconfig.pulsestim_signal, 
                                     clk=clock)
            print 'I just sent a stim pulse during encoding'
            log.logMessage('PULSE_SENT_ENCODING', ts_stim)
            clock.delay(trialconfig.PauseBtwnStimAndFirstWordStart)            
            

        # start them words
        for n, eachWord in enumerate(thisSessList):
            wordconfig = trialconfig.sequence(n)   
            timestamp = flashStimulus(Text(eachWord, size = trialconfig.wordHeight),
                                      clk=clock,duration=wordconfig.word_Duration)
            clock.delay(wordconfig.word_ISI, wordconfig.word_Jitter)            
            log.logMessage('WORD\t%s\t%d' % (eachWord, n), timestamp)
        (r,stoptime) = audio.stopRecording(t = clock.get())

        # do the math distract
        if trialconfig.doMathDistract:
            mathDistract(clk = clock,
                         mathlog = mathlog,
                         numVars = trialconfig.MATH_numVars,
                         maxProbs = trialconfig.MATH_maxProbs,
                         plusAndMinus = trialconfig.MATH_plusAndMinus,
                         minDuration = trialconfig.MATH_minDuration,
                         textSize = trialconfig.MATH_textSize)

        # Pause before recall
        clock.delay(trialconfig.PauseBeforeRecall, trialconfig.JitterBeforeRecall)
        
        # show the recall start indicator
        startText = video.showCentered(Text(trialconfig.recallStartText, 
                                                size = trialconfig.wordHeight))
        video.updateScreen(clock)
        startBeep.present(clock)
	
        # hide rec start text
        video.unshow(startText)
        video.updateScreen(clock)

        # Record responses
        (rec,timestamp) = audio.startRecording(str(state.trialNum), t=clock)
        log.logMessage('REC_START' % (),timestamp)
        
        # start STIM HERE
        if thisIsAStimTrial and not doEncodingStim:
            ts_stim = eeg.timedPulse(trialconfig.pulsestim_length,
                                     trialconfig.pulsestim_label, 
                                     trialconfig.pulsestim_signal, 
                                     clk=clock)
            print 'I just sent a stim pulse during retrieval'
            log.logMessage('PULSE_SENT_RETRIEVAL', ts_stim)
            clock.delay(trialconfig.PauseBtwnStimAndFirstWordStart)            

        
        # wait for key with max dur=recallDuration
        button,bc_time = bc.waitWithTime(maxDuration=trialconfig.recallDuration, clock=clock)

        # stop recording
        (r,stoptime) = audio.stopRecording(t = clock.get())        
        #stopBeep.present(clock)
        log.logMessage('REC_STOP' % (),stoptime)
              
        # update the trialnum state
        state.trialNum += 1

        # save the state after each trial
        exp.saveState(state)
        
    # save the state when the session is finished
    exp.saveState(state, trialNum = 0, sessionNum = state.sessionNum + 1)

    # Done
    timestamp = waitForAnyKey(clock,Text("Thank you!\nYou have completed the session."))
    log.logMessage('SESS_END',timestamp)

    # Catch up
    clock.wait()

# only do this if the experiment is run as a stand-alone program (not imported as a library)
if __name__ == "__main__":
    # make sure we have the min pyepl version
    checkVersion(MIN_PYEPL_VERSION)
    
    # start PyEPL, parse command line options, and do subject housekeeping
    exp = Experiment(use_eeg=False) # we need to configure EEGTrack manually
    exp.parseArgs()
    exp.setup()
    
    # allow users to break out of the experiment with escape-F1 (the default key combo)
    exp.setBreak()
    
    # get subj. config
    config = exp.getConfig()

    # if there was no saved state, run the prepare function
    if not exp.restoreState():
        print "\n\n\n\t***Running Prep***\n\n\n"
        prepare(exp, config)		
        print "\n\n\n\t***Done prep***\n\n\n"

    # now run the subject
    run(exp, config)
    
