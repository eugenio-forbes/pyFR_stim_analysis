jfburke 02-2013:
I modified pyFR_stim for a few reasons:
(1) stimulation needs to happen for the whole list, not indiviual
    words on the list.  Therefore, lists needed to be shorter.
    Therefore, we needed a distraction task.  
(2) Dr. Sperling wants manual control of stimulation.

Instructions:
(1) run pyFR_stim with the config set to baseline
    - in the config, set "isPrepSession = TRUE"
    - python pyFR_stim.py -s subj_ID    
(2) get psychometric function and select the list length
    - run analysis/getAccuracybyLL_stimPreSession.m
(3) run the stimulation version with the selected list length.
    - in the config, set "isPrepSession = FALSE"
    - select list length in the config

Make events:
(1) for no stim.... extractPYFRevents_noSTIM.m
(2) for stim....... extractPYFRevents_STIM.m

 
To-do:
(1) make tables... right now this just spits out a .txt file from the 
    pyepl code, which is very low-level.  Ideally, the program should spit 
    out a nice latex formatted table so we can keep up in the patient's room. 
(2) cos theta
