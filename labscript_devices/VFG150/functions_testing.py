# VFG functions
import numpy as np
import math
import random
import ctypes as ct
# from labscript import LabscriptError

from python_toolbox import caching


def freq(frequency):
    ''' freq(F) Generates a VFG-150 sequence that sets the frequency to F
        (given in Hz). '''
    ticks = round(frequency * (2**32) / 200e6) # external clock
    ret = [ (7*16 + 4), (math.floor(ticks / 2**24))%(2**8),\
        (math.floor(ticks / 2**16))%(2**8),\
        (math.floor(ticks / 2**8))%(2**8),\
        (ticks)%(2**8) ]
    return ret


def ampl(amplitude):
    ''' ampl(A) Generates a VFG-150 sequence that sets the amplitude to A. The
        range of A is from 0 (no output) to 1 (full output). The actual amplitude
        is frequency dependent and falls from -4 dBm at a frequency of 1 MHz to
        about -7 dBm at a frequency of 150 MHz.'''
    ticks = round(amplitude * (2**16 - 1))
    ret = [ (11*16 + 10),\
        math.floor(ticks / 2**8),\
        (ticks)%(2**8) ]
    
    # uint8 format
    ret = [max(min(x, 255), 0) for x in ret]
    return ret


def ampl_slope(slope):
    ''' AMPL_SLOPE(S) Generates a VFG-150 sequence that sets the amplitude slope
        to S. The slope is given in full scale per second: an amplitude slope
        of 1 raises the amplitude in 1 second from 0 to 1.
        This command generates an error if the amplitude required is too
        steep (greater than 1 per 5 ns).
        modified by Roman, 6/1/2011:
        - the maximum slope is 1 per 10ns, not per 5ns!
        - the range checking was not done properly '''
    # ticks = round(slope * (2**24-1) / 200E6)
    ticks = round(slope * 2**24 / 200E6)   # changed by Roman, 6/1/2011 (makes more sense this way)
    if ticks < -2**23:
        print('ampl_slope: amplitude too steep (set to -1e8)')
        ticks = -2**23
    
    if ticks > 2**23-1:
        print('ampl_slope: amplitude too steep (set to +1e8)')
        ticks = 2**23-1;
    
    if ticks < 0:
        ticks = ticks + 2**24

    ret = [ (14*16 + 12),\
        (math.floor(ticks / 2**16))%(2**8),\
        (math.floor(ticks / 2**8))%(2**8),\
        (ticks)%(2**8) ]
    return ret


def aux(port, state):
    '''aux(P, S) Generates a VFG-150 sequence that sets the auxiliary outs
       of the VFG-150.
       The VFG-150 has four auxiliary outputs (AUX_OUT1 to AUX_OUT4). This
       command accepts a bit mask P and a bit value S. To change the state
       of an output the bit at the corresponding position in the bit mask has
       to be set. If it is set, then the new state of the output is
       determined by the value of S at the corresponding bit position.
          
       Use S = 1 for 'ON'
           S = 0 for 'OFF'
       For example to set the first auxiliary output:
         aux_out(1,1)
       To reset the first auxiliary output:
         aux_out(1,0) '''
    
    port = 2**(port-1)
    state = port * state
    ret = [ 255, (port)%(16) * 16 + (state)%(16) ]
    return ret


def interval(time):
    '''interval(T) generates a VFG-150 sequence that inserts a pause of length
    T. The duration T is given in seconds.
    Pauses are used to structure the sequence that the VFG-150 shall generate.
    All commands without a pause between them will be executed at once.
 
    The pause generated with this command will have a dither of the order
    of 5 ns. This is done to prevent a systematic bias in cases where the
    requested delays round systematically to a multiple of 5 ns. '''
    # We don't want random numbers
    # ticks = round(time * 200E6 + random.random() - 0.5) # external clock
    # ticks = round(time * 192E6 + random.random() - 0.5) # internal clock
    ticks = round(time * 200E6) # external clock
    
    if ticks >= 2**32:
        print('interval: duration too long (21.474s maximum)')
	
    ret = [ (16*3 + 0),\
        (math.floor(ticks/2**24))%(2**8),\
        (math.floor(ticks/2**16))%(2**8),\
        (math.floor(ticks/2**8))%(2**8),\
        (ticks)%(2**8)]
    
    # Replacing zero or negative delays by NOPs
    if ticks <= 0:
        ret = [1,1,1,1,1]
    
    return ret


def phase(phi):
    '''phase(phi) Generates a VFG-150 sequence that sets the phase to PHI. PHI
       is given in radians.
       The notion of phase with this command is that of an absolute
       offset phase. So the phase jump that the generated waveform experiences
       is given by the difference between the current phase value and the new
       phase value that is set with this command. The default value for the
       offset phase at the beginning of each sequence is zero. '''
       
    ticks = round(phi * 2**16 / 2 / np.pi)
    ret = [ (9*16 + 8), (math.floor(ticks/2**8))%(2**8), (ticks)%(2**8) ]
    return ret


def reset_timebase():
    ''' reset_timebase() Generates a VFG-150 sequence that resets the time base
        register. This is the register which is incremented at each time step to
        track the absolute phase of the generated signal.'''
        
    ret = [ 4 ]
    return ret

def set_phase_continuous():
    ''' set_phase_continuous() generates a VFG-150 sequence that sets the
        synthesizer to phase continuous frequency switching as opposed to
        phase coherent switching. The default at the beginning of each
        sequence is phase continuous switching. '''
    
    ret = [ 2 ]
    return ret


def trigger():
    ''' trigger generates a VFG-150 sequence that inserts a "wait for trigger"
        event at the corresponding position in the sequence. A "wait for
        trigger" events stops the processing of the sequence until a trigger
        is registered at the trigger input of the VFG-150. The trigger input
        is found at the front panel and labeled "Trigger In".'''
        
    ret = [ 0, 0]
    return ret


def newInterval(t0, ph0, a0, a1, a2, err):
    ''' given a starting time and a starting phase, find the optimal end point of
        the next interval such that the phase error is within the prescribed
        bound '''
    if a2 > 0:
        f = a1 + math.sqrt(2*a2*(a0 + err - ph0))
    else:
        f = a1 - math.sqrt(2*a2*(a0 - err - ph0))
        
    t1 = t0 + 2*math.sqrt(err/abs(a2))+(f - a1)/a2
    ph1 = ph0 + f*(t1 - t0)
    j = np.array([[t1, ph1]])
    return j



def generalRamp(phaseFunc, phaseFuncD1, phaseFuncD2, amplitudeFunc,\
    pulseDuration, maxPhaseError):

    # TODO: uncomment
    # if len(sys.arg) != 6:
    #     raise LabscriptError('Incorrect number of input arguments');
    
    # if pulseDuration <= 0:
    #     raise LabscriptError('Pulse duration must be positive');
    
    # if maxPhaseError <= 0:
    #     raise LabscriptError('Maximum phase error must be positive')

    
    # build the list of times and accumulated phases
    L = np.array([[0,0]])
    while L[-1,0] < pulseDuration:
        add = newInterval(L[-1,0],L[-1,1],\
        phaseFunc(L[-1,0]),\
        phaseFuncD1(L[-1,0]),\
        phaseFuncD2(L[-1,0]),\
        maxPhaseError)
        L = np.concatenate((L,add))
    
    L[-1,0] = pulseDuration
    L[-1,1] = phaseFunc(pulseDuration)
    

    # add amplitude information to the list (third column)
    L = np.c_[L,amplitudeFunc(0) * np.ones(len(L[:,0]))]

    # for the beginning of each interval, make a list of frequency,
    # start amplitude, amplitude change rate, and duration
    intervals = L[1:,0] - L[0:-1,0]
    frequencies = ((L[1:,1] - L[0:-1,1])/intervals)/(2*np.pi);
    startamplitudes = L[0:-1,2]
    amplituderates = (L[1:,2]-L[0:-1,2])/intervals
    sequenceVec = np.c_[frequencies, startamplitudes, amplituderates, intervals]

    # make a sequence for the VFG-150
    sequence = []
    for i in range(len(intervals)):
        sequence = sequence + freq(frequencies[i]) +\
            ampl(startamplitudes[i]) +\
            ampl_slope(amplituderates[i]) +\
            interval(intervals[i])

    # at the end the sequence finishes on a constant note
    sequence = sequence + freq(phaseFuncD1(pulseDuration)/(2*np.pi)) +\
        ampl(amplitudeFunc(pulseDuration)) +\
        ampl_slope(0)
        
    sequenceVec = np.concatenate((sequenceVec, np.array([[phaseFuncD1(pulseDuration)/(2*np.pi),\
        amplitudeFunc(pulseDuration), 0, np.inf]])))

    return sequence, sequenceVec


# @caching.cache()
def exponentialFrequencyRamp(initialFrequency, finalFrequency,\
    dropTime, pulseDuration, maxPhaseError = math.pi/2, amplitude=1):
    print('calculating')
    # if pulseDuration <= 0:
    #     raise LabscriptError('Pulse duration must be positive')
    
    # if maxPhaseError <= 0:
    #     raise LabscriptError('Maximum phase error must be positive')

    # the phase function and its first and second derivatives:
    e1 = math.exp(pulseDuration/dropTime)
    e2 = math.exp(-pulseDuration/dropTime)
    h1 = 2*np.pi*(initialFrequency-finalFrequency)/(e1-1)
    h2 = 2*np.pi*(initialFrequency-finalFrequency)/(1-e2)
    phaseFunc = lambda t: (2*np.pi*finalFrequency-h1)*t \
        + dropTime*h2*(1-math.exp(-t/dropTime))
    phaseFuncD1 = lambda t: 2*np.pi*finalFrequency \
        + h2*(math.exp(-t/dropTime)-e2)
    phaseFuncD2 = lambda t: -h2*math.exp(-t/dropTime)/dropTime
    
    # the amplitude function:
    amplitudeFunc = lambda t: amplitude
    
    sequence, sequenceVec = generalRamp(phaseFunc, phaseFuncD1,\
        phaseFuncD2, amplitudeFunc,\
        pulseDuration, maxPhaseError)
        
    # sesequenceVec only used for debugging
    
    # remove decimals
    sequence = [int(round(i)) for i in sequence]
    return sequence














