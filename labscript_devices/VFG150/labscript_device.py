
import math
import numpy as np
from labscript import Device, LabscriptError, set_passed_properties

from labscript_devices.VFG150.functions import VFG_cached_class
from python_toolbox import caching
import functools

# functions for programming vfg
def f_freq(frequency):
    ''' freq(F) Generates a VFG-150 sequence that sets the frequency to F
        (given in Hz). '''
    ticks = round(frequency * (2**32) / 200e6) # external clock
    ret = [ (7*16 + 4), (math.floor(ticks / 2**24))%(2**8),\
        (math.floor(ticks / 2**16))%(2**8),\
        (math.floor(ticks / 2**8))%(2**8),\
        (ticks)%(2**8) ]
        
    return ret

def f_freq_dense(frequencyold,frequencynew):
    ''' f_freq_dense(Fo, Fn) overwrites an old frequency (given in Hz) with a new,
        only setting the bits thath change (dense programming, see VFG150 manual pp. 47)'''
    ro = np.array(f_freq(frequencyold))
    rn = np.array(f_freq(frequencynew))
    
    ind = np.argwhere(np.subtract(ro,rn))
    
    if len(ind) == 0:
        ret=[68, rn[4]]
    else:
        ret = [(7+(1-np.amin(ind)))*16 + 4+(4-np.amax(ind))]
        for i in range(np.amin(ind),np.amax(ind)+1):
            ret.append(rn[i])
    return ret 

def f_ampl(amplitude):
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

def f_ampl_slope(slope):
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

# One could fix the problem with inverted outputs with few lines of: if port ==1: port = 4 etc...
def f_aux(port, state):
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

def f_interval(time):
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
        ret = [1, 1, 1, 1, 1]
    
    return ret

### This does not work yet (no idea whz)... (tried to call it in generalramp)
def f_interval_dense(timeold, timenew):
    '''f_interval_dense(To, Tn) inserts a pause of length Tn
       by only setting bits that change (dense programming, see VFG150 manual pp. 47)  '''
    if timeold == 0:
        return f_interval(timenew)
    else:
        ro = np.array(f_interval(timeold))
        rn = np.array(f_interval(timenew))
        
        ind = np.argwhere(np.subtract(ro,rn))
        
        if len(ind) == 0:
            ret=[17, rn[3]]
        else:
            ret = [(3 + (1-np.amin(ind)))*16 + 0 + (4-np.amax(ind))]
            for i in range(np.amin(ind), np.amax(ind) + 1):
                ret.append(rn[i])
        return ret 

def f_phase(phi):
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

def f_reset_timebase():
    ''' reset_timebase() Generates a VFG-150 sequence that resets the time base
        register. This is the register which is incremented at each time step to
        track the absolute phase of the generated signal.'''
        
    ret = [ 4 ]
    
    return ret

def f_set_phase_continuous():
    ''' set_phase_continuous() generates a VFG-150 sequence that sets the
        synthesizer to phase continuous frequency switching as opposed to
        phase coherent switching. The default at the beginning of each
        sequence is phase continuous switching. '''
    
    ret = [ 2 ]
    
    return ret

def f_set_phase_coherent():
    ''' set_phase_coherent() generates a VFG-150 sequence that sets the
        synthesizer to phase coherent frequency switching as opposed to
        phase continuous switching. The default at the beginning of each
        sequence is phase continuous switching. '''
    
    ret = [ 3 ]
    
    return ret

def f_trigger():
    ''' trigger() generates a VFG-150 sequence that inserts a "wait for trigger"
        event at the corresponding position in the sequence. A "wait for
        trigger" events stops the processing of the sequence until a trigger
        is registered at the trigger input of the VFG-150. The trigger input
        is found at the front panel and labeled "Trigger In".'''
        
    ret = [48, 0, 0, 0, 0]
    
    return ret

def f_newInterval(t0, ph0, a0, a1, a2, err):
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


def f_generalRamp(phaseFunc, phaseFuncD1, phaseFuncD2, amplitudeFunc,\
    pulseDuration, maxPhaseError):
    
    if pulseDuration <= 0:
        raise LabscriptError('Pulse duration must be positive');
    
    if maxPhaseError <= 0:
        raise LabscriptError('Maximum phase error must be positive')

    # build the list of times and accumulated phases
    L = np.array([[0,0]])
    while L[-1,0] < pulseDuration:
        add = f_newInterval(L[-1,0],L[-1,1],\
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
    frequencies = ((L[1:,1] - L[0:-1,1])/intervals)/(2*np.pi)
    startamplitudes = L[0:-1,2]
    amplituderates = (L[1:,2]-L[0:-1,2])/intervals
    
    sequenceVec = np.c_[frequencies, startamplitudes, amplituderates, intervals]

    # make a sequence for the VFG-150
    sequence = []
    for i in range(len(intervals)):
        sequence = sequence + f_freq(frequencies[i]) +\
            f_ampl(startamplitudes[i]) +\
            f_ampl_slope(amplituderates[i]) +\
            f_interval(intervals[i])

    # at the end the sequence finishes on a constant note
    sequence = sequence + f_freq(phaseFuncD1(pulseDuration)/(2*np.pi)) +\
        f_ampl(amplitudeFunc(pulseDuration)) +\
        f_ampl_slope(0)
        
    sequenceVec = np.concatenate((sequenceVec, np.array([[phaseFuncD1(pulseDuration)/(2*np.pi),\
        amplitudeFunc(pulseDuration), 0, np.inf]])))

    return sequence, sequenceVec


def f_generalRamp_no_ampl_change(phaseFunc, phaseFuncD1, phaseFuncD2, amplitudeFunc,\
    pulseDuration, maxPhaseError):
    # to be used for frequency ramp with const ampl (more efficient)
    
    if pulseDuration <= 0:
        raise LabscriptError('Pulse duration must be positive');
    
    if maxPhaseError <= 0:
        raise LabscriptError('Maximum phase error must be positive')

    # build the list of times and accumulated phases
    L = np.array([[0,0]])
    while L[-1,0] < pulseDuration:
        add = f_newInterval(L[-1,0],L[-1,1],\
        phaseFunc(L[-1,0]),\
        phaseFuncD1(L[-1,0]),\
        phaseFuncD2(L[-1,0]),\
        maxPhaseError)
        L = np.concatenate((L,add))
    
    L[-1,0] = pulseDuration
    L[-1,1] = phaseFunc(pulseDuration)

    # for the beginning of each interval, make a list of frequency,
    # start amplitude, amplitude change rate, and duration
    intervals = L[1:,0] - L[0:-1,0]
    frequencies = ((L[1:,1] - L[0:-1,1])/intervals)/(2*np.pi)
    
    # make a sequence for the VFG-150
    sequence = []
    sequence += f_ampl_slope(0) + f_ampl(amplitudeFunc(0))
    
    for i in range(len(intervals)):
        sequence += f_freq(frequencies[i]) +\
                    f_interval(intervals[i])

    # at the end the sequence finishes on a constant note
    sequence += f_freq(phaseFuncD1(pulseDuration)/(2*np.pi)) +\
        f_ampl(amplitudeFunc(pulseDuration)) +\
        f_ampl_slope(0)
    
    # for debugging
    sequenceVec = np.c_[frequencies, intervals]
    sequenceVec = np.concatenate((sequenceVec, np.array([[phaseFuncD1(pulseDuration)/(2*np.pi), 0]])))

    return sequence, sequenceVec


def f_generalRamp_no_ampl_change_dense(phaseFunc, phaseFuncD1, phaseFuncD2, amplitudeFunc,\
    pulseDuration, maxPhaseError):
    # to be used for frequency ramp with const ampl (more efficient)
    
    if pulseDuration <= 0:
        raise LabscriptError('Pulse duration must be positive');
    
    if maxPhaseError <= 0:
        raise LabscriptError('Maximum phase error must be positive')

    # build the list of times and accumulated phases
    L = np.array([[0,0]])
    while L[-1,0] < pulseDuration:
        add = f_newInterval(L[-1,0],L[-1,1],\
        phaseFunc(L[-1,0]),\
        phaseFuncD1(L[-1,0]),\
        phaseFuncD2(L[-1,0]),\
        maxPhaseError)
        L = np.concatenate((L,add))
    
    L[-1,0] = pulseDuration
    L[-1,1] = phaseFunc(pulseDuration)

    # for the beginning of each interval, make a list of frequency,
    # start amplitude, amplitude change rate, and duration
    intervals = L[1:,0] - L[0:-1,0]
    frequencies = ((L[1:,1] - L[0:-1,1])/intervals)/(2*np.pi)
    
    # make a sequence for the VFG-150
    sequence = []
    sequence += f_ampl_slope(0) + f_ampl(amplitudeFunc(0))
    
    sequence +=f_freq(frequencies[0])+f_interval(intervals[0])
     
    for i in range(1,len(intervals)):
        sequence += f_freq_dense(frequencies[i-1],frequencies[i]) +\
                    f_interval(intervals[i])
                    # f_interval_dense(intervals[i-1],intervals[i])
                    
    # at the end the sequence finishes on a constant note
    sequence += f_freq(phaseFuncD1(pulseDuration)/(2*np.pi)) +\
        f_ampl(amplitudeFunc(pulseDuration)) +\
        f_ampl_slope(0)
    
    # for debugging
    sequenceVec = np.c_[frequencies, intervals]
    sequenceVec = np.concatenate((sequenceVec, np.array([[phaseFuncD1(pulseDuration)/(2*np.pi), 0]])))

    return sequence, sequenceVec

# most dense (until now): if we write freq byte 4, directly add interval information
def f_generalRamp_dense(phaseFunc, phaseFuncD1, phaseFuncD2, amplitudeFunc,\
    pulseDuration, maxPhaseError):
    # to be used for frequency ramp with const ampl (more efficient)
    
    if pulseDuration <= 0:
        raise LabscriptError('Pulse duration must be positive');
    
    if maxPhaseError <= 0:
        raise LabscriptError('Maximum phase error must be positive')

    # build the list of times and accumulated phases
    L = np.array([[0,0]])
    while L[-1,0] < pulseDuration:
        add = f_newInterval(L[-1,0],L[-1,1],\
        phaseFunc(L[-1,0]),\
        phaseFuncD1(L[-1,0]),\
        phaseFuncD2(L[-1,0]),\
        maxPhaseError)
        L = np.concatenate((L,add))
    
    L[-1,0] = pulseDuration
    L[-1,1] = phaseFunc(pulseDuration)

    # for the beginning of each interval, make a list of frequency,
    # start amplitude, amplitude change rate, and duration
    intervals = L[1:,0] - L[0:-1,0]
    frequencies = ((L[1:,1] - L[0:-1,1])/intervals)/(2*np.pi)
    
    # make a sequence for the VFG-150
    sequence = []
    # set slope and amplitude at beginning of ramp
    sequence += f_ampl_slope(0) + f_ampl(amplitudeFunc(0))
    # set first freq and interval
    sequence += f_freq(frequencies[0]) + f_interval(intervals[0])
     
    for i in range(1,len(intervals)):
        # if we write the last bit of freq, attach the interval information directly
        if f_freq_dense(frequencies[i-1], frequencies[i])[0]%16 == 4:
            sequence += [f_freq_dense(frequencies[i-1],frequencies[i])[0]-4] + f_freq_dense(frequencies[i-1],frequencies[i])[1:]
            sequence += f_interval(intervals[i])[1:] # no need for first bit!
            
        # if not, write them seperatly
        else:
            sequence += f_freq_dense(frequencies[i-1], frequencies[i])
            sequence += f_interval(intervals[i])
                    
    # at the end the sequence finishes on a constant note
    sequence += f_freq(phaseFuncD1(pulseDuration)/(2*np.pi)) +\
        f_ampl(amplitudeFunc(pulseDuration)) +\
        f_ampl_slope(0)
    
    # for debugging
    sequenceVec = np.c_[frequencies, intervals]
    sequenceVec = np.concatenate((sequenceVec, np.array([[phaseFuncD1(pulseDuration)/(2*np.pi), 0]])))

    return sequence, sequenceVec

@functools.lru_cache(maxsize=10)
def f_exponentialFrequencyRamp(initialFrequency, finalFrequency, dropTime, pulseDuration, maxPhaseError = np.pi/2, amplitude = 1):
    print('Calculating...')
    if pulseDuration <= 0:
        raise LabscriptError('Pulse duration must be positive')
    
    if maxPhaseError <= 0:
        raise LabscriptError('Maximum phase error must be positive')

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
    # could pass this as a simple arg as well..
    amplitudeFunc = lambda t: amplitude
    
    # sesequenceVec only used for debugging
    sequence, sequenceVec = f_generalRamp_dense(phaseFunc, phaseFuncD1,
        phaseFuncD2, amplitudeFunc,
        pulseDuration, maxPhaseError)

    # remove decimals
    sequence = [int(round(i)) for i in sequence]
    
    return sequence

def f_vfg_pulse(old_pulse_dict, new_pulse_dict):
    S = []
    for key in old_pulse_dict:
        if old_pulse_dict[key] != new_pulse_dict[key]:
            if key == 'aux1':
                S += f_aux(1, new_pulse_dict[key])
            elif key == 'aux2':
                S += f_aux(2, new_pulse_dict[key])
            elif key == 'aux3':
                S += f_aux(3, new_pulse_dict[key])
            elif key == 'aux4':
                S += f_aux(4, new_pulse_dict[key])   
            elif key == 'vfg_phase':
                S += f_phase(new_pulse_dict[key])
            elif key == 'vfg_frequency':
                S += f_freq(new_pulse_dict[key])
            elif key == 'vfg_amplitude':
                S += f_ampl(new_pulse_dict[key])
            elif key == 'vfg_slope':
                S += f_ampl_slope(new_pulse_dict[key])
    # always an interval
    S += f_interval(new_pulse_dict['dt'])
    
    return S

def f_first_vfg_pulse(pulse_dict):
    S = []
    S += f_aux(1, pulse_dict['aux1'])
    S += f_aux(2, pulse_dict['aux2'])
    S += f_aux(3, pulse_dict['aux3'])
    S += f_aux(4, pulse_dict['aux4'])
    S += f_phase(pulse_dict['vfg_phase'])
    S += f_freq(pulse_dict['vfg_frequency'])
    S += f_ampl(pulse_dict['vfg_amplitude'])
    S += f_ampl_slope(pulse_dict['vfg_slope'])
    S += f_interval(pulse_dict['dt'])
    
    return S


# used to generate the sequence S
# input is a list of pulse_dicts, created from pulse class in bec_utils
def f_pulse_seq(pulse_list):
    # always send
    S = f_set_phase_continuous() + f_trigger()
    
    # create first pulse
    S += f_first_vfg_pulse(pulse_list[0])
    
    # update and create further pulses
    if len(pulse_list) > 1:
        for i in range(len(pulse_list)-1):
            # create new pulses
            S += f_vfg_pulse(pulse_list[i], pulse_list[i+1])
    
    S += f_ampl(0) + f_ampl_slope(0)
    return S


class VFG150(Device):
    """A labscript_device for the VFG150 such that it can be controlled via python.
    """
    description = 'VFG150 via python'
    
    def __init__(self, name):
        Device.__init__(self, name, None, None)
        self.name = name
        self.BLACS_connection = 'VFG150' # This is needed to create the blacs_tab and will appear as [conn: self.BLACS_connection]
        self.S = []
        
    def generate_code(self, hdf5_file):
        # Code wich is executed by runmanager to create/store data in h5 file.
        # create group with the name of device
        grp = hdf5_file.create_group('/devices/'+self.name)
        # Write data to new dataset
        if self.S != []:
            grp.create_dataset('S', data = self.S)

    def freq(self,frequency):
        self.S +=f_freq(frequency)
    
    def ampl(self,amplitude):
        self.S +=f_ampl(amplitude)
        
    def ampl_slope(self,slope):
        self.S +=f_ampl_slope(slope)
        
    def aux(self,port, state):
        self.S +=f_aux(port, state)
    
    def interval(self,time):
        self.S +=f_interval(time)
        
    def phase(self,phi):
        self.S +=f_phase(phi)
    
    def reset_timebase(self):
        self.S +=f_reset_timebase()
        
    def set_phase_continuous(self):
        self.S +=f_set_phase_continuous()
        
    def set_phase_coherent(self):
        self.S +=f_set_phase_coherent()
    
    def trigger(self):
        self.S +=f_trigger()
    
    def exponentialFrequencyRamp(self,initialFrequency, finalFrequency, dropTime, pulseDuration, maxPhaseError = np.pi/2, amplitude = 1):
        self.S +=f_exponentialFrequencyRamp(initialFrequency, finalFrequency, dropTime, pulseDuration, maxPhaseError, amplitude)
       
    def vfg_pulse(self,old_pulse_dict, new_pulse_dict):
        self.S +=f_vfg_pulse(old_pulse_dict, new_pulse_dict)
    
    def first_vfg_pulse(self,pulse_dict):
        self.S +=f_first_vfg_pulse(pulse_dict)
    
    def pulse_seq(self,pulse_list):
        self.S +=f_pulse_seq(pulse_list)
            
            
  
        