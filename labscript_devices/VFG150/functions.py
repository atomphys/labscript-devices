import numpy as np
import math

from labscript import LabscriptError

from python_toolbox import caching

# VFG class for caching calculations
# this will only create one instance of VFG()
# all calculations will be memorised/cached
# (only by methods with the @caching.cache() decorator)

### NOTE ###
# calling a cached result will not execute the function and self.S
# is actually not updated....
# use local variable in labscript_file and save to labscript_device directly
# only use these functions to return the resultls!
#
# other option: use a save function to updated self.S after calling a cached result...
# this would require to reset self.S every run, since no new instance is created...
# otherwise this would scale to infinity..

class VFG_cached_class(metaclass=caching.CachedType):
    def __init__(self):
        self.S = []

    # not used!
    # def reset_S(self):
    #     self.S = []
    
    def freq(self, frequency):
        ''' freq(F) Generates a VFG-150 sequence that sets the frequency to F
            (given in Hz). '''
        ticks = round(frequency * (2**32) / 200e6) # external clock
        ret = [ (7*16 + 4), (math.floor(ticks / 2**24))%(2**8),\
            (math.floor(ticks / 2**16))%(2**8),\
            (math.floor(ticks / 2**8))%(2**8),\
            (ticks)%(2**8) ]
       
        # self.S += ret
        return ret
    
    
    def ampl(self, amplitude):
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
        
        # self.S += ret
        return ret
    
    def ampl_slope(self, slope):
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
            
        # self.S += ret
        return ret
    
    
    # One could fix the problem with inverted outputs with few lines of: if port ==1: port = 4 etc...
    def aux(self, port, state):
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
        
        # self.S += ret
        return ret
    
    def interval(self, time):
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
        
        # self.S += ret
        return ret
    
    def phase(self, phi):
        '''phase(phi) Generates a VFG-150 sequence that sets the phase to PHI. PHI
            is given in radians.
            The notion of phase with this command is that of an absolute
            offset phase. So the phase jump that the generated waveform experiences
            is given by the difference between the current phase value and the new
            phase value that is set with this command. The default value for the
            offset phase at the beginning of each sequence is zero. '''
           
        ticks = round(phi * 2**16 / 2 / np.pi)
        ret = [ (9*16 + 8), (math.floor(ticks/2**8))%(2**8), (ticks)%(2**8) ]
        
        # self.S += ret
        return ret
    
    def reset_timebase(self):
        ''' reset_timebase() Generates a VFG-150 sequence that resets the time base
            register. This is the register which is incremented at each time step to
            track the absolute phase of the generated signal.'''
            
        ret = [ 4 ]
        
        # self.S += ret
        return ret
    
    def set_phase_continuous(self):
        ''' set_phase_continuous() generates a VFG-150 sequence that sets the
            synthesizer to phase continuous frequency switching as opposed to
            phase coherent switching. The default at the beginning of each
            sequence is phase continuous switching. '''
        
        ret = [ 2 ]
        
        # self.S += ret
        return ret
    
    
    def set_phase_coherent(self):
        ''' set_phase_coherent() generates a VFG-150 sequence that sets the
            synthesizer to phase coherent frequency switching as opposed to
            phase continous switching. The default at the beginning of each
            sequence is phase continuous switching. '''
        
        ret = [ 3 ]
        
        # self.S += ret
        return ret
    
    
    def trigger(self):
        ''' trigger() generates a VFG-150 sequence that inserts a "wait for trigger"
            event at the corresponding position in the sequence. A "wait for
            trigger" events stops the processing of the sequence until a trigger
            is registered at the trigger input of the VFG-150. The trigger input
            is found at the front panel and labeled "Trigger In".'''
            
        ret = [48, 0, 0, 0, 0]
        
        # self.S += ret
        return ret
    
    def newInterval(self, t0, ph0, a0, a1, a2, err):
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
    

    def generalRamp(self, phaseFunc, phaseFuncD1, phaseFuncD2, amplitudeFunc,\
        pulseDuration, maxPhaseError):
        
        if pulseDuration <= 0:
            raise LabscriptError('Pulse duration must be positive');
        
        if maxPhaseError <= 0:
            raise LabscriptError('Maximum phase error must be positive')
    
        # build the list of times and accumulated phases
        L = np.array([[0,0]])
        while L[-1,0] < pulseDuration:
            add = self.newInterval(L[-1,0],L[-1,1],\
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
            sequence = sequence + self.freq(frequencies[i]) +\
                self.ampl(startamplitudes[i]) +\
                self.ampl_slope(amplituderates[i]) +\
                self.interval(intervals[i])
    
        # at the end the sequence finishes on a constant note
        sequence = sequence + self.freq(phaseFuncD1(pulseDuration)/(2*np.pi)) +\
            self.ampl(amplitudeFunc(pulseDuration)) +\
            self.ampl_slope(0)
            
        sequenceVec = np.concatenate((sequenceVec, np.array([[phaseFuncD1(pulseDuration)/(2*np.pi),\
            amplitudeFunc(pulseDuration), 0, np.inf]])))
    
        return sequence, sequenceVec
    

    def generalRamp_no_ampl_change(self, phaseFunc, phaseFuncD1, phaseFuncD2, amplitudeFunc,\
        pulseDuration, maxPhaseError):
        # to be used for frequency ramp with const ampl (more efficient)
        
        if pulseDuration <= 0:
            raise LabscriptError('Pulse duration must be positive');
        
        if maxPhaseError <= 0:
            raise LabscriptError('Maximum phase error must be positive')
    
        # build the list of times and accumulated phases
        L = np.array([[0,0]])
        while L[-1,0] < pulseDuration:
            add = self.newInterval(L[-1,0],L[-1,1],\
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
        sequence += self.ampl_slope(0) + self.ampl(amplitudeFunc(0))
        
        for i in range(len(intervals)):
            sequence += self.freq(frequencies[i]) +\
                        self.interval(intervals[i])
    
        # at the end the sequence finishes on a constant note
        sequence += self.freq(phaseFuncD1(pulseDuration)/(2*np.pi)) +\
            self.ampl(amplitudeFunc(pulseDuration)) +\
            self.ampl_slope(0)
        
        # for debugging
        sequenceVec = np.c_[frequencies, intervals]
        sequenceVec = np.concatenate((sequenceVec, np.array([[phaseFuncD1(pulseDuration)/(2*np.pi), 0]])))
    
        return sequence, sequenceVec
    
    
    @caching.cache()
    def exponentialFrequencyRamp(self, initialFrequency, finalFrequency,\
        dropTime, pulseDuration, maxPhaseError = np.pi/2, amplitude = 1):
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
        amplitudeFunc = lambda t: amplitude
        
        sequence, sequenceVec = self.generalRamp_no_ampl_change(phaseFunc, phaseFuncD1,\
            phaseFuncD2, amplitudeFunc,\
            pulseDuration, maxPhaseError)

        # sesequenceVec only used for debugging
        # remove decimals
        sequence = [int(round(i)) for i in sequence]        

        return sequence
    


    def vfg_pulse(self, old_pulse_dict, new_pulse_dict):
        S = []
        for key in old_pulse_dict:
            if old_pulse_dict[key] != new_pulse_dict[key]:
                if key == 'aux1':
                    S += self.aux(1, new_pulse_dict[key])
                elif key == 'aux2':
                    S += self.aux(2, new_pulse_dict[key])                
                elif key == 'aux3':
                    S += self.aux(3, new_pulse_dict[key])
                elif key == 'aux4':
                    S += self.aux(4, new_pulse_dict[key])   
                elif key == 'vfg_phase':
                    S += self.phase(new_pulse_dict[key])
                elif key == 'vfg_frequency':
                    S += self.freq(new_pulse_dict[key])
                elif key == 'vfg_amplitude':
                    S += self.ampl(new_pulse_dict[key])
                elif key == 'vfg_slope':
                    S += self.ampl_slope(new_pulse_dict[key])
        # always an interval
        S += self.interval(new_pulse_dict['dt'])
        
        return S
    
    def first_vfg_pulse(self, pulse_dict):
        S = []
        S += self.aux(1, pulse_dict['aux1'])
        S += self.aux(2, pulse_dict['aux2'])                
        S += self.aux(3, pulse_dict['aux3'])
        S += self.aux(4, pulse_dict['aux4'])   
        S += self.phase(pulse_dict['vfg_phase'])
        S += self.freq(pulse_dict['vfg_frequency'])
        S += self.ampl(pulse_dict['vfg_amplitude'])
        S += self.ampl_slope(pulse_dict['vfg_slope'])
        S += self.interval(pulse_dict['dt'])
        return S


    # used to generate the sequence S
    # to be loaded to the VFG
    # input is a list of pulse_dicts, created from pulse class in bec_utils
    def pulse_seq(self, pulse_list):
        # always send
        S = self.set_phase_coherent() + self.trigger()
        
        # create first pulse
        S += self.first_vfg_pulse(pulse_list[0])
        
        # update and create further pulses
        if len(pulse_list) > 1:
            for i in range(len(pulse_list)-1) :
                # replace None with last entry: only works if first pulse has everything defined
                for key in pulse_list[i+1]:
                    if pulse_list[i+1][key] == None:
                        pulse_list[i+1][key] = pulse_list[i][key]
                
                # create new pulses
                S += self.vfg_pulse(pulse_list[i], pulse_list[i+1])
        return S
                    
    
    
    
    

    
    
    
    
    
    
    
    
    
    
# =============================================================================
# 
#     ### TODO
#     # general pulse function for VFG
#     # do we need an additional interval? Wait before doing something...
#     def vfg_pulse(self, oldaux1, oldaux2, oldaux3, oldaux4, oldvfg_phase, oldvfg_slope, oldvfg_amplitude,
#                   oldvfg_frequency, newaux1, newaux2, newaux3, newaux4, newvfg_phase, newvfg_slope, newvfg_amplitude,
#                   newvfg_frequency, dt, **kwargs):
#         
#         S = []
#         
#         if oldaux1 != newaux1:
#             S += self.aux(1,newaux1)
#         if oldaux2 != newaux2:
#             S += self.aux(2,newaux2)
#         if oldaux3 != newaux3:
#             S += self.aux(3,newaux3)
#         if oldaux4 != newaux4:
#             S += self.aux(4,newaux4)
#             
#         if oldvfg_phase != newvfg_phase:
#             S += self.phase(newvfg_phase)
#             
#         if oldvfg_frequency != newvfg_frequency:
#             S += self.freq(newvfg_frequency)
#             
#         if oldvfg_amplitude != newvfg_amplitude:
#             S += self.ampl(newvfg_amplitude)
#     
#         if oldvfg_slope != newvfg_slope:
#             S += self.ampl_slope(newvfg_slope)
#         
#         S += self.interval(dt)
#     
#         return S
#     
#     ### TODO
#     def vfg_seq(self, **kwargs):
#         # always send
#         S = self.set_phase_continuous() + self.trigger()
#     
#         for i in range(len(kwargs['aux1'])):
#             if i == 0:
#                 # always set inital values (+1 to old)
#                 S += self.vfg_pulse(oldaux1 = kwargs['aux1'][i]+1,
#                                oldaux2 = kwargs['aux2'][i]+1,
#                                oldaux3 = kwargs['aux3'][i]+1,
#                                oldaux4 = kwargs['aux4'][i]+1,
#                                oldvfg_phase = kwargs['vfg_phase'][i]+1,
#                                oldvfg_slope = kwargs['vfg_slope'][i]+1,
#                                oldvfg_amplitude = kwargs['vfg_amplitude'][i]+1, 
#                                oldvfg_frequency = kwargs['vfg_frequency'][i]+1, 
#                                newaux1 = kwargs['aux1'][i],
#                                newaux2 = kwargs['aux2'][i],
#                                newaux3 = kwargs['aux3'][i],
#                                newaux4 = kwargs['aux4'][i],
#                                newvfg_phase = kwargs['vfg_phase'][i],
#                                newvfg_slope = kwargs['vfg_slope'][i],
#                                newvfg_amplitude = kwargs['vfg_amplitude'][i], 
#                                newvfg_frequency = kwargs['vfg_frequency'][i], 
#                                dt = kwargs['dt'][i])
#             else:
#                 S += self.vfg_pulse(oldaux1 = kwargs['aux1'][i-1],
#                                oldaux2 = kwargs['aux2'][i-1],
#                                oldaux3 = kwargs['aux3'][i-1],
#                                oldaux4 = kwargs['aux4'][i-1],
#                                oldvfg_phase = kwargs['vfg_phase'][i-1],
#                                oldvfg_slope = kwargs['vfg_slope'][i-1],
#                                oldvfg_amplitude = kwargs['vfg_amplitude'][i-1], 
#                                oldvfg_frequency = kwargs['vfg_frequency'][i-1], 
#                                newaux1 = kwargs['aux1'][i],
#                                newaux2 = kwargs['aux2'][i],
#                                newaux3 = kwargs['aux3'][i],
#                                newaux4 = kwargs['aux4'][i],
#                                newvfg_phase = kwargs['vfg_phase'][i],
#                                newvfg_slope = kwargs['vfg_slope'][i],
#                                newvfg_amplitude = kwargs['vfg_amplitude'][i], 
#                                newvfg_frequency = kwargs['vfg_frequency'][i], 
#                                dt = kwargs['dt'][i])    
#         return S
#     
#     
#     
# =============================================================================












