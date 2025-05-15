from blacs.tab_base_classes import Worker

class iq_modulatorWorker(Worker):
    def init(self):
        global properties; import labscript_utils.properties
        global h5_lock, h5py; import labscript_utils.h5_lock, h5py
        global np; import numpy as np
        global math; import math
        global pyvisa; import pyvisa
        global threading; import threading
        
        """Opens VISA connection to device."""
        rm = pyvisa.ResourceManager()
        self.iq_modulator = rm.open_resource(self.VISA_name)
        # set a higher default PyVisa timeout so that we can load
        # all the data we wnat onto te iq_modulator
        self.iq_modulator.timeout = 25000
        self.idn = self.iq_modulator.query('*IDN?')
        self.read = self.iq_modulator.read
        self.write = self.iq_modulator.write
        self.query = self.iq_modulator.query
        self.write_ascii_values = self.iq_modulator.write_ascii_values
        self.write_binary_values = self.iq_modulator.write_binary_values
        
        # check which device is connected
        manufacturer, model, sn, revision = self.idn.split(',')
        print('Connected to {} from {} (SN: {})'.format(model, manufacturer, sn))
        #print('Remember to change something in the function... Paolo said someting like: we need to add a math.floor somewhere.. Compare the latest Matlab script with the current version of Python functions :)')
        

    def program_ks(self):
        with h5py.File(self.h5_file, 'r') as hdf5_file:
            group = hdf5_file['/devices/'+self.device_name]
            # read out the dimension of the matrix
            # dim = group['dim'][()]
    
            # generate the structure of pulse_program
            # pulse_program = np.zeros((int(dim),7))
            # load the data from h5 file
            ### read_direct: old version
            # group['IQLoad'].read_direct(pulse_program)
            # print(pulse_program) ### testing
            
            ### testing: just load the data
            pulse_program = group['IQLoad'][()]
            print(pulse_program)
            # send functions to iq_modulator
            
            ### we need better input control
            ### how is it saved to h5 file?
            self.IQPlay(pulse_program.reshape(np.shape(pulse_program)[1], 7))

    def transition_to_buffered(self, device_name, h5file, initial_values, fresh):
        self.h5_file = h5file
        self.device_name = device_name
        with h5py.File(h5file, 'r') as hdf5_file:
            group = hdf5_file['/devices/'+self.device_name]
            # if there is no data assigned to iq_modulator: do nothing
            if 'IQLoad' in list(group.keys()):        
                # start the programming 0.1s into the experiment
                thread = threading.Timer(0.1, self.program_ks)
                thread.start()
            else:
                self.IQModulatorOutputOff()
                print('no IQLoad found for this run')

        return initial_values
                
 
    def IQModulatorSetup(self):
        # Channel 1 should be input into I, channel 2 into Q
        # setattr() used instead of ser(), new: replaced with higher timout
        sample_rate = 250e6
        max_sample_arb = 1e4                             # max. number of samples for a single-tone arb
        self.write('FORM:BORD SWAP')                     # correct memory format
        self.write('FUNC:ARB:SRAT ' + str(sample_rate))  # sample rate
        self.write('DATA:ARB2:FORM ABAB')                # data structure for 2-channel arb wforms
        self.write('FUNC:ARB:PTP 1')                     # scaling between given value and output voltage
        self.write('OUTP:TRIG OFF')                      # make sure the trigger output on the ext trig port (which we need to accept trigger input) is off
        self.write('TRIG:SOUR EXT')                      # receive triggers from the ext trig port
        self.write('TRIG:SLOP POS')                      # trigger on the negative slope (for inverted trigger signals)
        self.write('UNIT:ANGL RAD')                      # set angle units to radians
        return sample_rate, max_sample_arb
    
    def CheckInput(self, input):
        # Make sure that zero-amp. signals get zero freq. (avoid useless comput.)
        input_corr = input
        input_corr[:,0:2] = input[:,0:2]*(1*(input[:,2:4]!=0))
        
        # Make sure that only second tone doesn't get mistaken for two tones
        for i in range(0, len(input_corr[:,0])-1):
            if input_corr[i][0]==0 and input_corr[i][1]!=0:
                 input_corr[i,[0, 1, 2, 3, 4, 5]] = input_corr[i,[1, 0, 3, 2, 5, 4]]
        return input_corr
    
    def IQModulatorOutputOn(self):
        self.write('OUTP1 1')
        self.write('OUTP2 1')
        
    def IQModulatorOutputOff(self):
        self.write('OUTP1 0')
        self.write('OUTP2 0')
    
    def IQModulatorRead(self, command):
        self.write(command)
        message = self.read()
        return message
        
    def VoltToDAC(self, volt):
        dac_val = np.int16(32767 * volt)
        return dac_val
    
    def WaveForm(self, time, freq1, freq2, amp1, amp2, phase1, phase2, imbalance, phase_corr, dc_1, dc_2):
        if np.sqrt(amp1**2 + amp2**2) > 1:
            print('Total IQ amplitude beyond full scale.');
            wform = 0;
        else:
            wform = np.int16(np.zeros(2*np.size(time)))
            wform_I = self.VoltToDAC(imbalance[0] * amp1 * np.cos(2 * np.pi * freq1 * time + phase1 + phase_corr[0]) + imbalance[1] * amp2 * np.cos(2 * np.pi * freq2 * time + phase2 + phase_corr[1]) - dc_1)
            wform_Q = self.VoltToDAC(amp1 * np.sin(2 * np.pi * freq1 * time + phase1) + amp2 * np.sin(2 * np.pi * freq2 * time + phase2) - dc_2)
            wform[0::2] = wform_I
            wform[1::2] = wform_Q
        return wform
    
    def FreqApprox(self, freq, sample_rate, max_sample_arb):
        if freq == 0:
            freq_appr = 0
            sample_number = 8
            n_p = 1
        else:
            freq_sign = np.sign(freq)
            freq_abs = abs(freq)
            n_p_max = math.floor(freq_abs * max_sample_arb / sample_rate)        # max number of periods we can fit into max_sample_arb
            attempts = np.linspace(1, n_p_max, n_p_max)                          # try out all numbers of periods fitting in the maximum number of samples
            att_sample_numbers = np.round(sample_rate*attempts/freq_abs)
            att_freqs_appr = sample_rate * attempts / att_sample_numbers
            errors = abs(att_freqs_appr - freq_abs)/attempts                     # the error must be divided by the number of periods to really get the error per arb repetition
            min_value = np.min(errors)
            indices_min = np.where(errors == min_value)
            array_of_mins = indices_min[0]
            i_best = array_of_mins[0]
            n_p = attempts[i_best]
            sample_number = att_sample_numbers[i_best]
            freq_appr = freq_sign * att_freqs_appr[i_best]
        return freq_appr, sample_number, n_p
    
    def ImpairmentCalibration(self):
        freq_cal = 1e6* np.array([-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15])
        phase_corr_cal = 1e-4*np.array([-107,-99,-96.5,-90.25,-83.5,-76.75,-71,-63.5,-56.25,-49.25,-42.75,-36,-28.5,-21,-12.75,2.9,10.5,16.9,24.8,32.5,39.75,45,53.5,61.5,67.75,74.25,80,85.5,92.5,98.75])
        imbalance_cal = np.array([0.9976,0.9974,0.9974,0.99705,0.996725,0.99665,0.99635,0.9963,0.996225,0.9959,0.9958,0.99565,0.995575,0.995525,0.99535,0.99525,0.99535,0.9954,0.99545,0.995575,0.99565,0.99585,0.9961,0.996225,0.996425,0.996525,0.9968,0.9971,0.99735,0.99745])
    
        p_ph = np.polyfit(freq_cal,phase_corr_cal,3)
        p_imb = np.polyfit(freq_cal,imbalance_cal,4)
        
        return p_ph, p_imb
    
    def ImpairmentCalculate(self, input):
        # these get calculated in ImpairmentCalibration.m
        # is it needed or does it stay the same
        [p_ph, p_imb] = self.ImpairmentCalibration()
        phase_corr = np.polyval(p_ph, input[:,0:2])
        imbalance = np.polyval(p_imb, input[:,0:2])
        # The DC offsets seem to drift randomly, so I've just let them to zero...
        dc_1 = 0e-4 # dc offset to the signal of channel 1
        dc_2 = 0e-4 # same for channel 2
        return phase_corr, imbalance, dc_1, dc_2
    
    
    
    def IQPlay(self, input):
        AllTheTime = False

        # Make sure that the basic setting of the device are correct
        [sample_rate, max_sample_arb] = self.IQModulatorSetup()
        
        # Make sure that the input is formatted in a compatible way
        input = self.CheckInput(input)
        
        # Correct for the constant imperfections in the signal (IQ impairments)
        [phase_corr,imbalance,dc_1,dc_2] = self.ImpairmentCalculate(input)
        
        # Define the arrays that we'll need during the evaluation etc.
        freq_appr = np.zeros(len(input[:,0]))        # best approximation for the frequency
        sample_number = np.zeros(len(input[:,0]))    # best number of samples for each arb
        n_p = np.zeros(len(input[:,0]))              # number of periods executed in every arb
        n_rep = np.zeros(len(input[:,0]))            # number of repetitions of every arb during the sequence
        sample_rest = np.zeros(len(input[:,0]))      # number of samples still to be executed after the repetitions
        phase1 = np.zeros(len(input[:,0]))           # phases taking into account the timing
        phase2 = np.zeros(len(input[:,0]))
        times = np.insert(np.cumsum(input[:,6]),0,0) # time points (instead of time intervals)
    
        # the string that is being sent to the device to create the sequence
        if AllTheTime == True:
            seq_command = 'sequence,wform0,1,once,lowAtStart,4'
        else:
            seq_command = 'sequence,wform0,1,onceWaitTrig,lowAtStart,4,wform0,1,repeat,highAtStart,4' 
    
        # One additional zero waveform is added after trigger, otherwise channel 1
        # somehow doesn't go to zero, but is set to a few mV until trigger...
        # If this is not a problem, one can eliminate that additional waveform.
    
        # Turn off the outputs and clear the device's volatile memory before loading
        self.IQModulatorOutputOff()
        self.write('*CLS')
        self.write('SOUR1:DATA:VOL:CLE')
        self.write('SOUR2:DATA:VOL:CLE')
        
        # Initial zero arb to trick this fucker into working with triggers
        self.write_ascii_values('DATA:ARB2:DAC wform0,', np.zeros(16))
        
        # Go through the input 'vertically': for every iteration create the
        # waveform, load it onto the device, and update the sequence command
        # change the range of i due to matlab vs. python convention
        for i in range(0, len(input[:,0])):
            # The marker (e.g. to trigger an oscilloscope) will be on high only during the first waveform
            if i==0:
                marker = 'highAtStart'
                if AllTheTime == True:
                    repeat = ',repeatTilTrig,'
                else:
                    repeat = ',repeat,'
            # elif (i+1)%2 == 1:
            #     marker = 'highAtStart'
            #     repeat = ',repeat,'
            else:
                marker = 'lowAtStart'
                repeat = ',repeat,'
    
        
            # Distinguish between single- and two-tone parts
            # Single-tone parts can be sequenced with smaller arbs    
            if input[i][1] == 0: # single-tone parts
                # Approximate the frequencies to something suited to the sample rate
                # Since we need to sequence the signal, not having an 'integer number of samples' per arb would induce a phase jump at every repetition,
                # i.e. change the frequency of the signal.
                [freq_appr[i],sample_number[i],n_p[i]] = self.FreqApprox(input[i][0], sample_rate, max_sample_arb)
                print(freq_appr[i]-input[i][0])
                # Generate and upload the periodic arb signal
                # Generate DAC values:
                n_rep[i] = math.floor(sample_rate * input[i][6] / sample_number[i])
                phase1[i] = (input[i][4] - 2*np.pi*(freq_appr[i]*times[i]%(1)))%(2*np.pi)
                wform = self.WaveForm(np.linspace(1,int(sample_number[i]),int(sample_number[i]))/sample_rate, freq_appr[i], 0, input[i][2],0 ,phase1[i], 0, imbalance[i,:], phase_corr[i,:], dc_1, dc_2)
                # Configure and send arb file:
                header = 'DATA:ARB2:DAC wform' + str(i + 1) + ','
                self.write_binary_values(header, wform, datatype = 'h')
                # self.write('*WAI')
                # Update sequence command:
                seq_command = seq_command + ',wform' + str(i + 1) + ',' + str(int(n_rep[i]))  + repeat + marker + ',4'
                # Generate and upload rest arb signal, if necessary, and update the sequence
                sample_rest[i] = np.floor(np.remainder(input[i][6]*sample_rate, sample_number[i]))
                # print(input[i][6], sample_number)
                if sample_rest[i] > 0:
                    # Generate DAC values:
                    wform_rest = self.WaveForm(np.linspace(1,int(sample_number[i]),int(sample_number[i]))/sample_rate, freq_appr[i], 0, input[i][2], 0, phase1[i], 0, imbalance[i,:], phase_corr[i,:], dc_1, dc_2)
                    # Configure and send arb file:
                    header = 'DATA:ARB2:DAC wform' + str(i + 1) + 'rest,'
                    self.write_binary_values(header, wform_rest, datatype = 'h')
                    self.write('*WAI')
                    # Update sequence command:
                    seq_command = seq_command + ',wform' + str(i + 1) + 'rest,' + str(1) + ',repeat,' + marker + ',4'
            else: # two-tone parts
                # Generate and upload the arb signal
                # Generate DAC values:
                phase1[i] = (input[i][4] - 2*np.pi*(input[i][0]*times[i]%1))%(2*np.pi)
                phase2[i] = (input[i][5] - 2*np.pi*(input[i][1]*times[i]%1))%(2*np.pi)
                sample_number[i] = np.round(input[i,6]*sample_rate)
                wform = self.WaveForm(np.linspace(1,int(sample_number[i]),int(sample_number[i]))/sample_rate, input[i][0], input[i][1], input[i][2], input[i][3], phase1[i], phase2[i], imbalance[i,:], phase_corr[i,:], dc_1, dc_2)
                # Configure and send arb file:
                header = 'DATA:ARB2:DAC wform' + str(i + 1) + ','
                self.write_binary_values(header, wform, datatype = 'h')
                self.write('*WAI')
                # Update sequence command
                seq_command = seq_command + ',wform' + str(i + 1) + ',' + str(int(n_rep[i])) + ',repeat,' + marker + ',4'
        
        # Deterimine the lenth of the seq_command to manually do what binblockwrite does automatically
        length = len(seq_command)
        length_of_length = len(str(length))
        # Send sequence command
        self.write_ascii_values('SOUR1:DATA:SEQ #'+str(length_of_length)+str(length), seq_command, converter = 's', separator='')
        self.write('*WAI')
        
        # Play the sequence
        self.write('SOUR1:FUNC:ARB sequence')
        self.write('SOUR1:FUNC ARB')
        
        # Check that no errors occured before turning on the output
        raw_error = self.IQModulatorRead('SYST:ERR?')
        error_parts = raw_error.split(',')
        error_message = error_parts[1].rstrip('\n')
    
        # Turn on the outputs, if no errors occured
        if error_message == '"No error"':
            self.IQModulatorOutputOn()
            print('Loaded and turned on')
        else:
            while error_message != '"No error"':
                error_code = int(error_parts[0])
                print(error_code, error_message)
                raw_error = self.IQModulatorRead('SYST:ERR?')
                error_parts = raw_error.split(',')
                error_message = error_parts[1].rstrip('\n')
                
                
        #  freq_appr      # best approximation for the frequency
        #  sample_number  # best number of samples for each arb
        #  n_p            # number of periods executed in every arb
        #  n_rep          # number of repetitions of every arb during the sequence
        #  sample_rest    # number of samples still to be executed after the repetitions
        #  phase1         # phases taking into account the timing
        #  phase2
        #  times


    def check_remote_values(self):
        return None
    
    def program_manual(self,front_panel_values):
        return self.check_remote_values()
    
    def abort_transition_to_buffered(self):
        """Special abort shot configuration code belongs here.
        """
        return self.transition_to_manual(True)

    def abort_buffered(self):
        """Special abort shot code belongs here.
        """
        return self.transition_to_manual(True)
            
    def transition_to_manual(self,abort = False):
        """Simple transition_to_manual method where no data is saved."""         
        if abort:
            # If we're aborting the run, reset to original value
            self.program_manual(True)
        # If we're not aborting the run, stick with buffered value. Nothing to do really!
        return True
        
    def shutdown(self):
        """Closes VISA connection to device."""
        self.dev.close()










