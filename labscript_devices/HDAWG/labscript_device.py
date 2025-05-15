
import math
import numpy as np
from labscript import Device, LabscriptError, set_passed_properties
from labscript import config    
import logging
from logging.handlers import RotatingFileHandler
class HDAWG(Device):
    """A labscript_device for the HDAWG such that it can be controlled via python.
    """
    description = 'HDAWG via python'
    def __init__(self, name):
        Device.__init__(self, name, None, None)
        self.name = name
        self.BLACS_connection = 'HDAWG' # This is needed to create the blacs_tab and will appear as [conn: self.BLACS_connection]
        self.IQLoad = []
        self.dim = []
        self.log = logging.getLogger("HDAWG")
        if len(self.log.handlers) == 0:
            format = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
            stream_hndl = logging.StreamHandler()
            stream_hndl.setLevel(logging.INFO)
            stream_hndl.setFormatter(format)
            file_hndl = RotatingFileHandler("hdawg_pulses.log",mode='a',maxBytes=3*1024*1024,backupCount=2, encoding=None, delay=False)
            file_hndl.setLevel(logging.DEBUG)
            file_hndl.setFormatter(format)
            self.log.addHandler(stream_hndl)
            self.log.addHandler(file_hndl)
        self.log.setLevel(logging.DEBUG)
    def generate_code(self, hdf5_file):
        # Code wich is executed by runmanager to create/srote date in h5 file.
        if self.IQLoad != []:
            # create group with the name of device and its data
            grp = hdf5_file.create_group('/devices/'+self.name)
            grp.create_dataset('IQLoad', data = self.IQLoad)

    # takes list of pulses and loads IQSeq to h5 file
    def ks_pulse(self, pulse_dict):
        pulse_program = np.zeros(14)
        pulse_program[0] = pulse_dict['ks_frequency1']
        pulse_program[1] = pulse_dict['ks_frequency2']
        pulse_program[2] = pulse_dict['ks_amplitude1']
        pulse_program[3] = pulse_dict['ks_amplitude2']
        pulse_program[4] = pulse_dict['ks_phase1']
        pulse_program[5] = pulse_dict['ks_phase2']
        pulse_program[6] = pulse_dict['dt']
        pulse_program[7] = pulse_dict['aux1']
        pulse_program[8] = pulse_dict['aux2']
        pulse_program[9] = pulse_dict['aux3']
        pulse_program[10] = pulse_dict['aux4']
        pulse_program[11] = pulse_dict['RF_amp']
        pulse_program[12] = pulse_dict['RF_phase']
        pulse_program[13] = pulse_dict['RF_freq']
        #pulse_program[14] = pulse_dict['descr']
        return np.array(pulse_program)
    
    # takes list of pulses and loads IQSeq to h5 file
    def pulse_seq(self, pulse_list):
        # start with empyt list
        self.log.debug("entering pulse_seq")
        IQSeq = []
        # append the init pulse (first in sequence)
        IQSeq.append(self.ks_pulse(pulse_list[0]))
        # if we have more to do, check what to do
        if len(pulse_list) > 1:
            for i in range(0, len(pulse_list)-1):
                # for a similar pulse, add the time to last entry of IQSeq
                if pulse_list[i]['ks_frequency1'] == pulse_list[i+1]['ks_frequency1'] and\
                    pulse_list[i]['ks_frequency2'] == pulse_list[i+1]['ks_frequency2'] and\
                    pulse_list[i]['ks_amplitude1'] == pulse_list[i+1]['ks_amplitude1'] and\
                    pulse_list[i]['ks_amplitude2'] == pulse_list[i+1]['ks_amplitude2'] and\
                    pulse_list[i]['ks_phase1'] == pulse_list[i+1]['ks_phase1'] and\
                    pulse_list[i]['ks_phase1'] == pulse_list[i+1]['aux1'] and\
                    pulse_list[i]['ks_phase1'] == pulse_list[i+1]['aux2'] and\
                    pulse_list[i]['ks_phase1'] == pulse_list[i+1]['aux3'] and\
                    pulse_list[i]['ks_phase1'] == pulse_list[i+1]['aux4'] and\
                    pulse_list[i]['ks_phase1'] == pulse_list[i+1]['vfg_amplitude'] and\
                    pulse_list[i]['ks_phase1'] == pulse_list[i+1]['vfg_phase'] and\
                    pulse_list[i]['ks_phase1'] == pulse_list[i+1]['vfg_frequency'] and\
                    pulse_list[i]['ks_phase2'] == pulse_list[i+1]['ks_phase2']:
                    
                    IQSeq[-1][6] += pulse_list[i+1]['dt']
                # for a new pulse: append to IQSeq
                else:
                    IQSeq.append(self.ks_pulse(pulse_list[i+1]))
                    self.log.debug(pulse_list[i+1])
        
        # finally, load IQSeq to h5-file
        return self.IQLoadAndPlay(IQSeq)

    def IQLoadAndPlay(self, IQSeq):
        self.IQLoad.append(IQSeq)
        # self.dim.append(len(pulse_program))
            
  
        