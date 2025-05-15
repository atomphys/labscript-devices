import pyvisa
import numpy as np
from labscript import Device, LabscriptError, set_passed_properties
from labscript import config
import labscript_utils.h5_lock, h5py

class Keysight_Scope(Device):
    """A labscript_device for Agilent using a visa interface.
    
    define all the functions used by Agilent.
    """
    description = 'Keysight Scope via VISA'

    @set_passed_properties(property_names = {
        "device_properties":["VISA_name"]}
        )
    
    def __init__(self, name, VISA_name):
        Device.__init__(self, name, None, VISA_name)
        self.name = name
        self.BLACS_connection = VISA_name
        self.preseq=[]
        self.postseq=[]

    def generate_code(self, hdf5_file):
        grp = hdf5_file.create_group('/devices/'+self.name)
        # Device.generate_code(hdf5_file)
        # Create datasets with all actions to be done before (preseq) and after sequence.
        if self.preseq != []:
            grp.create_dataset('preseq',data=self.preseq)
        if self.postseq != []:
            grp.create_dataset('postseq',data=self.postseq)

    # Here follow some example functions
    # very specific, should be changed for individual use

    def channel_setup(self, chan_num, scale, offset, unit):
        # this is the full scale: multipy by 8
        scale = scale*8
        self.preseq.append(np.string_(':CHAN'+str(chan_num)+':RANG '+str(scale)+' '+unit))
        self.preseq.append(np.string_(':CHAN'+str(chan_num)+':OFFS '+str(offset)+' '+unit))

    def timescale(self, timescale, units='mS', delay=None):
        # cound also enter manually with timescale = 1e-9 or smthing like that
        # seems to be no optional input of units
        # maybe there is a more elegant way.. 
        if units == 'mS':
            timescale = timescale/1000
            delay = delay/1000
        if units == 'uS':
            timescale=timescale/1000000
            delay = delay /1000000
        if units == 'nS':
            timescale=timescale/1000000000
            delay = delay /1000000000
        self.preseq.append(np.string_('TIM:RANG '+str(timescale)))
                           
        if delay != None:
            self.preseq.append(np.string_(':TIM:DEL '+str(delay)))
    
    def trigger(self, chan, slope='POS', level=1):
        self.preseq.append(np.string_(':TRIG:SLOP '+slope))
        if type(chan) == int:
            self.preseq.append(np.string_(':TRIG:SOUR CHAN'+str(chan)))
        else:
            self.preseq.append(np.string_(':TRIG:SOUR '+chan))
        self.preseq.append(np.string_(':TRIG:LEV '+str(level)))
    
    # only stat which channel to acquire, do rest in worker process
    def acquire(self, chan_num):
        self.postseq.append(np.string_(str(chan_num)))
            
    







        