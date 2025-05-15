import pyvisa
import numpy as np
from labscript import Device, LabscriptError, set_passed_properties
from labscript import config    
class Agilent(Device):
    """A labscript_device for Agilent using a visa interface.
    
    define all the functions used by Agilent.
    """
    description = 'Agilent via Visa'

    @set_passed_properties(property_names = {
        "device_properties":["VISA_name"]}
        )
    
    def __init__(self, name, VISA_name):
        Device.__init__(self, name, None, VISA_name)
        self.BLACS_connection = VISA_name
        self.preseq=[]
        

    def generate_code(self, hdf5_file):
        # Code wich is executed by runmanager to create/store data in h5 file.
        # create group with the name of device
        grp = hdf5_file.create_group('/devices/'+self.name)
        # Create datasets with all actions to be done before (preseq) and after sequence.
        if self.preseq != []: 
            grp.create_dataset('preseq',data=self.preseq)

    
    def setfreqpow(self, f_int, p_int):
        self.preseq.append(np.string_(':FREQ:MODE FIX'))
        self.preseq.append(np.string_(':FREQ '+str(f_int)+' MHZ'))
        self.preseq.append(np.string_(':OUTP OFF'))
        self.preseq.append(np.string_(':POW:ALC:SOUR INT'))
        self.preseq.append(np.string_(':POW '+str(p_int)+' DBM'))
        self.preseq.append(np.string_(':OUTP ON'))
        
        
    def setfreqpowALC(self, f_int, p_int):
        self.preseq.append(np.string_(':FREQ:MODE FIX'))
        self.preseq.append(np.string_(':FREQ '+str(f_int)+' MHZ'))
        self.preseq.append(np.string_(':OUTP OFF'))
        self.preseq.append(np.string_(':POW:ALC:SOUR:EXT:COUP -14.16DB'))
        self.preseq.append(np.string_(':POW:ALC:SOUR DIOD'))
        self.preseq.append(np.string_(':POW '+str(p_int)+' DBM'))
        self.preseq.append(np.string_(':OUTP ON'))
        








        