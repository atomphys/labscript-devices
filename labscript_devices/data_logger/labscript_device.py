import pyvisa
import numpy as np
from labscript import Device, LabscriptError, set_passed_properties
from labscript import config
import labscript_utils.h5_lock, h5py

class data_logger(Device):
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
        # Create datasets with all actions to be done before (preseq) and after sequence.
        if self.preseq != []:
            grp.create_dataset('preseq',data=self.preseq)
        if self.postseq != []:
            grp.create_dataset('postseq',data=self.postseq)

    # Here follow some example functions
    # very specific, should be changed for individual use
    
    # we read out all the data anyway, just append something :)
    def log_data(self, status = True):
        if status == True:
            self.postseq.append(np.string_('logging'))




        