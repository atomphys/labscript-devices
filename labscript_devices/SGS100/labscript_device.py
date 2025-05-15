import pyvisa
import numpy as np
from labscript import Device, LabscriptError, set_passed_properties
from labscript import config    
class SGS100(Device):
    """A labscript_device for SGS100 using a visa interface.
    
    define all the functions used by SGS100.
    """
    description = 'SGS100 via Visa'

    @set_passed_properties(property_names = {
        "device_properties":["VISA_name"]}
        )
    
    def __init__(self, name, VISA_name):
        Device.__init__(self, name, None, VISA_name)
        self.BLACS_connection = VISA_name   
        self.preseq = []
        

    def generate_code(self, hdf5_file):
        # Code wich is executed by runmanager to create/store data in h5 file.
        # create group with the name of device        
        grp = hdf5_file.create_group('/devices/'+self.name)
        # Create datasets with all actions to be done before sequence.
        if self.preseq != []:
            grp.create_dataset('preseq',data = self.preseq)


    def normal_operation(self):
        self.preseq.append(np.string_('SOUR:OPM NORM'))

    def reference(self, source):
        self.preseq.append(np.string_('SOUR:ROSC:SOUR '+source))
        
    def synchronisation_bandwidth(self, bandwidth):
        self.preseq.append(np.string_('ROSC:EXT:SBAN '+bandwidth))

    def set_frequency(self, f_int):
        self.preseq.append(np.string_('FREQ '+str(f_int)+' MHz'))

    def set_power(self, p_int):
        self.preseq.append(np.string_('POW '+str(p_int)+' DBM'))
        
    def iq_modulation(self, status):
        self.preseq.append(np.string_('SOUR:IQ:STAT '+status))
        
    def output(self, status):
        self.preseq.append(np.string_('OUTP:STAT '+status))


    def write_this(self, text):
        self.preseq.append(np.string_(text))







        