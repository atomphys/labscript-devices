from blacs.tab_base_classes import Worker
class SGS100Worker(Worker):
    def init(self, timeout = 2):
        global pyvisa; import pyvisa
        global properties; import labscript_utils.properties
        global h5_lock, h5py; import labscript_utils.h5_lock, h5py
        
        # create VISA connection
        rm = pyvisa.ResourceManager()
        self.dev = rm.open_resource(self.VISA_name)
        self.dev.timeout = 1000 * timeout
        # define actions
        self.idn = self.dev.query('*IDN?')
        self.read = self.dev.read
        self.write = self.dev.write
        self.query = self.dev.query
        
        # check which device is connected
        manufacturer, model, sn, revision = self.idn.split(',')
        print('Connected to {} from {} (SN: {})'.format(model, manufacturer, sn))


                
    def transition_to_buffered(self, device_name, h5file, initial_values, fresh):
        # Store some parameters for saving data later
        self.h5_file = h5file
        self.device_name = device_name
        
        # read h5 file and execute write commands to device
        with h5py.File(h5file, 'r') as hdf5_file:
            group = hdf5_file['/devices/'+self.device_name]
            # check if there is something to do:
            if list(group.keys()) != []:
                data = [entry.decode() for entry in group['preseq']]
                for i in range(0, len(data)):
                    self.write(data[i])

        # return something such that labscript can continue
        return initial_values
 
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
    










