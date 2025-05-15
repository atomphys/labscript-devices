import pyvisa
from blacs.tab_base_classes import Worker


class Keysight_ScopeWorker(Worker):
    def init(self, timeout=5):
        global properties; import labscript_utils.properties
        global h5_lock,h5py; import labscript_utils.h5_lock, h5py
        global np; import numpy as np

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
        self.h5_file = h5file
        self.device_name = device_name
        # read h5 file and execute write commands to device
        with h5py.File(h5file, 'r') as hdf5_file:
            group = hdf5_file['/devices/'+self.device_name]
            # check if there is something to do:
            # else: control manually
            if 'preseq' in list(group.keys()):
                # if we have something to do, set the triggermode to normal and single mode
                # for some reason, KEYSIGHT scopes don't like running mode for saving data
                self.write(':TRIG:SWE NORM')
                data = [entry.decode() for entry in group['preseq']]
                for i in range(0, len(data)):
                    self.write(data[i])
                if 'postseq' in list(group.keys()):
                    self.write(':SING')
                # if we are not acqiering the data, we can also be in running mode
                else:
                    self.write(':RUN')
                    
        # return something such that labscript can continue
        return initial_values
 

    # read out the selected traces
    # only works for single acquisation mode
    def acquire(self):
        traces = []
        preambles = []
        data = [entry.decode() for entry in self.group['postseq']]
        try:
            for i in range(0, len(data)):
                self.dev.write(':WAV:FORM ASC')
                self.dev.write(':WAV:POINTS MAX')
                self.dev.write(':WAV:SOUR CHAN'+data[i])
                preamble = self.dev.query(':WAV:PRE?')
                trace = self.dev.query(':WAV:DATA?')
                # formatting of string: second number indicates how many follow...
                trace = trace[2+int(trace[1]):] # could also be deleted during analysis..
                traces.append(np.string_(trace))
                preambles.append(np.string_(preamble))
        
            grp = self.hdf5_file.create_group('/data/'+self.device_name)
            print('Saving traces...')
            grp.create_dataset('traces', data=traces)
            grp.create_dataset('preambles', data=traces)
            print('Done')

        except:
            print('could not save traces, are the correct channels selected?\n'
                  'can only be saved if displaed on the scope!')
            
        if 'preseq' not in list(self.group.keys()):
            # prepare for next run: single mode
            self.write(':SING')
                
    def transition_to_manual(self,abort = False):
        if abort:
            # If we're aborting the run, reset to original value
            self.program_manual(True)
            return True
        
        with h5py.File(self.h5_file, 'r+') as self.hdf5_file:
            self.group = self.hdf5_file['/devices/'+self.device_name]
            # check if we want to acquire:
            if 'postseq' in list(self.group.keys()):
                # if so, do it!
                self.acquire()
            else:
                print('nothing saved!')

        return True
    
    def program_manual(self,front_panel_values):
        return None
    
    def abort_transition_to_buffered(self):
        """Special abort shot configuration code belongs here.
        """
        return self.transition_to_manual(True)

    def abort_buffered(self):
        """Special abort shot code belongs here.
        """
        return self.transition_to_manual(True)
    
    def shutdown(self):
        """Closes VISA connection to device."""
        self.dev.close()
    










