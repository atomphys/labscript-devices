import pyvisa
from blacs.tab_base_classes import Worker


class data_loggerWorker(Worker):
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
        
        #  check what we are doing
        status = self.dev.query('STATus:OPERation:CONDition?')
        status = int(status.rstrip('\n')) # just a number
        # specific values in manual
        if status == 48:
            return initial_values
            print('We are scanning with no USB device') ### TESTING
        elif status == 176 or status == 432:
            return initial_values
            print('We are scanning and saving to USB device') ### TESTING
        elif status == 128 or status == 384 or status == 0:
            print('We are not scanning, starting remotly!')
            self.dev.write('INIT')
        else:
            print('Not sure what we are doing!')
        
        # read h5 file and execute write commands to device
        with h5py.File(h5file, 'r') as hdf5_file:
            group = hdf5_file['/devices/'+self.device_name]
            # check if there is something to do:
            # else: control manually
            if 'preseq' in list(group.keys()):
                pass ### TODO currently, programming outside of labscript
                # data = [entry.decode() for entry in group['preseq']]
                # for i in range(0, len(data)):
                #     self.write(data[i])

                    
        # return something such that labscript can continue
        return initial_values
 
    def error_check(self):        
        raw_error = self.dev.query('SYST:ERR?')
        error_parts = raw_error.split(',')
        error_message = error_parts[1].rstrip('\n')
        if error_message == '"No error"':
            return True
        else:
            while error_message != '"No error"':
                error_code = int(error_parts[0])
                print(error_code, error_message)
                raw_error = self.dev.query('SYST:ERR?')
                error_parts = raw_error.split(',')
                error_message = error_parts[1].rstrip('\n')
            return True
                
        # TIME
        # dev.query('SYST:TIME?')
        
    def transition_to_manual(self, abort = False):
        if abort:
            # If we're aborting the run, reset to original value
            self.program_manual(True)
            return True
        
        with h5py.File(self.h5_file, 'r+') as self.hdf5_file:
            self.group = self.hdf5_file['/devices/'+self.device_name]
            # check if we want to acquire:
            if 'postseq' in list(self.group.keys()):
                status = self.error_check()
                if status == True:
                    # if so, acquire!
                    channel_number = [entry.decode() for entry in self.group['postseq']]
                    # read out channel_number of last datapoints, and delete them from reading memory
                    # data = self.dev.query('DATA:REMove? '+str(channel_number[0]))
                    # faster version: binblock format
                    # also R? does not return any error if we have nothing in reading memory
                    # theoretically this works without the optional channel_number input
                    # and will just read everything stored in reading memory....
                    self.write('FORMat:READing:TIME ON')
                    self.write('FORMat:READing:TIME:TYPE ABS')
                    self.write('FORMat:READing:UNIT OFF')
                    data = self.query('R?')  #+str(channel_number[0]))
                    data = data[2+int(data[1]):]
                    # if data was collected: save it
                    if len(data)>1:                    
                        ### TESTING: split the collected data into 3 parts, channel, date, timestamp
                        data = data.rsplit(',')
                        chan_number = int(len(data))/7 # 7 entries per channel
                        
                        # post processing
                        self.channel = []
                        self.date = []
                        self.timestamp = []
                        for i in range(0,int(chan_number)):
                            self.channel.append(data[i*7])
                            date = '.'.join(data[i*7+1:i*7+4])
                            self.date.append(date)
                            ### INFO
                            # what format do we want? change accordingly in dataset
                            timestamp = ''.join(data[i*7+4:i*7+7])
                            self.timestamp.append(timestamp)

                        grp = self.hdf5_file.create_group('/data/'+self.device_name)
                        # grp.create_dataset('data', data=np.string_(data))
                        # if we have many channels it may be of interest
                        # to check which ones were scanned and also store this data..
                        
                        ### INFO
                        # here, we can choose how we want to store the collected data
                        # do we want to convert the timestamps to numbers??
                        
                        grp.create_dataset('channels', data=np.array([float(i) for i in self.channel]))
                        grp.create_dataset('date', data=np.string_(self.date))
                        # grp.create_dataset('timestamp', data=np.string_(self.timestamp))
                        grp.create_dataset('timestamp', data=np.array([float(i) for i in self.timestamp]))
                        print('Saving collected data... Done')
                    else:
                        print('seems we did not have any values to store...\n'
                              'is the trigger connected, are we scanning?')
            else:
                print('nothing saved!\n'
                      'use the log_data() function!')
            
                            
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
    










