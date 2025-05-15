from blacs.tab_base_classes import Worker
import datetime

class matcamWorker(Worker):
    global properties; import labscript_utils.properties
    global h5_lock, h5py; import labscript_utils.h5_lock, h5py
    global socket; import socket
    global system; from os import system
    global sleep; from time import sleep
    #TODO import AlliedVision.universal # this needs to be installed
    global dirname, basename, isfile, join; from os.path import dirname, basename, isfile, join    
    global threading; import threading
    
    def init(self):
        # Start matcam
        # system('matlab -r "matcamRec;" -nosplash -nodesktop')
        
        # Start BLACS -> matcam communication
        self.socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

        # print(self.socket)
        # ip = '192.168.0.13' #the old matcam computor
        ip = '127.0.0.1' #point to itself
        
        port = 4012
        
        # Wait for matcam's responsiveness. 20s is enough when lucky, but can get much slower....
        success = False
        tries = 100
        while not(success): # Matcam takes a while to start, try to connect to it until successful.
            try:
                self.socket.connect((ip,port))
                success = True
            except:
                sleep(1)
                tries-=1
                print(str(tries)+' tries left')
                if tries == 0:
                    raise Exception('Matcam socket failed to connect.')
        
        # otherwise we create errors...
        self.socket.setblocking(0)
        
        ### Testing for plotting
        self.last_sequence_index = []
        self.last_scan_date = ''
        
        
        # TODO: No idea what this does/if needed
        # Setup camera-trigger properly since matcam modified it at startup.
        # if socket.gethostname() == 'phys-atom-23':
            # AlliedVision.universal.setTriggerMode()
        


    def arm(self, msg='', h5file=None):
        """ Sends arm message to matcam via sockets. """
        folderpath = dirname(h5file).split(':')[1]
        filename = basename(h5file).split('.')[0]
        
        # here we can decide where to store the shot! maybe change this part for custom implementation?
        with h5py.File(self.h5file, 'r') as hdf5_file:
            scan = hdf5_file.attrs['run number']
            
            # repetition = hdf5_file.attrs['run repeat'] # only saved if repeat is turned on
            # check if we want to send SetCurrentBaseDir
            if self.last_sequence_index != hdf5_file.attrs['sequence_index'] or self.last_scan_date != hdf5_file.attrs['sequence_date']:
                self.last_sequence_index = hdf5_file.attrs['sequence_index']
                self.last_scan_date = hdf5_file.attrs['sequence_date']
                self.msg = '||SetCurrentBaseDir '+folderpath+'\ ||SetNextFilename '+filename #+str(msg)
            
            # if we are in the same run, only send SetNextFilename
            else:
                self.msg = '||SetNextFilename '+filename #+str(msg)
        armtime=datetime.datetime.now().strftime("%d-%b-%Y-%H-%M-%S-%f")
        self.msg=self.msg+'||armtime '+str(armtime)+str(str(msg.decode()[:-1])+'EOD\n')
        #self.msg=self.msg+'EOD'                  
        print(self.msg)
        
        # encode used to send binary format: b'msg'
        self.socket.send(self.msg.encode())
       
    
    def transition_to_buffered(self, device_name, h5file, initial_values, fresh):
        self.execthreads=[]
        self.h5file = h5file
        msg = ''
        # Read h5 file and prepare the threads.
        with h5py.File(self.h5file, 'r') as hdf5_file:
            group = hdf5_file['/devices/'+device_name]
            if 'msg' in list(group.keys()):
                data = group['msg']
                msg = data[()]
            if 'tArmMatcam' in list(group.keys()):
                timer = group['tArmMatcam']
                tArmMatcam = int(timer[()])
            
                print('arming in '+str(tArmMatcam) +' seconds!')
                thread = threading.Timer(tArmMatcam, self.arm, args = (msg, self.h5file))
                thread.start()
            else:
                print('tArmMatcam was not saved to h5 file {}, will not be armed!'.format(self.h5file))

        return {}
    
    
    # write status to recovery file 
    def testfun(self, s='random text'):
        f = open('C:\\Users\\BEC_control\\Desktop\\MatcamRecovery.txt','a')
        f.write(s+'\n')
        return f.close()

    def transition_to_manual(self):
        # we dont use the recovery!
        return True
        
        sleep(1) # matcam needs some additional delay to save data when repeating?
        # check wether matcam has saved an image correctly
        folderpath = dirname(self.h5file)
        # change this accordingly to arm() definition
        imgfile = join(folderpath,"matcamdata_"+basename(self.h5file)+".mat")
        if not(isfile(imgfile)): # an image has not been saved correctly, try to recover
            tries = 50
            recovered = False
            self.testfun('Image not yet stored correctly...')
            while((tries > 0) and not(recovered)):
                tries-=1
                self.socket.send(b'status?')
                self.testfun('status on %s? sent'%imgfile)
                sleep(2)
                try:
                    ret = self.socket.recv(1000) # matcam should return its status
                except:
                    ret='Unknown'
                self.testfun('Matcam status:'+ret+' -- Tries remaining: %i -- %s'%(tries,imgfile))
                if isfile(imgfile): # after some time, the image is saved correctly
                    self.testfun('Finally found the image file %s.\n ### SUCCESS ###'%imgfile)
                    recovered = True
                    return True
                if 'Active' in ret:
                    self.testfun('Matcam is finally active -- image file might still not exist though. %s\n### FAIL ###'%imgfile)
                    recovered = True
                    return True ### TODO in this case, the last experiment should be repeated since no image was saved?
            if not(recovered):
                print("Error: matcam didn't save the image correctly and automatic recovery failed. Matcam status: %s. %s"%(ret,imgfile))
        return True
       
    def abort_buffered(self):
        return self.abort()
        
    def abort_transition_to_buffered(self):
        return self.abort()
    
    def abort(self):
        return True
    
    def program_manual(self, values):
        return {}

    def shutdown(self):
        self.socket.shutdown()
        self.socket.close()
        return
        









