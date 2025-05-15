from labscript import Device, LabscriptError, set_passed_properties
from labscript import config

class matcam(Device):
    description = 'interface for matcam'
        
    def __init__(self, name, parent_device, connection, trigger_edge_type='rising', **kwargs):   
        Device.__init__(self, name, parent_device, connection, **kwargs)       
        self.BLACS_connection = 'Matcam' # appears in BLACS tab
        self.msg = ''
        self.timer = []
        
    def generate_code(self, hdf5_file):
        # Code wich is executed by runmanager to create/srote data in h5 file.
        if self.msg != '':
            # create group with the name of device 
            grp = hdf5_file.create_group('/devices/'+self.name)
            # store the message
            grp.create_dataset('msg', data = self.msg)
        # store timer inf. if not empty
        if self.timer != []:
            grp.create_dataset('tArmMatcam', data = self.timer)
    
    
    # what is the default msg?
    # do not forget the '\n' at the end, otherwise matlab will not read the message!
    
    # expects tuples as input
    def arm(self, *args):
        msg = "||SetLastRunDeleteable false"
        for a in args:
            msg+="||"+a[0]+" "+str(a[1])
        self.msg = msg+"\n"
    
    # expects dictionary as input
    def arm_dict(self, t, dic):
        # add timer information
        self.timer.append(t)
        
        # add message information
        msg = "||SetLastRunDeleteable false"
        for k in dic:
            msg+="||"+k+" "+str(dic[k])
        self.msg = msg+"\n"


    # def save_time(self, timer):
    #     self.timer.append(timer)