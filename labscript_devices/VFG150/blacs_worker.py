from blacs.tab_base_classes import Worker

class VFG150Worker(Worker):
    
    def init(self):
        global properties; import labscript_utils.properties
        global h5_lock, h5py; import labscript_utils.h5_lock, h5py
        global np; import numpy as np
        global LabscriptError; from labscript import LabscriptError
        global ct; import ctypes as ct
        global cdll; from ctypes import cdll
        global threading; import threading
        global slp; from time import sleep as slp
        
        global ctime; from time import ctime as ctime
        ###
        import usb.core
        import usb.util
        
        # find the vfg-150
        self.vfg = usb.core.find(idVendor=0x0bd7, idProduct=0xa035)
        if self.vfg is None:
            raise ValueError('VFG not found')
        else:
            print('connected to '+ self.vfg.product + ' from '+self.vfg.manufacturer)

        
        
# =============================================================================
#         ### Old version (like Matlab)
#         # define the path of .dll file (save .h file in same location)
#         dll_path = "C:\\Users\\BEC_control\\Desktop\\VFG150_Matlab_original\\usbdrvd.dll"
#         
#         # dll_path = "C:\\Users\\BEC_control\\Desktop\\DLL\\DLL used by Matlab and Goodtime\\usbdrvd.dll"
#         # load the library
#         # self.vfg150 = cdll.LoadLibrary(dll_path)
#         self.vfg150 = ct.WinDLL(dll_path)
#         # define arg types for BulkWrite
#         # needed???
#         self.vfg150.USBDRVD_BulkWrite.argtypes= [ct.c_ulong, ct.c_ulong, ct.c_void_p, ct.c_ulong]
#         
#         # check if a device is connected, else raise error
#         self.dev_count = self.vfg150.USBDRVD_GetDevCount()
#         if self.dev_count == 0:
#             raise LabscriptError('Could not find a suitable usb device!')
#         
# =============================================================================
     
    def new_program_VFG(self):
        with h5py.File(self.h5_file, 'r') as hdf5_file:
            group = hdf5_file['/devices/'+self.device_name]
            # load the data from h5 file
            S = group['S'][()]
            # add an end-of-sequence marker to the sequence (will not be executed otherwise)
            S = np.append(S, [48, 255, 255, 255, 255])
            # send the data
            #print('reset and sending data')
            self.vfg.reset()
            # print('sending data '+str(len(S)))
            buffer = self.vfg.write(2, S, timeout = 5000)  # The endpoint to write to is 2 (see VFG Manual from Toptica page 45), timeout 5s
                
            if buffer != len(S):
                print(ctime()+' : '+str(buffer) + ' out of '+ str(len(S)))
                self.vfg.reset()
                # try to reset
                self.vfg.write(2, [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]+[48, 255, 255, 255, 255], timeout = 5000)
                
            ### Matlab approach: sending in pieces
            # bytestowrite = len(S);
            # sendtotal = 0
            # for i in np.arange(0, bytestowrite, 1024):
            #     write_size = min(1024, bytestowrite)
            #     Send = S[i:i+write_size]
            #     buffer = self.vfg.write(2, Send,timeout = 5000)
            #     bytestowrite = bytestowrite  - write_size
            #     sendtotal += buffer
            #     print('sent '+(str(sendtotal)+'/'+str(len(S)-bytestowrite)))
                
            # buffer = self.vfg.write(2, [48, 255, 255, 255, 255])
            # print('sent '+(str(sendtotal)+'/'+str(len(S))))
                   
            # buffer = self.vfg.write(2, S,5000) # The endpoint to write to is 2 (see VFG Manual from Toptica page 45), timeout 5s
            # self.vfg.reset()
            
            
            
    def program_VFG(self):
        # open the device handle
        self.handle = self.vfg150.USBDRVD_OpenDevice(ct.c_ulong(1))
        if self.handle == 0:
            raise LabscriptError('Could not get a handle for the device')            
            
        # read the data from h5 file and execute it
        print('connected and loading')
        with h5py.File(self.h5_file, 'r') as hdf5_file:
            group = hdf5_file['/devices/'+self.device_name]
            S = group['S'][()]
        
            # add an end-of-sequence marker to the sequence
            S = np.append(S, [48, 255, 255, 255, 255])
            bytestowrite = len(S)
            
            ### resembles Matlab vfg150 approach
            ### sending in pieces
            # print(np.arange(0,bytestowrite,2000))
            # for i in np.arange(0,bytestowrite,2000):
            #     print(i)
            #     write_size = min(2000,bytestowrite)
            #     print(write_size)
            #     c_array = ct.c_int8 * write_size
            #     print(c_array)
            #     Send = c_array(*np.uint8(S[i:i+write_size]))
            #     print(S[i:i+write_size])
            #     result = self.vfg150.USBDRVD_BulkWrite(self.handle, 0, Send, len(Send))
            #     bytestowrite = bytestowrite - write_size
            #########################

           
            # create empty c_int8 array
            c_array = ct.c_int8 * len(S)
            # store data on the array
            Send = c_array(*np.uint8(S))

            # mydata = bytearray(S)
            # buffer_for_S = ct.create_string_buffer(bytes(mydata), len(mydata))
            
            # send data to vfg150, if unsuccessful: raise error
            # result = self.vfg150.USBDRVD_BulkWrite(self.handle, 0, buffer_for_S, len(mydata))
            result = self.vfg150.USBDRVD_BulkWrite(self.handle, 0, Send, len(Send))
            print(str(result) + ' out of ' + str(len(Send)))   
            if result == 0:
                # close device connection
                # self.vfg150.USBDRVD_CloseDevice(1)
                raise LabscriptError('Could not write to the pipe')
        print('loaded and closing')
        slp(1) ### not needed?
        # close device connection  
        self.vfg150.USBDRVD_CloseDevice(1)
        print('closed connection')
        
        
    def transition_to_buffered(self, device_name, h5file, initial_values, fresh):
        self.h5_file = h5file
        self.device_name = device_name

        # do we want to program the VFG in a thread?
        # check if we need the VFG
        with h5py.File(h5file, 'r') as hdf5_file:
            group = hdf5_file['/devices/'+self.device_name]
            if 'S' in list(group.keys()):
                # self.new_program_VFG()
          
                # self.thread = threading.Timer(0.1, self.program_VFG)
                self.thread = threading.Timer(0.1, self.new_program_VFG)

                # self.thread.name='VFGThread'
                self.thread.start()
                # print(threading.enumerate())
                
            else:
                self.thread = None



        # ############## Original ######################################################
        # # TODO: could this be done already in init? never close the device connection?)
        # # open the device handle
        # self.handle = self.vfg150.USBDRVD_OpenDevice(1)
        # print('open')
        # if self.handle == 0:
        #     raise LabscriptError('Could not get a handle for the device')
        #    
        # # read the data from h5 file and execute it
        # with h5py.File(h5file, 'r') as hdf5_file:
        #     group = hdf5_file['/devices/'+self.device_name]
        #     data = group['S']
        #     S = data[()]
        #     # TODO there is a matlab version (currently used?) which uses 49 instead of 48... why?
        #     # add an end-of-sequence marker to the sequence
        #     # S = np.append(S, [48, 255, 255, 255, 255])
        #     # create empty c_int8 array
        #     c_array = ct.c_int8 * len(S)
        #     # store data on the array
        #     S = c_array(*np.uint8(S))
        #     # send data to vfg150, if unsuccessful: raise error
        #     result = self.vfg150.USBDRVD_BulkWrite(self.handle, 0, S, len(S))
        #     print(result)
        #     if result == 0:
        #         # close device connection
        #         self.vfg150.USBDRVD_CloseDevice(1)
        #         raise LabscriptError('Could not write to the pipe')
        # # close device connection
        # print('closing')
        # self.vfg150.USBDRVD_CloseDevice(1)
        # print('closed')
        # ################################################################################
        
        return initial_values

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
            
    def transition_to_manual(self, abort = False):
        """Simple transition_to_manual method where no data is saved."""         
        # If we're not aborting the run, stick with buffered value. Nothing to do really!
        # print('Waiting for thread to finish')
        if self.thread != None:
            self.thread.join()
        # print('Thread finished')
        return True
        
    def shutdown(self):
        print('shutdown')
        """Closes vfg150 handle and frees library"""
        # close device connection
        # self.vfg150.USBDRVD_CloseDevice(1)
        # Unload the DLL so that it can be rebuilt
        libHandle = self.vfg150._handle
        del self.vfg150
        
        # free the library
        self.free_library(libHandle)

    def free_library(self, handle):
        kernel32 = ct.WinDLL('kernel32', use_last_error=True)
        kernel32.FreeLibrary.argtypes = [ct.wintypes.HMODULE]
        kernel32.FreeLibrary(handle)








