#####################################################################
#                                                                   #
# /labscript_devices/SpinnakerCameraAlex/blacs_workers.py           #
#                                                                   #
# Copyright 2019, Monash University and contributors                #
#                                                                   #
# This file is part of labscript_devices, in the labscript suite    #
# (see http://labscriptsuite.org), and is licensed under the        #
# Simplified BSD License. See the license.txt file in the root of   #
# the project for the full license.                                 #
#                                                                   #
#####################################################################




# The first part of this file is a copy of the class Spinnaker_Camera from
# /labscript_devices/SpinnakerCamera/blacs_workers.py,
# the second part of this file is a copy of the class IMAQdxCameraWorker from 
# /labscript_devices/IMAQdxCamera/blacs_workers.py,
# where we renamed the class IMAQdxCameraWorker by SpinnakerCameraWorker and used 
# interface_class = Spinnaker_Camera instead of interface_class = IMAQdx_Camera.
# Of course, in the beginning we also imported all importations stated in the "original" files.
# Additionally we imported Image and ImageChops from PIL and did the following things in the class SpinnakerCameraWorker:
# 1) Adding self.calculations = None, self.n_calculated_images = None and self.image_calculations = [] in the init-function. 
# 2) Editing the get_camera-function. 
# 3) Adding new code towards the end of the transition_to_manual-function. [most relevant part!]

# Moreover in the class SPinnakerCameraObject in the __init__-function we set 
# self.timeout = 50000 [ms] instead of 5000.


# Original imaqdx_camera server by dt, with modifications by rpanderson and cbillington.
# Refactored as a BLACS worker by cbillington
# PsSpin implementation by spe



import sys
from time import perf_counter
from blacs.tab_base_classes import Worker
import threading
import numpy as np
from labscript_utils import dedent
import labscript_utils.h5_lock
import h5py
import labscript_utils.properties
import zmq   

from labscript_utils.ls_zprocess import Context
from labscript_utils.shared_drive import path_to_local
from labscript_utils.properties import set_attributes

from enum import IntEnum
from time import sleep, perf_counter

from PIL import Image, ImageChops  # added (for image editing)


# Don't import nv yet so as not to throw an error, allow worker to run as a dummy
# device, or for subclasses to import this module to inherit classes without requiring
# nivision
nv = None



class Spinnaker_Camera(object):
    def __init__(self, serial_number):
        """Initialize Spinnaker API camera.

        Serial number should be of string(?) type."""
        global PySpin
        import PySpin
        
        self.system = PySpin.System.GetInstance()

        ver = self.system.GetLibraryVersion()
        min_ver = (1,23,0,27) # first release with python 3.6 support
        if (ver.major, ver.minor, ver.type, ver.build) < min_ver:
            raise RuntimeError(f"PySpin version {ver} must be >= {min_ver}")

        camList = self.system.GetCameras()
        numCams = camList.GetSize()

        if numCams==0:
            raise ValueError('No cameras found!')

        if isinstance(serial_number, int):
            self.camera = camList.GetBySerial('%d' % serial_number)
        else:
            self.camera = camList.GetBySerial(serial_number)
        self.camera.Init()
        camList.Clear()

        # Set the timeout to 50 s:
        self.timeout = 50000 # in ms

        # Set the abort acquisition thingy:
        self._abort_acquisition = False
        self.exception_on_failed_shot = True

    def get_attribute_names(self, visibility):
        names = []
        def get_node_names_in_category(node_category, prefix=''):
            for node_feature in node_category.GetFeatures():
                # Ensure node is available and readable
                if (not PySpin.IsAvailable(node_feature) or not
                    PySpin.IsReadable(node_feature)):
                    continue

                # Get the feature name:
                feature_name = node_feature.GetName()

                # Category nodes must be dealt with separately in order to retrieve subnodes recursively.
                if node_feature.GetPrincipalInterfaceType() == PySpin.intfICategory:
                    get_node_names_in_category(PySpin.CCategoryPtr(node_feature),
                                               prefix=feature_name + '::')
                else:
                    names.append(prefix + feature_name)

        node = self.camera.GetNodeMap()
        get_node_names_in_category(PySpin.CCategoryPtr(node.GetNode('Root')))

        return names

    def get_attribute(self, name, stream_map=False):
        """Return current values dictionary of attribute of the given name"""
        #print('Getting attribute %s.' % name)
        name = name.split('::')

        if stream_map:
            nodemap = self.camera.GetTLStreamNodeMap()
        else:
            nodemap = self.camera.GetNodeMap()
        node = nodemap.GetNode(name[-1])

        if PySpin.IsAvailable(node) and PySpin.IsReadable(node):
            if node.GetPrincipalInterfaceType() == PySpin.intfIInteger:
                return PySpin.CIntegerPtr(node).GetValue()
            elif node.GetPrincipalInterfaceType() == PySpin.intfIFloat:
                return PySpin.CFloatPtr(node).GetValue()
            elif node.GetPrincipalInterfaceType() == PySpin.intfIBoolean:
                return PySpin.CBooleanPtr(node).GetValue()
            else:
                return PySpin.CValuePtr(node).ToString()

    def set_attributes(self, attr_dict):
        for k, v in attr_dict.items():
            self.set_attribute(k, v)

    def set_stream_attribute(self, name, value):
        self.set_attribute(name, value, stream_map=True)

    def set_attribute(self, name, value, stream_map=False):
        #print('Setting attribute %s.' % name)
        name = name.split('::')

        if stream_map:
            nodemap = self.camera.GetTLStreamNodeMap()
        else:
            nodemap = self.camera.GetNodeMap()
        node = nodemap.GetNode(name[-1])

        if PySpin.IsAvailable(node) and PySpin.IsWritable(node):
            if node.GetPrincipalInterfaceType() == PySpin.intfIInteger:
                 PySpin.CIntegerPtr(node).SetValue(value)
            elif node.GetPrincipalInterfaceType() == PySpin.intfIFloat:
                 PySpin.CFloatPtr(node).SetValue(value)
            elif node.GetPrincipalInterfaceType() == PySpin.intfIBoolean:
                PySpin.CBooleanPtr(node).SetValue(value)
            else:
                PySpin.CValuePtr(node).FromString(value)

            sleep(0.05)
            # Sometimes this doesn't work, so let's check and print warnings if it
            # fails:
            name = '::'.join(name)
            return_value = self.get_attribute(name, stream_map=stream_map)
            if return_value != value:
                print('WARNING: setting attribute %s to %s failed. '%(name, str(value)) +
                      'Returned value %s instead'%str(return_value))
            else:
                print('Successfully set %s to %s.'%(name, str(return_value)))
        else:
            print('WARNING: not capable of writing attribute %s.'%'::'.join(name))


    def snap(self):
        """Acquire a single image and return it"""
        self.configure_acquisition(continuous=False, bufferCount=1)
        #self.trigger()
        image = self.grab()
        self.camera.EndAcquisition()
        return image

    def grab(self):
        """Grab and return single image during pre-configured acquisition."""
        #print('Grabbing...')
        image_result = self.camera.GetNextImage(self.timeout)
        img = self._decode_image_data(image_result.GetData())
        image_result.Release()
        return img

    def grab_multiple(self, n_images, images):
        """Grab n_images into images array during buffered acquistion."""
        print(f"Attempting to grab {n_images} images.")
        for i in range(n_images):
            if self._abort_acquisition:
                print("Abort during acquisition.")
                self._abort_acquisition = False
                return

            images.append(self.grab())
            print(f"Got image {i+1} of {n_images}.")
        print(f"Got {len(images)} of {n_images} images.")

    def trigger(self):
        """Execute software trigger"""
        nodemap = self.camera.GetNodeMap()
        trigger_cmd = PySpin.CCommandPtr(nodemap.GetNode('TriggerSoftware'))
        if not PySpin.IsAvailable(trigger_cmd) or not PySpin.IsWritable(trigger_cmd):
            print('WARNING: Unable to execute trigger. Aborting...')
        else:
            trigger_cmd.Execute()

    def configure_acquisition(self, continuous=True, bufferCount=10):
        self.pix_fmt = self.get_attribute('PixelFormat')
        self.height = self.get_attribute('Height')
        self.width = self.get_attribute('Width')

        # Unless the camera settings are set properly, in cntinuous mode
        # the camera will generally move faster than BLACS, and so the buffer
        # will fill up.  With a Flea3, I was unable to solve the prolem
        # easily.  It really is quite annoying.
        if continuous:
            self.set_stream_attribute('StreamBufferCountMode', 'Manual')
            self.set_stream_attribute('StreamBufferCountManual', bufferCount)
            self.set_stream_attribute('StreamBufferHandlingMode', 'NewestFirst')
            self.set_attribute('AcquisitionMode', 'Continuous')
        elif bufferCount == 1:
            # The StreamBufferCountMode originally was set to 'Auto', but this feature was depreciated by Spinnaker version 3.0.0.118
            self.set_stream_attribute('StreamBufferCountMode', 'Manual')
            self.set_stream_attribute('StreamBufferCountManual', 1)
            self.set_stream_attribute('StreamBufferHandlingMode', 'OldestFirst')
            self.set_attribute('AcquisitionMode', 'SingleFrame')
        else:
            self.set_stream_attribute('StreamBufferCountMode', 'Manual')
            self.set_stream_attribute('StreamBufferCountManual', bufferCount)
            self.set_stream_attribute('StreamBufferHandlingMode', 'OldestFirst')
            self.set_attribute('AcquisitionMode', 'MultiFrame')
            self.set_attribute('AcquisitionFrameCount', bufferCount)

        self.camera.BeginAcquisition()

    def _decode_image_data(self, img):
        """Spinnaker image buffers require significant formatting.
        This returns what one would expect from a camera.
        configure_acquisition must be called first to set image format parameters."""
        if self.pix_fmt.startswith('Mono'):
            if self.pix_fmt.endswith('8'):
                dtype = 'uint8'
            else:
                dtype = 'uint16'
            image = np.frombuffer(img, dtype=dtype).reshape(self.height, self.width)
        else:
            msg = """Only MONO image types currently supported.
            To add other image types, add conversion logic from returned
            uint8 data to desired format in _decode_image_data() method."""
            raise ValueError(dedent(msg))
        return image.copy()

    def stop_acquisition(self):
        print('Stopping acquisition...')
        self.camera.EndAcquisition()

        # This is supposed to provide debugging info, but as with most things
        # in PySpin, it appears to be completely useless:.
        num_frames=self.get_attribute('StreamTotalBufferCount', stream_map=True)
        failed_frames=self.get_attribute('StreamFailedBufferCount', stream_map=True)
        underrun_frames=self.get_attribute('StreamBufferUnderrunCount', stream_map=True)
        print('Stream info: %s frames acquired, %s failed, %s underrun' %
              (str(num_frames), str(failed_frames), str(underrun_frames)))

    def abort_acquisition(self):
        print('Stopping acquisition...')
        self._abort_acquisition = True

    def close(self):
        print('Closing down the camera...')
        self.camera.DeInit()
        self.camList.Clear()
        self.system.ReleaseInstance()
        
    

class SpinnakerCameraWorker(Worker):
    # Subclasses may override this if their interface class takes only the serial number
    # as an instantiation argument, otherwise they may reimplement get_camera():
    interface_class = Spinnaker_Camera
   
    def init(self):
        self.camera = self.get_camera()
        print("Setting attributes... (This statement comes from SpinnakerCameraAlex-blacs_workers)")
        self.smart_cache = {}
        self.set_attributes_smart(self.camera_attributes)
        self.set_attributes_smart(self.manual_mode_camera_attributes)
        print("Initialisation complete")
        self.images = None
        self.n_images = None
        self.attributes_to_save = None
        self.exposures = None
        self.calculations = None  # added [but seems to be not necessary]
        self.n_calculated_images = None  # added [but seems to be not necessary]    
        self.image_calculations = []  # added [but not seems to be not necessary]
        self.acquisition_thread = None
        self.h5_filepath = None
        self.stop_acquisition_timeout = None
        self.exception_on_failed_shot = None
        self.continuous_stop = threading.Event()
        self.continuous_thread = None
        self.continuous_dt = None
        self.image_socket = Context().socket(zmq.REQ)
        self.image_socket.connect(
            f'tcp://{self.parent_host}:{self.image_receiver_port}'
        )

    def get_camera(self):
        """Return an instance of the camera interface class. Subclasses may override
        this method to pass required arguments to their class if they require more
        than just the serial number."""
        # New: Take MockCamera not into account anymore (other option: keep everything as it was 
        # and copy class MockCamera from /labscript_devices/IMAQdxCamera/blacs_workers.py) 
        # [this is probably not that relevant because it also worked without editing]
        # if self.mock:
        #     return MockCamera()
        # else:
        #     return self.interface_class(self.serial_number)
        return self.interface_class(self.serial_number) 

    def set_attributes_smart(self, attributes):
        """Call self.camera.set_attributes() to set the given attributes, only setting
        those that differ from their value in, or are absent from self.smart_cache.
        Update self.smart_cache with the newly-set values"""
        uncached_attributes = {}
        for name, value in attributes.items():
            if name not in self.smart_cache or self.smart_cache[name] != value:
                uncached_attributes[name] = value
                self.smart_cache[name] = value
        self.camera.set_attributes(uncached_attributes)

    def get_attributes_as_dict(self, visibility_level):
        """Return a dict of the attributes of the camera for the given visibility
        level"""
        names = self.camera.get_attribute_names(visibility_level)
        attributes_dict = {name: self.camera.get_attribute(name) for name in names}
        return attributes_dict

    def get_attributes_as_text(self, visibility_level):
        """Return a string representation of the attributes of the camera for
        the given visibility level"""
        attrs = self.get_attributes_as_dict(visibility_level)
        # Format it nicely:
        lines = [f'    {repr(key)}: {repr(value)},' for key, value in attrs.items()]
        dict_repr = '\n'.join(['{'] + lines + ['}'])
        return self.device_name + '_camera_attributes = ' + dict_repr

    def snap(self):
        """Acquire one frame in manual mode. Send it to the parent via
        self.image_socket. Wait for a response from the parent."""
        image = self.camera.snap()
        self._send_image_to_parent(image)

    def _send_image_to_parent(self, image):
        """Send the image to the GUI to display. This will block if the parent process
        is lagging behind in displaying frames, in order to avoid a backlog."""
        metadata = dict(dtype=str(image.dtype), shape=image.shape)
        self.image_socket.send_json(metadata, zmq.SNDMORE)
        self.image_socket.send(image, copy=False)
        response = self.image_socket.recv()
        assert response == b'ok', response

    def continuous_loop(self, dt):
        """Acquire continuously in a loop, with minimum repetition interval dt"""
        while True:
            if dt is not None:
                t = perf_counter()
            image = self.camera.grab()
            self._send_image_to_parent(image)
            if dt is None:
                timeout = 0
            else:
                timeout = t + dt - perf_counter()
            if self.continuous_stop.wait(timeout):
                self.continuous_stop.clear()
                break

    def start_continuous(self, dt):
        """Begin continuous acquisition in a thread with minimum repetition interval
        dt"""
        assert self.continuous_thread is None
        self.camera.configure_acquisition()
        self.continuous_thread = threading.Thread(
            target=self.continuous_loop, args=(dt,), daemon=True
        )
        self.continuous_thread.start()
        self.continuous_dt = dt
        
    def stop_continuous(self, pause=False):
        """Stop the continuous acquisition thread"""
        assert self.continuous_thread is not None
        self.continuous_stop.set()
        self.continuous_thread.join()
        self.continuous_thread = None
        self.camera.stop_acquisition()
        # If we're just 'pausing', then do not clear self.continuous_dt. That way
        # continuous acquisition can be resumed with the same interval by calling
        # start(self.continuous_dt), without having to get the interval from the parent
        # again, and the fact that self.continuous_dt is not None can be used to infer
        # that continuous acquisiton is paused and should be resumed after a buffered
        # run is complete:
        if not pause:
            self.continuous_dt = None

    def transition_to_buffered(self, device_name, h5_filepath, initial_values, fresh):
        if getattr(self, 'is_remote', False):
            h5_filepath = path_to_local(h5_filepath)
        if self.continuous_thread is not None:
            # Pause continuous acquistion during transition_to_buffered:
            self.stop_continuous(pause=True)
        with h5py.File(h5_filepath, 'r') as f:
            group = f['devices'][self.device_name]
            if not 'EXPOSURES' in group:
                return {}
            self.h5_filepath = h5_filepath
            self.exposures = group['EXPOSURES'][:]
            self.n_images = len(self.exposures)

            # Get the camera_attributes from the device_properties
            properties = labscript_utils.properties.get(
                f, self.device_name, 'device_properties'
            )        
            camera_attributes = properties['camera_attributes']
            self.stop_acquisition_timeout = properties['stop_acquisition_timeout']
            self.exception_on_failed_shot = properties['exception_on_failed_shot']
            saved_attr_level = properties['saved_attribute_visibility_level']
            self.camera.exception_on_failed_shot = self.exception_on_failed_shot
        # Only reprogram attributes that differ from those last programmed in, or all of
        # them if a fresh reprogramming was requested:
        if fresh:
            self.smart_cache = {}
        self.set_attributes_smart(camera_attributes)
        # Get the camera attributes, so that we can save them to the H5 file:
        if saved_attr_level is not None:
            self.attributes_to_save = self.get_attributes_as_dict(saved_attr_level)
        else:
            self.attributes_to_save = None
        print(f"Configuring camera for {self.n_images} images.")
        self.camera.configure_acquisition(continuous=False, bufferCount=self.n_images)
        self.images = []
        self.acquisition_thread = threading.Thread(
            target=self.camera.grab_multiple,
            args=(self.n_images, self.images),
            daemon=True,
        )
        self.acquisition_thread.start()
        return {}

    def transition_to_manual(self):
        if self.h5_filepath is None:
            print('No camera exposures in this shot.\n')
            return True
        assert self.acquisition_thread is not None
        self.acquisition_thread.join(timeout=self.stop_acquisition_timeout)
        if self.acquisition_thread.is_alive():
            msg = """Acquisition thread did not finish. Likely did not acquire expected
                number of images. Check triggering is connected/configured correctly"""
            if self.exception_on_failed_shot:
                self.abort()
                raise RuntimeError(dedent(msg))
            else:
                self.camera.abort_acquisition()
                self.acquisition_thread.join()
                print(dedent(msg), file=sys.stderr)
        self.acquisition_thread = None

        print("Stopping acquisition.")
        self.camera.stop_acquisition()

        print(f"Saving {len(self.images)}/{len(self.exposures)} images.")

        with h5py.File(self.h5_filepath, 'r+') as f:
            # Use orientation for image path, device_name if orientation unspecified
            if self.orientation is not None:
                image_path = 'images/' + self.orientation
            else:
                image_path = 'images/' + self.device_name
            image_group = f.require_group(image_path)
            image_group.attrs['camera'] = self.device_name
            
            # Save camera attributes to the HDF5 file:
            if self.attributes_to_save is not None:
                set_attributes(image_group, self.attributes_to_save)

            # Whether we failed to get all the expected exposures:
            image_group.attrs['failed_shot'] = len(self.images) != len(self.exposures)

            # key the images by name and frametype. Allow for the case of there being
            # multiple images with the same name and frametype. In this case we will
            # save an array of images in a single dataset.
            images = {
                (exposure['name'], exposure['frametype']): []
                for exposure in self.exposures
            }

            # Iterate over expected exposures, sorted by acquisition time, to match them
            # up with the acquired images:
            self.exposures.sort(order='t')
            for image, exposure in zip(self.images, self.exposures):
                images[(exposure['name'], exposure['frametype'])].append(image)

            # Save images to the HDF5 file:
            for (name, frametype), imagelist in images.items():
                data = imagelist[0] if len(imagelist) == 1 else np.array(imagelist)
                print(f"Saving frame(s) {name}/{frametype}.")
                group = image_group.require_group(name)
                dset = group.create_dataset(
                    frametype, data=data, dtype='uint16', compression='gzip'
                )
                # Specify this dataset should be viewed as an image
                dset.attrs['CLASS'] = np.string_('IMAGE')
                dset.attrs['IMAGE_VERSION'] = np.string_('1.2')
                dset.attrs['IMAGE_SUBCLASS'] = np.string_('IMAGE_GRAYSCALE')
                dset.attrs['IMAGE_WHITE_IS_ZERO'] = np.uint8(0)

        # If the images are all the same shape, send them to the GUI for display:
        try:
            image_block = np.stack(self.images)
        except ValueError:
            print("Cannot display images in the GUI, they are not all the same shape")
        else:
            self._send_image_to_parent(image_block)

        self.images = None
        self.n_images = None
        self.attributes_to_save = None
        self.exposures = None
        self.stop_acquisition_timeout = None
        self.exception_on_failed_shot = None
        print("Setting manual mode camera attributes.")
        self.set_attributes_smart(self.manual_mode_camera_attributes)
        if self.continuous_dt is not None:
            # If continuous manual mode acquisition was in progress before the bufferd
            # run, resume it:
            self.start_continuous(self.continuous_dt)
       
        # NEW CODE ______________________________________________________________________________________________________________________________________________________________________        
        with h5py.File(self.h5_filepath, 'r+') as f:
            # Don't do anything more if the calculate function is not called in the experiment file; 
            # if the calculate function is called at least once get information about parameters and how often it is called
            group = f['devices'][self.device_name]
            if not 'CALCULATIONS' in group:
                 self.h5_filepath = None
                 return True           
            self.calculations = group['CALCULATIONS'][:]
            self.n_calculated_images = len(self.calculations)
          
            # Use orientation for image path, device_name if orientation unspecified and create place for saving calculated images and N 
            if self.orientation is not None:
                image_path = 'images/' + self.orientation
            else:
                image_path = 'images/' + self.device_name
            image_group = f.require_group(image_path)
            calculated_image_group = image_group.require_group("calculated_images")
            
            # Key the calculated images by the stated parameters
            calculated_images = {
                (calculation['calculated_image_name'], calculation['image_name1'], calculation['frametype1'], calculation['image_name2'], calculation['frametype2'],
                 calculation['box_left'], calculation['box_top'], calculation['box_right'], calculation['box_bottom'])
                for calculation in self.calculations
            }
                   
            # Get calculated images and calculate N           
            print(f"Attempting to get {self.n_calculated_images} calculated image(s).")         
            self.image_calculations = []  # needed at this position!            
            i = 0           
            for (calculated_image_name, image_name1, frametype1, image_name2, frametype2, box_left, box_top, box_right, box_bottom) in calculated_images:                 
                # Create arrays of desired images according to given names and frametypes
                #array1 = f["/images/" + self.device_name + "/" + image_name1.decode('utf-8') + "/" + frametype1.decode('utf-8')][:]
                array1 = f["/images/" + self.device_name + "/" + image_name1 + "/" + frametype1][:]
                #array2 = f["/images/" + self.device_name + "/" + image_name2.decode('utf-8') + "/" + frametype2.decode('utf-8')][:]
                array2 = f["/images/" + self.device_name + "/" + image_name2 + "/" + frametype2][:]
                
                # # OPTION 1 for calculating new image from "image1" and "image2" 
                # # Calculate new image from "image1" and "image2"
                # calculatedArray = array1 - array2   # in general: calculatedArray = f(array1,array2)                                                 
                # # Create calculated image from calculatedArray (to be able to cut it) 
                # calculatedImage = Image.fromarray(calculatedArray.astype('uint16'))                            
                                
                # OPTION 2 for calculating new image from "image1" and "image2"  
                # (this is (in comparison to just subtracting to uint-arrays) the better option to see the "real" difference between two images!)
                # Create the images from the arrays and show them [when using uint16 then the images won't be shown correctly by the
                #                                                  .show()-command and also ImageChops.subtract doesn't work!]
                image1 = Image.fromarray(array1.astype('uint8'))
                image2 = Image.fromarray(array2.astype('uint8'))
                #image1.show()
                #image2.show() 
                # Calculate difference between image1 and image2
                calculatedImage = ImageChops.subtract(image1,image2)
                
                # Cut calculated image according to the given box-parameters
                box = [box_left, box_top, box_right, box_bottom]
                calculatedImageCut = calculatedImage.crop(box)     
                print(f"Got calculated image {i+1} of {self.n_calculated_images}.")
                  
                # Show the cropped calculated image
                # calculatedImageCut.show()                
                # print(f"Showed calculated image {i+1} of {self.n_calculated_images}.")
                
                # Create array from the cropped calculated image (for calculating N) 
                # (larger uint-type than before because we want to sum all entries later!)
                calculatedArrayCut = np.asarray(calculatedImageCut, dtype='uint64')             
                                          
                # Sum all entries of calculatedArrayCut
                N = sum(sum(calculatedArrayCut))                        
                #print(calculated_image_name.decode('utf-8'), "_N = ", N)
                print(calculated_image_name+"_N =", N)
                
                # Save calculated image to the HDF5 file
                print(f"Saving calculated image {calculated_image_name}.")
                dset = calculated_image_group.create_dataset(
                    name = calculated_image_name, data=calculatedArrayCut, dtype='uint16', compression='gzip'
                )
                # Specify this dataset should be viewed as an image
                dset.attrs['CLASS'] = np.string_('IMAGE')
                dset.attrs['IMAGE_VERSION'] = np.string_('1.2')
                dset.attrs['IMAGE_SUBCLASS'] = np.string_('IMAGE_GRAYSCALE')
                dset.attrs['IMAGE_WHITE_IS_ZERO'] = np.uint16(0)   
                dset.attrs['atom_number'] = np.uint64(N)
                
                # # Save single N-value to the HDF5 file
                # dset2_name = calculated_image_name.decode('utf-8')+"_N"
                # dset2 = calculated_image_group.create_dataset(
                #     name=dset2_name, data=N, dtype='uint64'
                # )
                
                # This is needed for being able to save all N-values in an array (see below)
                self.image_calculations.append((calculated_image_name, N))  
                                
                i += 1
                   
            # Save all N-values in an array with corresponding calculated_image_name to the HDF5 file
            vlenstr = h5py.special_dtype(vlen=str)
            table_dtypes_N = [
                ('calculated_image_name', vlenstr),
                ('atom_number', 'uint64'),
            ]
            data_N = np.array(self.image_calculations, dtype=table_dtypes_N)
            dset_N = calculated_image_group.create_dataset(
                name='calculated_images_N', data=data_N, dtype=table_dtypes_N
            )  
            
                   
        print() # just to have an empty line in the blacs terminal output       
        self.calculations = None
        self.n_calculated_images = None
        # END NEW CODE _____________________________________________________________________________________________________________________________________________________________
        
        self.h5_filepath = None
        return True

    def abort(self):
        if self.acquisition_thread is not None:
            self.camera.abort_acquisition()
            self.acquisition_thread.join()
            self.acquisition_thread = None
            self.camera.stop_acquisition()
        self.camera._abort_acquisition = False
        self.images = None
        self.n_images = None
        self.attributes_to_save = None
        self.exposures = None
        self.acquisition_thread = None
        self.h5_filepath = None
        self.stop_acquisition_timeout = None
        self.exception_on_failed_shot = None
        # Resume continuous acquisition, if any:
        if self.continuous_dt is not None and self.continuous_thread is None:
            self.start_continuous(self.continuous_dt)
        return True

    def abort_buffered(self):
        return self.abort()

    def abort_transition_to_buffered(self):
        return self.abort()

    def program_manual(self, values):
        return {}

    def shutdown(self):
        if self.continuous_thread is not None:
            self.stop_continuous()
        self.camera.close()
        
