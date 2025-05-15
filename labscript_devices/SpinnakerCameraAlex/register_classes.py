#####################################################################
#                                                                   #
# /labscript_devices/SpinnakerCameraAlex/register_classes.py        #
#                                                                   #
# Copyright 2019, Monash University and contributors                #
#                                                                   #
# This file is part of labscript_devices, in the labscript suite    #
# (see http://labscriptsuite.org), and is licensed under the        #
# Simplified BSD License. See the license.txt file in the root of   #
# the project for the full license.                                 #
#                                                                   #
#####################################################################




# This is a copy of /labscript_devices/SpinnakerCamera/register_classes.py,
# where we substituted SpinnakerCamera by SpinnakerCameraAlex 
# (to have unique names in the "general" labscript_devices-file (which is necessary!))
# and used BLACS_tab='labscript_devices.SpinnakerCameraAlex.blacs_tabs.SpinnakerCameraTab'
# instead of BLACS_tab='labscript_devices.SpinnakerCamera.blacs_tabs.SpinnakerCameraTab'.




from labscript_devices import register_classes

register_classes(
    'SpinnakerCameraAlex',
    BLACS_tab='labscript_devices.SpinnakerCameraAlex.blacs_tabs.SpinnakerCameraTab',
    runviewer_parser=None,
)
