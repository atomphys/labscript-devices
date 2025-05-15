
"""
Sets which BLACS_tab belongs to each labscript device.
"""

import labscript_devices

labscript_devices.register_classes(
    'matcam',
    BLACS_tab='labscript_devices.matcam.blacs_tab.matcamTab',
    runviewer_parser=None
)
