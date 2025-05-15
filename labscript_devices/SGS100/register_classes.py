
"""
Sets which BLACS_tab belongs to each labscript device.
"""

import labscript_devices

labscript_devices.register_classes(
    'SGS100',
    BLACS_tab='labscript_devices.SGS100.blacs_tab.SGS100Tab',
    runviewer_parser=None
)
